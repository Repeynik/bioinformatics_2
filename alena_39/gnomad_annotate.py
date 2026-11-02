#!/usr/bin/env python3
import csv
import json
import os
import sys
import time
import urllib.request
from datetime import date

GNOMAD_API = 'https://gnomad.broadinstitute.org/api'
DATASET = 'gnomad_r4'
BATCH = 10
DELAY = 3
FIELDS = '''
    variant_id rsids
    exome { ac an af homozygote_count filters }
    genome { ac an af homozygote_count filters }
    joint {
      ac an homozygote_count filters
      populations { id ac an }
    }
'''
COLUMNS = [
    'id', 'gnomad_variant_id', 'in_gnomad', 'rsids',
    'exome_ac', 'exome_an', 'exome_af', 'exome_filters',
    'genome_ac', 'genome_an', 'genome_af', 'genome_filters',
    'joint_ac', 'joint_an', 'joint_af', 'joint_homozygotes', 'joint_filters',
    'grpmax_af', 'grpmax_group',
]
GRPMAX_GROUPS = ('afr', 'amr', 'eas', 'nfe', 'sas')


def read_vcf(path):
    records = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom, pos, vid, ref, alt = line.rstrip('\n').split('\t')[:5]
            records.append({'chrom': chrom, 'pos': int(pos), 'id': vid, 'ref': ref, 'alt': alt})
    return records


def gnomad_id(rec):
    return f"{rec['chrom']}-{rec['pos']}-{rec['ref']}-{rec['alt']}"


def query_batch(recs):
    parts = [f'v{i}: variant(variantId: "{gnomad_id(r)}", dataset: {DATASET}) {{ {FIELDS} }}'
             for i, r in enumerate(recs)]
    body = json.dumps({'query': '{ ' + ' '.join(parts) + ' }'}).encode()
    req = urllib.request.Request(GNOMAD_API, data=body, headers={'Content-Type': 'application/json'})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                data = json.load(resp)
            return [data.get('data', {}).get(f'v{i}') for i in range(len(recs))]
        except urllib.error.HTTPError as e:
            if e.code != 429:
                raise
            print(f"  rate limit, waiting {10 * (attempt + 1)} s ...")
            time.sleep(10 * (attempt + 1))
    raise RuntimeError('gnomad api: too many retries')


def af(block):
    if not block or not block.get('an'):
        return None
    return block['ac'] / block['an']


def fmt(x):
    if x is None:
        return '-'
    if isinstance(x, float):
        return f'{x:.3g}'
    return str(x)


def variant_row(rec, v):
    row = {'id': rec['id'], 'gnomad_variant_id': gnomad_id(rec), 'in_gnomad': 'yes' if v else 'no'}
    if v:
        ex, ge, jo = v.get('exome'), v.get('genome'), v.get('joint')
        row['rsids'] = ','.join(v.get('rsids') or [])
        for name, block in (('exome', ex), ('genome', ge), ('joint', jo)):
            block = block or {}
            row[f'{name}_ac'] = block.get('ac')
            row[f'{name}_an'] = block.get('an')
            row[f'{name}_af'] = af(block)
            row[f'{name}_filters'] = ','.join(block.get('filters') or []) or ('PASS' if block else None)
        if jo:
            row['joint_homozygotes'] = jo.get('homozygote_count')
            pops = [p for p in jo.get('populations') or [] if p['id'] in GRPMAX_GROUPS and p['an']]
            if pops:
                best = max(pops, key=lambda p: p['ac'] / p['an'])
                row['grpmax_af'] = best['ac'] / best['an']
                row['grpmax_group'] = best['id']
    return {k: fmt(row.get(k)) for k in COLUMNS}


def save_settings(path, n_variants):
    with open(path, 'w', encoding='utf-8') as f:
        f.write("gnomAD annotation settings\n")
        f.write("==========================\n")
        f.write(f"date: {date.today().isoformat()}\n")
        f.write(f"source: gnomAD GraphQL API ({GNOMAD_API})\n")
        f.write(f"dataset: {DATASET} (gnomAD v4, GRCh38, exomes + genomes)\n")
        f.write(f"input variants: {n_variants}\n")
        f.write("\nquery:\n")
        f.write("  variant(variantId: '<chrom>-<pos>-<ref>-<alt>', dataset: gnomad_r4)\n")
        f.write("  fields: exome / genome / joint ac, an, homozygote_count, filters,\n")
        f.write("          joint populations\n")
        f.write("\nderived values:\n")
        f.write("  *_af = ac / an\n")
        f.write("  joint_af = allele frequency over exomes + genomes (used in the final table)\n")
        f.write("  grpmax_af = max af over afr, amr, eas, nfe, sas (gnomAD grpmax groups)\n")
        f.write("  filters: PASS if the variant passed all gnomAD filters\n")


def main():
    if len(sys.argv) < 2:
        print("usage: python gnomad_annotate.py <file.vcf> [out_dir]")
        sys.exit(1)
    out_dir = sys.argv[2] if len(sys.argv) > 2 else 'gnomad'
    os.makedirs(out_dir, exist_ok=True)
    records = read_vcf(sys.argv[1])

    print(f"  querying gnomad ({DATASET}) for {len(records)} variants ...")
    raw = {}
    for i in range(0, len(records), BATCH):
        chunk = records[i:i + BATCH]
        for rec, v in zip(chunk, query_batch(chunk)):
            raw[rec['id']] = v
        time.sleep(DELAY)
    rows = [variant_row(rec, raw[rec['id']]) for rec in records]

    with open(os.path.join(out_dir, 'gnomad_output.json'), 'w', encoding='utf-8') as f:
        json.dump(raw, f, indent=1)
    with open(os.path.join(out_dir, 'gnomad_results.tsv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    save_settings(os.path.join(out_dir, 'gnomad_settings.txt'), len(records))
    print(f"  variants found in gnomad: {sum(r['in_gnomad'] == 'yes' for r in rows)}")
    print(f"  saved to {out_dir}/")
    return 0


if __name__ == '__main__':
    sys.exit(main())
