#!/usr/bin/env python3
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from datetime import date

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
TOOL = 'nsu_bioinf_task4'
DELAY = 0.4
COLUMNS = [
    'id', 'chrom', 'pos', 'ref', 'alt', 'in_clinvar', 'variation_id', 'accession', 'title',
    'classification', 'review_status', 'last_evaluated', 'conditions', 'molecular_consequence',
    'other_records_at_position',
]


def read_vcf(path):
    records = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom, pos, vid, ref, alt = line.rstrip('\n').split('\t')[:5]
            records.append({'chrom': chrom, 'pos': int(pos), 'id': vid, 'ref': ref, 'alt': alt})
    return records


def eutils(endpoint, **params):
    params.update(db='clinvar', retmode='json', tool=TOOL)
    url = f'{EUTILS}/{endpoint}.fcgi?{urllib.parse.urlencode(params)}'
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = json.load(resp)
    time.sleep(DELAY)
    return data


def search_position(chrom, pos):
    term = f'{chrom}[CHR] AND {pos}[CPOS]'
    return eutils('esearch', term=term, retmax=100)['esearchresult']['idlist']


def fetch_summaries(uids, batch_size=50):
    summaries = {}
    for i in range(0, len(uids), batch_size):
        result = eutils('esummary', id=','.join(uids[i:i + batch_size]))['result']
        for uid in result.get('uids', []):
            summaries[uid] = result[uid]
    return summaries


def spdi_of(summary):
    for vs in summary.get('variation_set', []):
        spdi = vs.get('canonical_spdi', '')
        if spdi:
            return spdi.split(':')
    return None


def matches(summary, rec):
    spdi = spdi_of(summary)
    if not spdi or len(spdi) != 4:
        return False
    _, pos0, ref, alt = spdi
    return int(pos0) == rec['pos'] - 1 and ref == rec['ref'] and alt == rec['alt']


def record_row(rec, matched, others):
    row = {k: rec[k] for k in ('id', 'chrom', 'pos', 'ref', 'alt')}
    row['in_clinvar'] = 'yes' if matched else 'no'
    if matched:
        s = matched[0]
        gc = s.get('germline_classification', {})
        row.update({
            'variation_id': s['uid'],
            'accession': s.get('accession_version', ''),
            'title': s.get('title', ''),
            'classification': gc.get('description', ''),
            'review_status': gc.get('review_status', ''),
            'last_evaluated': gc.get('last_evaluated', '').split(' ')[0],
            'conditions': '; '.join(t['trait_name'] for t in gc.get('trait_set', [])),
            'molecular_consequence': ', '.join(s.get('molecular_consequence_list', [])),
        })
    row['other_records_at_position'] = ' | '.join(
        f"{s.get('accession', '')} {s.get('title', '')}: "
        f"{s.get('germline_classification', {}).get('description') or 'no classification'}"
        for s in others)
    for k in COLUMNS:
        row.setdefault(k, '-')
        if row[k] == '':
            row[k] = '-'
    return row


def save_settings(path, n_variants, db_info):
    with open(path, 'w', encoding='utf-8') as f:
        f.write("ClinVar annotation settings\n")
        f.write("===========================\n")
        f.write(f"date: {date.today().isoformat()}\n")
        f.write(f"source: NCBI E-utilities ({EUTILS})\n")
        f.write(f"clinvar last update: {db_info.get('lastupdate', 'N/A')}\n")
        f.write(f"clinvar records: {db_info.get('count', 'N/A')}\n")
        f.write(f"input variants: {n_variants}\n")
        f.write("\nquery:\n")
        f.write("  esearch db=clinvar term='<chrom>[CHR] AND <pos>[CPOS]' (CPOS = GRCh38 position)\n")
        f.write("  esummary db=clinvar id=<found ids>\n")
        f.write("\nmatching:\n")
        f.write("  a record is accepted only if canonical SPDI equals <chrom>:<pos-1>:<REF>:<ALT>\n")
        f.write("  other records at the same position (other alleles, indels, inversions)\n")
        f.write("  are listed in other_records_at_position and are not used for the variant\n")
        f.write("  classification = germline_classification (aggregate VCV level)\n")


def main():
    if len(sys.argv) < 2:
        print("usage: python clinvar_annotate.py <file.vcf> [out_dir]")
        sys.exit(1)
    out_dir = sys.argv[2] if len(sys.argv) > 2 else 'clinvar'
    os.makedirs(out_dir, exist_ok=True)
    records = read_vcf(sys.argv[1])

    print(f"  searching clinvar for {len(records)} positions ...")
    hits = {}
    for rec in records:
        hits[rec['id']] = search_position(rec['chrom'], rec['pos'])
    uids = sorted({u for ids in hits.values() for u in ids})
    print(f"  found {len(uids)} clinvar records at these positions")
    summaries = fetch_summaries(uids)
    db_info = eutils('einfo')['einforesult']['dbinfo'][0]

    rows = []
    raw = {}
    for rec in records:
        found = [summaries[u] for u in hits[rec['id']] if u in summaries]
        matched = [s for s in found if matches(s, rec)]
        others = [s for s in found if not matches(s, rec)]
        rows.append(record_row(rec, matched, others))
        raw[rec['id']] = {'query_hits': hits[rec['id']], 'matched': matched, 'other': others}

    with open(os.path.join(out_dir, 'clinvar_output.json'), 'w', encoding='utf-8') as f:
        json.dump(raw, f, indent=1, ensure_ascii=False)
    with open(os.path.join(out_dir, 'clinvar_results.tsv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    save_settings(os.path.join(out_dir, 'clinvar_settings.txt'), len(records), db_info)
    print(f"  variants with clinvar record: {sum(r['in_clinvar'] == 'yes' for r in rows)}")
    print(f"  saved to {out_dir}/")
    return 0


if __name__ == '__main__':
    sys.exit(main())
