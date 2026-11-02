#!/usr/bin/env python3
import csv
import json
import os
import sys
import time
import urllib.request
from datetime import date

REST = 'https://rest.ensembl.org'
VEP_URL = f'{REST}/vep/homo_sapiens/region'
VEP_PARAMS = {
    'canonical': 1,
    'mane': 1,
    'hgvs': 1,
    'numbers': 1,
    'protein': 1,
    'variant_class': 1,
    'appris': 1,
    'tsl': 1,
    'ccds': 1,
    'uniprot': 1,
    'AlphaMissense': 1,
    'REVEL': 1,
    'CADD': 1,
    'SpliceAI': 1,
}
TSV_COLUMNS = [
    'Uploaded_variation', 'Location', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene',
    'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp',
    'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons',
    'Existing_variation', 'DISTANCE', 'STRAND', 'CANONICAL', 'MANE_SELECT', 'TSL', 'APPRIS',
    'CCDS', 'ENSP', 'SWISSPROT', 'SIFT', 'PolyPhen', 'AlphaMissense', 'REVEL', 'CADD_PHRED',
    'SpliceAI_DS_max', 'gnomADe_AF', 'gnomADg_AF', 'CLIN_SIG', 'PUBMED',
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


def get_json(url):
    req = urllib.request.Request(url, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.load(resp)


def post_json(url, payload):
    body = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=body, headers={
        'Content-Type': 'application/json', 'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=300) as resp:
        return json.load(resp)


def run_vep(records, batch_size=100):
    lines = [f"{r['chrom']} {r['pos']} {r['id']} {r['ref']} {r['alt']} . . ." for r in records]
    results = []
    for i in range(0, len(lines), batch_size):
        payload = dict(VEP_PARAMS, variants=lines[i:i + batch_size])
        results.extend(post_json(VEP_URL, payload))
        time.sleep(1)
    return results


def range_str(start, end):
    if start is None:
        return '-'
    return str(start) if start == end else f'{start}-{end}'


def known_variant(result, allele):
    for cv in result.get('colocated_variants', []):
        if cv.get('id', '').startswith('rs'):
            freqs = cv.get('frequencies', {}).get(allele, {})
            return cv, freqs
    return {}, {}


def allele_clin_sig(cv, allele):
    for part in cv.get('clin_sig_allele', '').split(';'):
        if part.startswith(f'{allele}:'):
            return part.split(':', 1)[1]
    return ''


def transcript_row(result, tc):
    allele = tc.get('variant_allele', '')
    cv, freqs = known_variant(result, allele)
    spliceai = tc.get('spliceai') or {}
    ds = [spliceai[k] for k in ('DS_AG', 'DS_AL', 'DS_DG', 'DS_DL') if k in spliceai]
    am = tc.get('alphamissense') or {}
    sift = f"{tc['sift_prediction']}({tc['sift_score']})" if 'sift_prediction' in tc else '-'
    pp = f"{tc['polyphen_prediction']}({tc['polyphen_score']})" if 'polyphen_prediction' in tc else '-'
    return {
        'Uploaded_variation': result['id'],
        'Location': f"{result['seq_region_name']}:{result['start']}",
        'Allele': allele,
        'Consequence': ','.join(tc.get('consequence_terms', [])),
        'IMPACT': tc.get('impact', '-'),
        'SYMBOL': tc.get('gene_symbol', '-'),
        'Gene': tc.get('gene_id', '-'),
        'Feature_type': 'Transcript',
        'Feature': tc.get('transcript_id', '-'),
        'BIOTYPE': tc.get('biotype', '-'),
        'EXON': tc.get('exon', '-'),
        'INTRON': tc.get('intron', '-'),
        'HGVSc': tc.get('hgvsc', '-'),
        'HGVSp': tc.get('hgvsp', '-'),
        'cDNA_position': range_str(tc.get('cdna_start'), tc.get('cdna_end')),
        'CDS_position': range_str(tc.get('cds_start'), tc.get('cds_end')),
        'Protein_position': range_str(tc.get('protein_start'), tc.get('protein_end')),
        'Amino_acids': tc.get('amino_acids', '-'),
        'Codons': tc.get('codons', '-'),
        'Existing_variation': cv.get('id', '-'),
        'DISTANCE': tc.get('distance', '-'),
        'STRAND': tc.get('strand', '-'),
        'CANONICAL': 'YES' if tc.get('canonical') else '-',
        'MANE_SELECT': tc.get('mane_select', '-'),
        'TSL': tc.get('tsl', '-'),
        'APPRIS': tc.get('appris', '-'),
        'CCDS': tc.get('ccds', '-'),
        'ENSP': tc.get('protein_id', '-'),
        'SWISSPROT': ','.join(tc.get('swissprot', [])) or '-',
        'SIFT': sift,
        'PolyPhen': pp,
        'AlphaMissense': f"{am['am_class']}({am['am_pathogenicity']})" if am else '-',
        'REVEL': tc.get('revel', '-'),
        'CADD_PHRED': tc.get('cadd_phred', '-'),
        'SpliceAI_DS_max': max(ds) if ds else '-',
        'gnomADe_AF': freqs.get('gnomade', '-'),
        'gnomADg_AF': freqs.get('gnomadg', '-'),
        'CLIN_SIG': allele_clin_sig(cv, allele) or '-',
        'PUBMED': ','.join(str(p) for p in cv.get('pubmed', [])) or '-',
    }


def save_settings(path, n_variants):
    software = get_json(f'{REST}/info/software')
    assembly = get_json(f'{REST}/info/assembly/homo_sapiens')
    sources = get_json(f'{REST}/info/variation/homo_sapiens')
    wanted = ('dbSNP', 'ClinVar', 'gnomAD', '1000Genomes', 'COSMIC')
    with open(path, 'w', encoding='utf-8') as f:
        f.write("Ensembl VEP annotation settings\n")
        f.write("===============================\n")
        f.write(f"date: {date.today().isoformat()}\n")
        f.write(f"endpoint: POST {VEP_URL}\n")
        f.write(f"ensembl release: {software['release']}\n")
        f.write(f"assembly: {assembly['assembly_name']} ({assembly['assembly_accession']})\n")
        f.write(f"geneset update: {assembly.get('genebuild_last_geneset_update', 'N/A')}\n")
        f.write("transcript set: Ensembl/GENCODE (all transcripts, no --pick filtering)\n")
        f.write(f"input variants: {n_variants}\n")
        f.write("\nparameters:\n")
        for k, v in VEP_PARAMS.items():
            f.write(f"  {k} = {v}\n")
        f.write("  sift/polyphen = default (prediction + score)\n")
        f.write("  colocated variants + allele frequencies = default\n")
        f.write("\nannotation sources:\n")
        for s in sources:
            if s.get('name', '').startswith(wanted):
                f.write(f"  {s['name']}: {s.get('version')}\n")


def main():
    if len(sys.argv) < 2:
        print("usage: python vep_annotate.py <file.vcf> [out_dir]")
        sys.exit(1)
    out_dir = sys.argv[2] if len(sys.argv) > 2 else 'vep'
    os.makedirs(out_dir, exist_ok=True)
    records = read_vcf(sys.argv[1])
    print(f"  sending {len(records)} variants to ensembl vep rest ...")
    results = run_vep(records)
    print(f"  got {len(results)} results")

    with open(os.path.join(out_dir, 'vep_output.json'), 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=1)
    rows = [transcript_row(r, tc) for r in results for tc in r.get('transcript_consequences', [])]
    with open(os.path.join(out_dir, 'vep_output.tsv'), 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    save_settings(os.path.join(out_dir, 'vep_settings.txt'), len(records))
    print(f"  transcript rows: {len(rows)}")
    print(f"  saved to {out_dir}/")
    return 0


if __name__ == '__main__':
    sys.exit(main())
