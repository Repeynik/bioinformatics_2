#!/usr/bin/env python3
import csv
import json
import sys
from collections import Counter

SEVERITY = [
    'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained',
    'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification',
    'feature_elongation', 'feature_truncation', 'inframe_insertion', 'inframe_deletion',
    'missense_variant', 'protein_altering_variant', 'splice_donor_5th_base_variant',
    'splice_region_variant', 'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant', 'start_retained_variant', 'stop_retained_variant',
    'synonymous_variant', 'coding_sequence_variant', 'mature_miRNA_variant',
    '5_prime_UTR_variant', '3_prime_UTR_variant', 'non_coding_transcript_exon_variant',
    'intron_variant', 'NMD_transcript_variant', 'non_coding_transcript_variant',
    'coding_transcript_variant', 'upstream_gene_variant', 'downstream_gene_variant',
    'TFBS_ablation', 'TFBS_amplification', 'TF_binding_site_variant',
    'regulatory_region_ablation', 'regulatory_region_amplification', 'regulatory_region_variant',
    'intergenic_variant', 'sequence_variant',
]
RANK = {term: i for i, term in enumerate(SEVERITY)}
COLUMNS = [
    'input_id', 'location', 'ref_alt', 'gene', 'vep_consequence', 'impact', 'transcript',
    'refseq_mane', 'transcript_rule', 'hgvsc', 'hgvsp', 'vep_most_severe_any',
    'most_severe_transcripts', 'sift', 'polyphen', 'alphamissense', 'cadd_phred',
    'rsid', 'gnomad_af', 'gnomad_grpmax_af', 'gnomad_filters',
    'clinvar', 'clinvar_classification', 'clinvar_review_status', 'clinvar_accession',
    'clinvar_other_records', 'vep_clin_sig',
]


def read_tsv(path):
    with open(path, 'r', encoding='utf-8') as f:
        return {row['id']: row for row in csv.DictReader(f, delimiter='\t')}


def severity(tc):
    return min(RANK.get(t, len(SEVERITY)) for t in tc.get('consequence_terms', ['sequence_variant']))


def choose_transcript(tcs):
    mane = [t for t in tcs if t.get('mane_select')]
    if mane:
        return min(mane, key=severity), 'MANE Select'
    canonical = [t for t in tcs if t.get('canonical') and t.get('biotype') == 'protein_coding']
    if canonical:
        return min(canonical, key=severity), 'Ensembl canonical'
    return min(tcs, key=severity), 'most severe (no MANE/canonical)'


def short_hgvs(value):
    return value.split(':', 1)[1] if value and ':' in value else '-'


def pred(tc, name):
    if f'{name}_prediction' not in tc:
        return '-'
    return f"{tc[f'{name}_prediction']}({tc[f'{name}_score']})"


def rs_colocated(result):
    for cv in result.get('colocated_variants', []):
        if cv.get('id', '').startswith('rs'):
            return cv
    return {}


def vep_clin_sig(cv, allele):
    values = [part.split(':', 1)[1] for part in cv.get('clin_sig_allele', '').split(';')
              if part.startswith(f'{allele}:')]
    return ','.join(values) if values else '-'


def most_severe_note(result, chosen):
    msc = result['most_severe_consequence']
    if msc in chosen.get('consequence_terms', []):
        return '-'
    hits = [t for t in result['transcript_consequences'] if msc in t.get('consequence_terms', [])]
    by_gene = Counter(t.get('gene_symbol', '?') for t in hits)
    total = Counter(t.get('gene_symbol', '?') for t in result['transcript_consequences'])
    parts = []
    for gene, n in by_gene.items():
        ids = [t['transcript_id'] for t in hits if t.get('gene_symbol', '?') == gene]
        shown = ','.join(ids) if n <= 2 else f'{n} transcripts'
        parts.append(f'{gene} {n}/{total[gene]} ({shown})')
    return '; '.join(parts)


def build_row(result, clinvar, gnomad):
    tcs = result['transcript_consequences']
    tc, rule = choose_transcript(tcs)
    allele = tc.get('variant_allele', '')
    ref, alt = result['allele_string'].split('/')
    cv = rs_colocated(result)
    am = tc.get('alphamissense') or {}
    cl = clinvar.get(result['id'], {})
    gn = gnomad.get(result['id'], {})
    grpmax = gn.get('grpmax_af', '-')
    if grpmax != '-':
        grpmax = f"{grpmax} ({gn.get('grpmax_group')})"
    others = [o for o in cl.get('other_records_at_position', '-').split(' | ') if o != '-']
    return {
        'input_id': result['id'],
        'location': f"{result['seq_region_name']}:{result['start']}",
        'ref_alt': f'{ref}>{alt}',
        'gene': tc.get('gene_symbol', '-'),
        'vep_consequence': ','.join(tc.get('consequence_terms', [])),
        'impact': tc.get('impact', '-'),
        'transcript': tc.get('transcript_id', '-'),
        'refseq_mane': tc.get('mane_select', '-'),
        'transcript_rule': rule,
        'hgvsc': short_hgvs(tc.get('hgvsc')),
        'hgvsp': short_hgvs(tc.get('hgvsp')),
        'vep_most_severe_any': result['most_severe_consequence'],
        'most_severe_transcripts': most_severe_note(result, tc),
        'sift': pred(tc, 'sift'),
        'polyphen': pred(tc, 'polyphen'),
        'alphamissense': f"{am['am_class']}({am['am_pathogenicity']})" if am else '-',
        'cadd_phred': tc.get('cadd_phred', '-'),
        'rsid': cv.get('id', '-'),
        'gnomad_af': gn.get('joint_af', '-'),
        'gnomad_grpmax_af': grpmax,
        'gnomad_filters': gn.get('joint_filters', '-'),
        'clinvar': cl.get('in_clinvar', 'no'),
        'clinvar_classification': cl.get('classification', '-'),
        'clinvar_review_status': cl.get('review_status', '-'),
        'clinvar_accession': cl.get('accession', '-'),
        'clinvar_other_records': str(len(others)) if others else '-',
        'vep_clin_sig': vep_clin_sig(cv, allele),
    }


def af_bin(value):
    if value == '-':
        return 'not in gnomAD'
    x = float(value)
    if x >= 0.05:
        return 'common (AF >= 5%)'
    if x >= 0.01:
        return 'low frequency (1% <= AF < 5%)'
    if x >= 0.001:
        return 'rare (0.1% <= AF < 1%)'
    return 'very rare (AF < 0.1%)'


def print_summary(rows):
    blocks = [
        ('consequence (chosen transcript)', Counter(r['vep_consequence'] for r in rows)),
        ('most severe consequence (any transcript)', Counter(r['vep_most_severe_any'] for r in rows)),
        ('transcript rule', Counter(r['transcript_rule'] for r in rows)),
        ('gnomAD joint AF', Counter(af_bin(r['gnomad_af']) for r in rows)),
        ('ClinVar classification', Counter(r['clinvar_classification'] for r in rows)),
    ]
    for title, counter in blocks:
        print(f"=== {title} ===")
        for key, n in counter.most_common():
            print(f"  {key}: {n}")
        print()
    diff = [r for r in rows if r['most_severe_transcripts'] != '-']
    print("=== most severe consequence differs from chosen transcript ===")
    for r in diff:
        print(f"  {r['input_id']} {r['gene']}: {r['vep_consequence']} -> "
              f"{r['vep_most_severe_any']} in {r['most_severe_transcripts']}")
    print()
    mismatch = [r for r in rows if (r['clinvar'] == 'yes') != (r['vep_clin_sig'] != '-')]
    print("=== ClinVar (NCBI) vs VEP clin_sig ===")
    for r in mismatch:
        print(f"  {r['input_id']}: clinvar={r['clinvar_classification']}, vep={r['vep_clin_sig']}")


def main():
    if len(sys.argv) < 5:
        print("usage: python build_table.py <vep_output.json> <clinvar_results.tsv> "
              "<gnomad_results.tsv> <out.tsv>")
        sys.exit(1)
    with open(sys.argv[1], 'r', encoding='utf-8') as f:
        vep = json.load(f)
    clinvar = read_tsv(sys.argv[2])
    gnomad = read_tsv(sys.argv[3])
    rows = [build_row(r, clinvar, gnomad) for r in vep]
    with open(sys.argv[4], 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  rows: {len(rows)}, saved to {sys.argv[4]}")
    print()
    print_summary(rows)
    return 0


if __name__ == '__main__':
    sys.exit(main())
