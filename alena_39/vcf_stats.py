#!/usr/bin/env python3
import sys
from collections import Counter

TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}


def read_vcf(path):
    header = []
    records = []
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('##'):
                header.append(line)
                continue
            if line.startswith('#'):
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = line.split('\t')[:8]
            records.append({
                'chrom': chrom, 'pos': int(pos), 'id': vid, 'ref': ref, 'alt': alt,
                'qual': qual, 'filter': filt, 'info': info,
            })
    return header, records


def main():
    if len(sys.argv) < 2:
        print("usage: python vcf_stats.py <file.vcf>")
        sys.exit(1)
    header, records = read_vcf(sys.argv[1])
    chroms = Counter(r['chrom'] for r in records)
    snv = [r for r in records if len(r['ref']) == 1 and len(r['alt']) == 1]
    ts = sum(1 for r in snv if (r['ref'], r['alt']) in TRANSITIONS)
    tv = len(snv) - ts
    spectrum = Counter(f"{r['ref']}>{r['alt']}" for r in snv)

    print("=== vcf header ===")
    for line in header:
        print(f"  {line}")
    print()
    print("=== records ===")
    print(f"  total records:  {len(records)}")
    print(f"  snv:            {len(snv)}")
    print(f"  multiallelic:   {sum(1 for r in records if ',' in r['alt'])}")
    print(f"  with rsID:      {sum(1 for r in records if r['id'].startswith('rs'))}")
    print(f"  QUAL = '.':     {sum(1 for r in records if r['qual'] == '.')}")
    print(f"  FILTER = '.':   {sum(1 for r in records if r['filter'] == '.')}")
    print(f"  INFO = '.':     {sum(1 for r in records if r['info'] == '.')}")
    print()
    print("=== chromosomes ===")
    for chrom in sorted(chroms, key=lambda c: (len(c), c)):
        print(f"  chr{chrom}: {chroms[chrom]}")
    print()
    print("=== transitions / transversions ===")
    print(f"  transitions:    {ts}")
    print(f"  transversions:  {tv}")
    print(f"  ts/tv:          {ts / tv:.2f}" if tv else "  ts/tv:          N/A")
    print()
    print("=== substitutions ===")
    for sub in sorted(spectrum):
        kind = 'ts' if tuple(sub.split('>')) in TRANSITIONS else 'tv'
        print(f"  {sub} ({kind}): {spectrum[sub]}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
