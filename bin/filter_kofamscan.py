#!/usr/bin/env python3
import argparse
import os
import glob

def parse_kofamscan(filepath):
    hits = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            if parts[0] == '*':
                significant = True
                gene, ko = parts[1], parts[2]
                threshold, score, evalue = float(parts[3]), float(parts[4]), float(parts[5])
                definition = ' '.join(parts[6:])
            else:
                significant = False
                gene, ko = parts[0], parts[1]
                threshold, score, evalue = float(parts[2]), float(parts[3]), float(parts[4])
                definition = ' '.join(parts[5:])

            hits.append({
                'gene':        gene,
                'ko':          ko,
                'threshold':   threshold,
                'score':       score,
                'evalue':      evalue,
                'definition':  definition,
                'significant': significant,
                'score_ratio': score / threshold if threshold > 0 else 0
            })
    return hits

def classify_hit(hit, evalue_cutoff=1e-5, min_putative_score=50):
    """
    Tier 1 — confirmed:  passes KOFamScan threshold AND E-value ≤ cutoff
    Tier 2 — putative:   below threshold but score ≥ min_putative_score AND E-value ≤ cutoff
                         (captures divergent homologs common in environmental metagenomics)
    Tier 3 — weak:       everything else
    """
    if hit['significant'] and hit['evalue'] <= evalue_cutoff:
        return 'confirmed'
    elif (not hit['significant']
          and hit['score'] >= min_putative_score
          and hit['evalue'] <= evalue_cutoff):
        return 'putative'
    else:
        return 'weak'

def filter_kofamscan(input_dir, output_file,
                     min_tier='putative',
                     evalue_cutoff=1e-5,
                     min_putative_score=50):

    tier_rank = {'confirmed': 0, 'putative': 1, 'weak': 2}
    min_rank  = tier_rank[min_tier]

    all_hits = []
    for filepath in sorted(glob.glob(os.path.join(input_dir, '*.out'))):
        bin_name = os.path.splitext(os.path.basename(filepath))[0]
        for hit in parse_kofamscan(filepath):
            tier = classify_hit(hit, evalue_cutoff, min_putative_score)
            if tier_rank[tier] <= min_rank:
                hit['bin']  = bin_name
                hit['tier'] = tier
                all_hits.append(hit)

    with open(output_file, 'w') as out:
        out.write('bin\tgene\tko\tthreshold\tscore\tscore_ratio\tevalue\ttier\tdefinition\n')
        for hit in sorted(all_hits, key=lambda x: (x['bin'], x['gene'], x['ko'])):
            out.write(
                f"{hit['bin']}\t{hit['gene']}\t{hit['ko']}\t"
                f"{hit['threshold']}\t{hit['score']}\t{hit['score_ratio']:.3f}\t"
                f"{hit['evalue']}\t{hit['tier']}\t{hit['definition']}\n"
            )

    # Print a summary to stdout so it shows up in the Nextflow log
    confirmed = sum(1 for h in all_hits if h['tier'] == 'confirmed')
    putative  = sum(1 for h in all_hits if h['tier'] == 'putative')
    print(f"Total hits written: {len(all_hits)}")
    print(f"  Confirmed (above threshold, E ≤ {evalue_cutoff}): {confirmed}")
    print(f"  Putative  (score ≥ {min_putative_score}, E ≤ {evalue_cutoff}):  {putative}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Filter KOFamScan output by confidence tier.'
    )
    parser.add_argument('-i', '--input_dir',   required=True,
                        help='Directory of KOFamScan .out files')
    parser.add_argument('-o', '--output',      required=True,
                        help='Output filtered TSV')
    parser.add_argument('--min_tier',          default='putative',
                        choices=['confirmed', 'putative'],
                        help='Minimum tier to include [default: putative]')
    parser.add_argument('--evalue_cutoff',     type=float, default=1e-5,
                        help='Maximum E-value to retain [default: 1e-5]')
    parser.add_argument('--min_putative_score', type=float, default=50,
                        help='Minimum score for putative hits [default: 50]')
    args = parser.parse_args()

    filter_kofamscan(
        args.input_dir,
        args.output,
        args.min_tier,
        args.evalue_cutoff,
        args.min_putative_score
    )
