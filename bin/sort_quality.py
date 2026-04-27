#!/usr/bin/env python3
"""
Split CheckM2 quality_report.tsv into MIMAG-style tiers.

Outputs three TSVs into --outdir:
  high_quality.tsv    completeness >= 90  AND contamination < 5
  medium_quality.tsv  completeness >= 50  AND contamination < 10  (and not HQ)
  low_quality.tsv     everything else

Each TSV keeps the original header so downstream awk `{print $1}` still grabs
the bin name.
"""
import argparse
import os


HQ_COMPLETENESS = 90.0
HQ_CONTAMINATION = 5.0
MQ_COMPLETENESS = 50.0
MQ_CONTAMINATION = 10.0


def tier_for(completeness: float, contamination: float) -> str:
    if completeness >= HQ_COMPLETENESS and contamination < HQ_CONTAMINATION:
        return "high"
    if completeness >= MQ_COMPLETENESS and contamination < MQ_CONTAMINATION:
        return "medium"
    return "low"


def split_quality(input_file: str, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)

    with open(input_file) as f:
        header = next(f)
        rows = f.readlines()

    buckets = {"high": [], "medium": [], "low": []}
    for line in rows:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 3:
            continue
        try:
            completeness = float(cols[1])
            contamination = float(cols[2])
        except ValueError:
            continue
        buckets[tier_for(completeness, contamination)].append(line)

    for tier, lines in buckets.items():
        out_path = os.path.join(outdir, f"{tier}_quality.tsv")
        with open(out_path, "w") as out:
            out.write(header)
            out.writelines(lines)
        print(f"  {tier}: {len(lines)} bin(s) -> {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split CheckM2 quality_report.tsv into HQ/MQ/LQ TSVs."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Path to CheckM2 quality_report.tsv")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Directory to write {high,medium,low}_quality.tsv into")
    args = parser.parse_args()

    split_quality(args.input, args.outdir)
