#!/usr/bin/env python3
import os
import re
import argparse

def extract_contig_info(bins_dir, output_file):
    header_pattern = re.compile(r">(\S+)")
    concoct_suffix = re.compile(r'\.concoct_part_\d+$')

    with open(output_file, 'w') as out_f:
        for fasta_file in os.listdir(bins_dir):
            if fasta_file.endswith(".fa") or fasta_file.endswith(".fasta"):
                fasta_path = os.path.join(bins_dir, fasta_file)
                # Bin name is the filename without extension, e.g. "1" from "1.fa"
                bin_name = os.path.splitext(fasta_file)[0]

                with open(fasta_path, 'r') as f:
                    for line in f:
                        match = header_pattern.match(line)
                        if match:
                            raw_contig = match.group(1)
                            # Strip .concoct_part_N suffix to match other binners
                            clean_contig = concoct_suffix.sub('', raw_contig)
                            out_f.write(f"{clean_contig}\t{bin_name}_metabinner\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate contig-to-bin TSV from MetaBinner output bins."
    )
    parser.add_argument("input_dir", help="Path to MetaBinner bins directory (containing .fa files)")
    parser.add_argument("output_file", help="Path to output TSV file")
    args = parser.parse_args()

    extract_contig_info(args.input_dir, args.output_file)
