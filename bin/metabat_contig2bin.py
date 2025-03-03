import os
import argparse
import re

def extract_contig_info(fasta_dir, output_file):
    header_pattern = re.compile(r">(\S+)")  # Regex to capture the contig name
    
    with open(output_file, 'w') as out_f:
        for fasta_file in os.listdir(fasta_dir):
            if fasta_file.endswith(".fa") or fasta_file.endswith(".fasta"):
                fasta_path = os.path.join(fasta_dir, fasta_file)
                genome_name = os.path.splitext(os.path.basename(fasta_file))[0]  # Remove file extension
                
                with open(fasta_path, 'r') as f:
                    for line in f:
                        match = header_pattern.match(line)
                        if match:
                            contig_name = match.group(1)
                            out_f.write(f"{contig_name}\t{genome_name}\n")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract contig names and genome sources from FASTA files.")
    parser.add_argument("input_dir", help="Path to the directory containing FASTA files.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    args = parser.parse_args()
    
    extract_contig_info(args.input_dir, args.output_file)
