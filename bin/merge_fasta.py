from collections import defaultdict
import sys

def process_fasta(input_fasta, output_fasta):
    sequences = defaultdict(str)  # Dictionary to store concatenated sequences
    
    with open(input_fasta, 'r') as infile:
        header = None
        sequence = []
        
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] += ''.join(sequence)  # Store previous sequence
                
                header = line.split(".concoct_part_")[0]  # Remove suffix
                sequence = []  # Reset sequence buffer
            else:
                sequence.append(line)
        
        if header:
            sequences[header] += ''.join(sequence)  # Store last sequence
    
    with open(output_fasta, 'w') as outfile:
        for header, seq in sequences.items():
            outfile.write(f"{header}\n")
            outfile.write(f"{seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output.fasta")
    else:
        process_fasta(sys.argv[1], sys.argv[2])

