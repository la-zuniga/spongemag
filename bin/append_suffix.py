import sys

def append_suffix(input_file, suffix, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                sequence_name, genome_number = parts
                outfile.write(f"{sequence_name}\t{genome_number}_{suffix}\n")
            else:
                print(f"Skipping malformed line: {line.strip()}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <suffix> <output_file>")
    else:
        _, input_file, suffix, output_file = sys.argv
        append_suffix(input_file, suffix, output_file)

