import sys

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith(">"):
                # Split the header at the first space and keep only the first part
                header = line.split(" ")[0]
                outfile.write(header + "\n")
            else:
                # Write the sequence line as it is
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_fasta> <output_fasta>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_fasta(input_file, output_file)

