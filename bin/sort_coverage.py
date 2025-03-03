import argparse

# Function to filter coverage table
def filter_coverage_by_kmer(coverage_file, kmer_file, output_file):
    # Read kmer file and extract contig names into a set
    with open(kmer_file, 'r') as kmer:
        kmer_contigs = set(line.split(',')[0] for line in kmer)  # Split by comma and extract the first column (contig name)

    # Read coverage table and filter based on kmer contigs
    with open(coverage_file, 'r') as coverage, open(output_file, 'w') as output:
        # Write the header to the output file
        output.write("contig\tcov_mean_sample\n")
        # Process each line in the coverage file
        for line in coverage:
            contig_name = line.split('\t')[0]  # Assuming the first column is the contig name in TSV
            if contig_name in kmer_contigs:
                output.write(line)

# Set up argument parsing
def main():
    parser = argparse.ArgumentParser(description="Filter coverage table based on kmer file.")
    parser.add_argument('-c', '--coverage', required=True, help="Path to the coverage table (TSV file).")
    parser.add_argument('-k', '--kmer', required=True, help="Path to the kmer file (CSV file).")
    parser.add_argument('-o', '--output', required=True, help="Path to save the filtered output (TSV file).")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with the parsed arguments
    filter_coverage_by_kmer(args.coverage, args.kmer, args.output)

# Run the script if it is executed directly
if __name__ == "__main__":
    main()
