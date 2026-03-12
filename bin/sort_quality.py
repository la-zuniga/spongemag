#!/usr/bin/env python3
import argparse

def filter_quality(input_file, output_file):
    with open(input_file, 'r') as file:
        header = next(file)
        data = file.readlines()

    filtered_lines = []
    for line in data:
        columns = line.strip().split('\t')
        if len(columns) < 3:
            continue
        name = columns[0]
        completeness = float(columns[1])
        contamination = float(columns[2])
        if completeness >= 95.0 and contamination < 5.0:
            filtered_lines.append(line)

    with open(output_file, 'w') as outfile:
        outfile.write(header)
        outfile.writelines(filtered_lines)

    print(f"Filtered results saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter quality report based on completeness and contamination.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input TSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output filtered TSV file")
    args = parser.parse_args()

    filter_quality(args.input, args.output)






















