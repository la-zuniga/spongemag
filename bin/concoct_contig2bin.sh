#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

input_file="$1"
output_file="$2"

# Convert commas to tabs and save to output file
sed 's/,/\t/g' "$input_file" > "$output_file"

echo "Conversion complete: $output_file"
