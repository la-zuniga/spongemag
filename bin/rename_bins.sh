#!/bin/bash

# Define the annotation directories
for tool in concoct metabat; do
  annotation_dir="prokka_${tool}_annotation"

  # Loop through each genome folder
  for genome in "$annotation_dir"/*/; do
    genome_name=$(basename "$genome")

    # Loop through each Prokka output file and rename it
    for file in "$genome"PROKKA_*; do
      ext="${file##*.}"
      mv "$file" "$genome/${genome_name}_${tool}.${ext}"
    done
  done
done
