#!/bin/bash

# Define the directories
for tool in prokka_concoct_annotation prokka_metabat_annotation; do
  if [ -d "$tool" ]; then
    # Extract the binning tool name (remove 'prokka_' and '_annotation')
    binning_tool="${tool#prokka_}"
    binning_tool="${binning_tool%_annotation}"

    for folder in "$tool"/*; do
      if [ -d "$folder" ]; then
        bin_id=$(basename "$folder") # Extract genome number
        for file in "$folder"/*; do
          ext="${file##*.}"                           # Extract file extension
          new_name="${bin_id}_${binning_tool}.${ext}" # Format: <genome_number>_<binning_tool>.<ext>
          mv "$file" "$folder/$new_name"
        done
      fi
    done
  fi
done
