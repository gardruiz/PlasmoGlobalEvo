#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory_path> <output_directory"
    exit 1
fi

# Assign the directory path to a variable
directory="$1"
output_directory="$2"

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Directory not found: $directory"
    exit 1
fi

# Loop through all .vcf.gz files in the specified directory
for file in "$directory"/*.vcf.gz; do
    if [ -f "$file" ]; then
        echo "Processing file: $file"
        python Scripts/plot_density.py "$file" "$output_directory"
    else
        echo "No matching files found in directory: $directory"
        exit 1
    fi
done

