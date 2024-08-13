#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 directory_path output_file"
    exit 1
fi

directory_path="$1"
output_file="$2"

# Check if the output file already exists
if [ -e "$output_file" ]; then
    echo "Output file already exists. Removing it."
    rm "$output_file"
fi

# List .gz files in the specified directory and save full file paths
ls "$directory_path"/*.vcf.gz > "$output_file"



