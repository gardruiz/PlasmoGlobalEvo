#!/bin/bash

# Check if folder path is provided as argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <folder_path>"
    exit 1
fi

# Get folder path from argument
folder="$1"

# Check if the folder exists
if [ ! -d "$folder" ]; then
    echo "Folder '$folder' does not exist."
    exit 1
fi

# Create 'clean' directory if it doesn't exist
clean_dir="$folder/clean"
if [ ! -d "$clean_dir" ]; then
    mkdir "$clean_dir"
fi

# Loop through each file in the folder
for file in "$folder"/*.out; do
    # Extract input file name without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"

    # Output file name
    output_file="$clean_dir/${filename_no_ext}.tab"

    # Run awk command for the current file
    awk '/seq\. seq./,/LWL85, LPB93 & LWLm methods/ {
        if (!/LWL85, LPB93 & LWLm methods/) 
            print
    }' "$file" > "$output_file"
done


