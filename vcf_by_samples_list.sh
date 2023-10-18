#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <vcf_directory/vcf_file> <Sample_lists_directory>"
    exit 1
fi

vcf_path="$1"
parent_dir="$2"

# Check if the parent directory exists
if [ ! -d "$parent_dir" ]; then
    echo "Parent directory $parent_dir does not exist."
    exit 1
fi

# Extract the VCF directory and file name
vcf_dir=$(dirname "$vcf_path")
vcf_file=$(basename "$vcf_path")

# Check if the VCF file exists
if [ ! -f "$vcf_dir/$vcf_file" ]; then
    echo "VCF file $vcf_dir/$vcf_file does not exist."
    exit 1
fi

# Loop through all subdirectories (sample list folders) in the parent directory
for sample_list_dir in "$parent_dir"/*/; do
    # Check if it's a directory
    if [ -d "$sample_list_dir" ]; then
        # Get the name of the samples directory (the last part of the path)
        samples_dir=$(basename "$sample_list_dir")

        # Loop through all sample list files in the current sample list directory
        for sample_list_file in "$sample_list_dir"/*.txt; do
            # Extract the sample file name without extension
            sample_name=$(basename "$sample_list_file" .txt)

            # Define the output directory based on the VCF file directory and the samples directory name
            output_dir="$vcf_dir/$samples_dir"

            # Create the output directory if it doesn't exist
            mkdir -p "$output_dir"

            # Define the output file path
            output_file="$output_dir/$sample_name.vcf"

            # Run the script with the VCF file, sample list file, and output file
            Scripts/vcf_by_sample.sh "$vcf_dir/$vcf_file" "$sample_list_file" "$output_file"

            # Message to indicate completion
            echo "Processed $sample_list_file"
        done
    fi
done
    
