#!/bin/bash

# Check if the folder argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <folder_path>"
    exit 1
fi

# Specify the folder containing your VCF files
vcf_folder="$1"


# Loop through all VCF files in the folder
for vcf_file in "$vcf_folder"/*.vcf*; do
    # Define the full path of the output file with "_filtered.vcf.gz" appended
    output_dir=$(dirname "$vcf_file")
    base_name=$(basename "$vcf_file")
    base_name_no_spaces="${base_name// /_}"
    base_name_no_extension="${base_name_no_spaces%.vcf.gz}"
    filtered_file="${output_dir}/${base_name_no_extension}_filtered.vcf.gz"

    # Use bcftools view to filter by FILTER=PASS and save to the compressed file
    bcftools view --output-type z -i 'FILTER="PASS"' -o "$filtered_file" "$vcf_file"

    # Replace the original file with the filtered file
    mv "$filtered_file" "$vcf_file"
    
    echo "Filtered and replaced: $vcf_file"
done

