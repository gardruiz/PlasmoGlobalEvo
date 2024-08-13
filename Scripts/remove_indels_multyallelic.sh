#!/bin/bash

# Check if the directory argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <vcf_directory>"
    exit 1
fi

# Get the directory containing VCF files
vcf_dir="$1"

# Create the "SNPs" folder if it doesn't exist
mkdir -p "$vcf_dir/SNPs"

# Iterate over each VCF file in the directory
for vcf_file in "$vcf_dir"/*.vcf.gz; do
    # Extract filename without extension
    filename=$(basename -- "$vcf_file")
    filename_no_ext="${filename%.vcf.gz}"

    # Run bcftools view to filter multiallelic sites and indels
    bcftools view -e 'N_ALT > 1 || TYPE="indel"' "$vcf_file" -o "$vcf_dir/SNPs/$filename_no_ext.vcf.gz"
done

