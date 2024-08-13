#!/bin/bash

# Directory containing the VCF files
vcf_dir="/media/Daniel/ADATA HM8001/Malaria/Variants/"

# Check if the directory exists
if [ ! -d "$vcf_dir" ]; then
  echo "Directory does not exist: $vcf_dir"
  exit 1
fi

# Change to the directory
cd "$vcf_dir"

# List all VCF files in the directory
vcf_files=(*.vcf.gz)

# Number of CPU cores available for parallel processing
num_cores=$(nproc)

# Function to index a VCF file
index_vcf() {
  tabix -p vcf "$1"
  echo "Indexed $1"
}

# Iterate through the VCF files and index them in parallel
for vcf_file in "${vcf_files[@]}"; do
  index_vcf "$vcf_file" &
done

# Wait for all indexing processes to finish
wait

echo "Indexing of VCF files in $vcf_dir completed."

