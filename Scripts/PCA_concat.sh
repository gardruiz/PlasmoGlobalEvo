#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_folder> <output_directory>/<output_base_name>"
    exit 1
fi

# Set the input folder and output base name from the command line arguments
input_folder="$1"
out_base="$2"

# Concatenate VCF files in the input folder using bcftools
bcftools concat -o "${out_base}_concatenated.vcf.gz" "$input_folder"/*.vcf.gz

# Check if the concatenation was successful
if [ $? -eq 0 ]; then
    echo "VCF concatenation successful. Proceeding with PCA."

    # Perform linkage pruning and PCA using the concatenated VCF file
    plink2 --memory 29000 \ 
          --vcf "${out_base}_concatenated.vcf.gz" \
          --double-id \
          --allow-extra-chr \
          --set-missing-var-ids @:# \
          --max-alleles 2 \
          --make-bed \
          --indep-pairwise 50 10 0.1 \
          --mind 0.1 \
          --rm-dup \
          --genotyping-rate \
          --out "$out_base"

    plink2 --bfile "${out_base}" \
          --double-id \
          --allow-extra-chr \
          --set-missing-var-ids @:# \
          --max-alleles 2 \
          --exclude "${out_base}.prune.out" \
          --make-bed \
          --maf 0.05 \
          --geno 0.1 \
          --pca approx \
          --genotyping-rate \
          --out "$out_base"

else
    echo "VCF concatenation failed. Please check for errors."
fi

