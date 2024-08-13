#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_vcf> <output_directory>/<output_base_name>"
    exit 1
fi

# Set the input VCF file from the first command line argument
vcf="$1"

# Set the output directory and base name from the second command line argument
out_base="$2"
<<"COMENT"
# Perform linkage pruning - i.e. identify prune sites
plink2 --vcf $vcf \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --max-alleles 2 \
      --indep-pairwise 50 10 0.1 \
      --rm-dup \
      --genotyping-rate \
      --make-bed \
      --memory 28000 \
      --debug \
      --out "$out_base"
COMENT

# Use the output base name of the first script as input for --exclude in the second script
plink2 --bfile "${out_base}" \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --max-alleles 2 \
      --exclude $out_base.prune.out \
      --remove Config/exclude.txt \
      --mind 0.01\
      --geno 0.01 \
      --maf 0.01 \
      --pca approx allele-wts \
      --out "$out_base"




#      --remove exclude.txt \
