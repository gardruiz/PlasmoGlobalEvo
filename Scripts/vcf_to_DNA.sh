#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.vcf> <genomic_region>"
    exit 1
fi

reference="Genome/Pfalciparum.genome.fasta"
vcf="$1"
genomic_region="$2"

# Create a temporary file for subset_reference
subset_reference=$(samtools faidx "$reference" "$genomic_region")

# Extract VCF file name without extension
vcf_filename=$(basename -- "$vcf")
vcf_filename_noext="${vcf_filename%%.*}"  # Use %% to remove both .vcf and .vcf.gz extensions
vcf_dir=$(dirname -- "$vcf")

# Create the sequences directory inside the VCF file directory if it doesn't exist
output_dir="$vcf_dir"

# Check if the index file exists
if [ ! -e "${vcf}.tbi" ]; then
    # Generate an index file for the VCF
    bcftools index "$vcf"
    echo "Index file generated for $vcf"
else
    echo "Index file already exists for $vcf"
fi


# Path to the Python script
python_script="Scripts/remove_indels.py"

# Iterate through each sample in the VCF
for sample in $(bcftools query -l "$vcf"); do

    # Generate consensus sequence for the current sample using the filtered VCF
    bcftools consensus -H A --mark-ins lc --mark-del "*" -s "$sample" -f <(echo "$subset_reference") "$vcf" -o "$output_dir/${sample}_consensus.fa"

    # Log the length of the generated consensus sequence
    consensus_length=$(grep -v '>' "$output_dir/${sample}_consensus.fa" | tr -d '\n' | wc -c)
    echo "Length of consensus sequence for $sample: $consensus_length"

    # Run the Python script to process the consensus sequence
    python3 "$python_script" "$output_dir/${sample}_consensus.fa" "$subset_reference"  "$output_dir/${sample}_final_consensus.fas"

done

# Concatenate all individual consensus sequences into a single file
output_combined="$output_dir/$vcf_filename_noext.fasta"

# Remove the existing combined file if it exists
if [ -e "$output_combined" ]; then
    rm "$output_combined"
    echo "Removed existing combined consensus file: $output_combined"
fi

cat $output_dir/*.fas > "$output_combined"
rm $output_dir/*.fas
rm $output_dir/*.fa
