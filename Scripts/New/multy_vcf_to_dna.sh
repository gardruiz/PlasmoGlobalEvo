#!/bin/bash

# Check if both the VCF folder and genomic region are provided as arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <vcf_dir> <genomic_region>"
    exit 1
fi

vcf_dir="$1"
genomic_region="$2"
sequences_dir="$vcf_dir/sequences"

# Create the sequences directory if it doesn't exist
mkdir -p "$sequences_dir"

# Assuming your script is in the same folder as your .vcf files
script_path="./Scripts/vcf_to_DNA.sh"

# Loop through all .vcf files in the folder
for vcf_file in "$vcf_dir"/*.vcf.gz; do
    if [ -f "$vcf_file" ]; then
        # Normalize the VCF file using bcftools
        normalised_vcf_file="$sequences_dir/$(basename "${vcf_file%.vcf}_normalized.vcf")"
        bcftools norm -a -f Genome/Pfalciparum.genome.fasta -o "$normalised_vcf_file" "$vcf_file"

        # Filter the normalized VCF file using awk
        filtered_vcf_file="$sequences_dir/$(basename "${vcf_file%.vcf}_filtered.vcf")"
        awk '/^#/ || (length($4)==1 && length($5)==1) {print $0}' "$normalised_vcf_file" > "$filtered_vcf_file"
	
	# Compress the filtered VCF file to VCF.gz
        compressed_vcf_file="$sequences_dir/$(basename "${vcf_file%.vcf}_filtered.vcf.gz")"
        bgzip -c "$filtered_vcf_file" > "$compressed_vcf_file"
        tabix -p vcf "$compressed_vcf_file"	
        
	# Run the vcf_to_DNA.sh script with the filtered VCF file and genomic region
        echo "Processing $vcf_file for genomic region $genomic_region..."
        bash "$script_path" "$compressed_vcf_file" "$genomic_region"

        # Delete intermediate files after using them
        rm "$normalised_vcf_file" "$filtered_vcf_file" "$compressed_vcf_file" "$compressed_vcf_file.csi" "$compressed_vcf_file.tbi"
    fi
done


