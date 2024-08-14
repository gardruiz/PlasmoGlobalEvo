#!/bin/bash

# Check if both the VCF folder and genomic region are provided as arguments
if [ "$#" -ne 3 ]; then
	echo "Usage: $0 <vcf_dir> <genomic_region> <bed_file>"
    exit 1
fi

vcf_dir="$1"
genomic_region="$2"
bed_file="$3"
sequences_dir="$vcf_dir/sequences"

# Create the sequences directory if it doesn't exist
mkdir -p "$sequences_dir"

# Path to scripts
dna_script="Scripts/vcf_to_DNA.sh"
cds_script="Scripts/CDS.py"

# Loop through all .vcf files in the folder
for vcf_file in "$vcf_dir"/*.vcf.gz; do
    if [ -f "$vcf_file" ]; then

    	# Run bcftools and count the lines
	variants=$(bcftools view -H "$vcf_file" | wc -l)

	# Print the result to the screen
	echo "Number of variants in the VCF file: $variants"
	
	# no_indels_vcf_file="$sequences_dir/$(basename "${vcf_file%.vcf.gz}_noindels.vcf")"
	filtered_vcf_file="$sequences_dir/$(basename "${vcf_file%.vcf.gz}_filtered.vcf")"

        # Filter based on reference length <= 3 using awk
        bcftools view -h "$vcf_file" > tmp/header.vcf
        bcftools view -H "$vcf_file" | awk 'length($4) <= 3' > tmp/filtered_variants.vcf
        cat tmp/header.vcf tmp/filtered_variants.vcf | bcftools view -Oz -o "$filtered_vcf_file"


        # Normalize the VCF file using bcftools
        normalised_vcf_file="$sequences_dir/$(basename "${vcf_file%.vcf.gz}_normalized.vcf")"
        bcftools norm -a -f Genome/Pfalciparum.genome.fasta -o "$normalised_vcf_file" "$filtered_vcf_file"

	echo "Number of variants after filtering indels and multy-allelic variants after normalization: $no_ambi"

        # Compress the filtered VCF file to VCF.gz
        compressed_vcf_file="$sequences_dir/$(basename "${normalised_vcf_file%_normalized.vcf}.vcf.gz")"
        bgzip -c "$normalised_vcf_file" > "$compressed_vcf_file"
        tabix -p vcf "$compressed_vcf_file"     
        
        # Run the vcf_to_DNA.sh script with the filtered VCF file and genomic region
        echo "Processing $vcf_file for genomic region $genomic_region..."
        bash "$dna_script" "$compressed_vcf_file" "$genomic_region"

	# Run the CDS.py script with the output of vcf_to_DNA.sh
        echo "Obtaining CDS from DNA..."
	cds_file="$sequences_dir/$(basename "${vcf_file%.vcf.gz}_CDS.fasta")"
	dna_output="$sequences_dir/$(basename "${vcf_file%.vcf.gz}.fasta")"
	
	python3 "$cds_script" "$dna_output" "$bed_file" "$cds_file"

        # Delete intermediate files after using them
	rm  "$normalised_vcf_file" "$compressed_vcf_file" "$compressed_vcf_file.tbi" "$dna_output" "$filtered_vcf_file" "tmp/header.vcf" "tmp/filtered_variants.vcf"   
    fi
done



