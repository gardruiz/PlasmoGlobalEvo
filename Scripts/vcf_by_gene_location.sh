#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: $0 <gene_locations_file>"
  exit 1
fi

gene_locations_file="$1"
gene_locations_output_file="${gene_locations_file%.*}_locations.txt"
gene_locations_output_directory="/home/Daniel/Malaria/Genome/Genes"

# Create the output directory if it doesn't exist
mkdir -p "$gene_locations_output_directory"

while read word; do
  result=$(grep -w "Name=$word" /home/Daniel/Malaria/Genome/Pfalciparum_replace_Pf3D7_MIT_v3_with_Pf_M76611.gff | cut -f1,4,5)
  
  # Append the $word as the first term in each line of the result
  while read line; do
    echo "$word $line" >> "$gene_locations_output_directory/$gene_locations_output_file"
  done <<< "$result"

  # Loop through each line in the input file
  while IFS=$'\t' read -r field1 field2 field3; do
    # Construct the argument using the fields
    arg="$field1:$field2-$field3"

    # Construct the input file name
    input_file="$field1.pf7.vcf.gz"

    # Create the output directory 
    vcf_output_directory="/home/Daniel/Malaria/Variants/vcf"

    # Construct the output file name
    vcf_output_file="vcf_output_directory/$word.vcf"
    
    # Run your command with the constructed argument
    bcftools view -Oz --threads 8  -R "$arg"  $input_file -o $vcf_output_file

    # Check if vcftools completed successfully
    if [ $? -eq 0 ]; then
	    echo "Information for specified samples was successfully extracted and saved to $output_file."
    else
	    echo "Error: vcftools failed to extract information for $word."
	    exit 1
    fi

  done <<< "$result"
done < "$gene_locations_file"




