#!/bin/bash

if [ $# -ne 5 ]; then
  echo "Usage: $0 <gene_names_file>  <sample_list_directory> <vcf_file> <GFF_file> <use_gene_name>" 
  exit 1
fi

gene_locations_file="$1"
sample_list_parent_dir="$2"
input_file="$3"
gff_file="$4"
gene_locations_output_file="${gene_locations_file%.*}_locations.txt"
gene_locations_output_directory=$(dirname "$1")
use_gene_name="$5"

# Create the output directory if it doesn't exist
mkdir -p "$gene_locations_output_directory"

# Check if the output file exists and handle it
if [ -e "$gene_locations_output_file" ]; then
    echo "File $gene_locations_output_file already exists. Replacing it..."
    rm "$gene_locations_output_file"
fi

while read word; do
     if [ $use_gene_name = "True" ]; then
        result=$(grep "Name=$word;" "$gff_file" | cut -f1,4,5)
  
    else
        result=$(grep "ID=$word;" "$gff_file" | cut -f1,4,5)
    fi

    # Append the $word as the first term in each line of the result
    while read line; do
    echo  -e "$word\t$line" >> "$gene_locations_output_file"
    done <<< "$result"

    # Loop through each line in the genes file
    while IFS=$'\t' read -r field1 field2 field3; do
    # Construct the argument using the fields
    arg="$field1:$field2-$field3"
    
    # Create the output directory 
    vcf_output_directory="variants/$word"
    mkdir -p "$vcf_output_directory"

    # Construct the output file name
    vcf_output_file="$vcf_output_directory/$word.vcf"
    
    
    # Check if the output file already exists
    if [ -f "$vcf_output_file" ]; then
        echo "File $vcf_output_file already exists. Skipping extraction for $word."
        continue  # Skip to the next iteration
    fi

    # Print the bcftools command for debugging
    echo "Running command: bcftools view --threads 12 -r \"$arg\" \"$input_file\" -o \"$vcf_output_file\""

    # Run the bcftools command with the constructed argument
    if ! bcftools view --threads 12 -r "$arg" "$input_file" -Oz -o "$vcf_output_file" ; then
        echo "bcftools failed to extract the information for gene: $word, skipping"
        continue
    fi

    # Check if bcftools completed successfully
    if [ $? -eq 0 ]; then
        echo "Information for specified region was successfully extracted and saved to $vcf_output_file"
        else
        echo "Error: bcftools failed to extract information for $word"
        exit 1
    fi
    
    # Extract the VCF directory and file name
    vcf_dir=$(dirname "$vcf_output_file")
    vcf_file=$"$vcf_output_file"

    # Check if the VCF file exists
    if [ ! -f "$vcf_file" ]; then
    echo "VCF file $vcf_file does not exist."
      exit 1
    fi

   # Loop through all subdirectories (sample list folders) in the parent directory
   for sample_list_dir in "$sample_list_parent_dir"/*/; do
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
              output_file="$output_dir/$sample_name.vcf.gz"

              # Run the script with the VCF file, sample list file, and output file
              Scripts/vcf_by_sample.sh "$vcf_file" "$sample_list_file" "$output_file"

              # Message to indicate completion                          
              echo "Processed $sample_list_file"
          done
      fi
    done

  done <<< "$result"
done < "$gene_locations_file"




