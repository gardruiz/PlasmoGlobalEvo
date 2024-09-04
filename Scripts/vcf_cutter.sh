#!/bin/bash

if [ $# -ne 4 ]; then
  echo "Usage: $0 <gene_names_file>  <sample_list_directory> <vcf_file> <GFF_file>" 
  exit 1
fi

gene_locations_file="$1"
sample_list_parent_dir="$2"
input_file="$3"
gff_file="$4"
gene_locations_output_file="${gene_locations_file%.*}_locations.txt"
gene_locations_output_directory="Genome/genes"

# Create the output directory if it doesn't exist
mkdir -p "$gene_locations_output_directory"

while read word; do
  result=$(grep -w "Name=$word" "$gff_file" | cut -f1,4,5)
  
  # Append the $word as the first term in each line of the result
  while read line; do
    echo "$word $line" >> "$gene_locations_output_file"
  done <<< "$result"

  # Loop through each line in the genes file
  while IFS=$'\t' read -r field1 field2 field3; do
    # Construct the argument using the fields
    arg="$field1:$field2-$field3"
    
    # Create the output directory 
    vcf_output_directory="Variants/$word"
    mkdir -p "$vcf_output_directory"

    # Construct the output file name
    vcf_output_file="$vcf_output_directory/$word.vcf"
  
    # Print the bcftools command for debugging
    echo "Running command: bcftools view --threads 8 -r \"$arg\" \"$input_file\" -o \"$vcf_output_file\""

    # Run your command with the constructed argument
    bcftools view  --threads 8  -r "$arg"  $input_file -o $vcf_output_file

    # Check if vcftools completed successfully
    if [ $? -eq 0 ]; then
	    echo "Information for specified region was successfully extracted and saved to $vcf_output_file."
    else
	    echo "Error: vcftools failed to extract information for $word."
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
              output_file="$output_dir/$sample_name.vcf"

              # Run the script with the VCF file, sample list file, and output file
              Scripts/vcf_by_sample.sh "$vcf_file" "$sample_list_file" "$output_file"

              # Message to indicate completion                          
              echo "Processed $sample_list_file"
          done
      fi
    done

  done <<< "$result"
done < "$gene_locations_file"




