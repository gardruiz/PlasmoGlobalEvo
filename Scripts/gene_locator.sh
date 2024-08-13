#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: $0 <gene_locations_file>"
  exit 1
fi

gene_locations_file="$1"
output_file="${gene_locations_file%.*}_locations.txt"

while read word; do
  grep -w "Name=$word" Genome/Genes/Pfalciparum_replace_Pf3D7_MIT_v3_with_Pf_M76611.gff | cut -f1,4,5 >> "$output_file"
done < "$gene_locations_file"

