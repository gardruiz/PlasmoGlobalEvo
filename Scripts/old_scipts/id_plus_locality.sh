#!/bin/bash

input_file="$1"
output_file="$2"
text_to_add="$3"

if [ $# -ne 3 ]; then
    echo "Usage: $0 input_file output_file text_to_add"
    exit 1
fi

awk -v text="$text_to_add" -F'\t' -v OFS='\t' '{print $0, text}' "$input_file" > "$output_file"

echo "Second column added. Result saved in $output_file."

