#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <folder> <output_file>"
    exit 1
fi

# Specify the directory containing the files
folder="$1"

# Specify the output file
output_file="$2"

# Loop over all files in the folder
for file in "$folder"/*; do
    # Extract country names from the filename
    filename=$(basename "$file")
    country1=$(echo "$filename" | cut -d '-' -f 1 | cut -d '_' -f 1)
    country2=$(echo "$filename" | cut -d '-' -f 2 | cut -d '_' -f 1)

    echo "Processing file: $file"
    echo "Country 1: $country1"
    echo "Country 2: $country2"


    # Calculate dN/dS starting from the third line and skip lines where column 7 == 99.0000
    result=$(awk 'NR>2 && $7 != 99.0000 {
        sumdN += $8;
        sumdS += $11;
        count++;
    } 
    END {
        if (count > 0) {
            print sumdN/count "\t" sumdS/count "\t" (sumdN/count)/(sumdS/count);
        } else {
            print "0\t0\t0";
        }
    }' "$file")

    echo "Result: $result"

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"
done


