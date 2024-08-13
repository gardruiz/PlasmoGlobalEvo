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

# Initialize variables to calculate total dN, dS, and number of records
total_dN=0
total_dS=0
total_NR=0

# Loop over all files to calculate total dN, dS, and number of records
for file in "$folder"/*; do
    # Debug: Print current file being processed
    echo "Processing file: $file"

    # Print header for debugging
    echo "Printing values of dN (column 8) and dS (column 11):"
    echo "---------------------------------------------------"

    # Calculate total dN, dS, and number of records for the current file
    values=$(awk 'NR>1 && NF>1 { 
        print $8, $11; 
        total_dN += $8; 
        total_dS += $11; 
        total_NR++ 
    } END { 
        print total_dN, total_dS, total_NR 
    }' "$file")

    # Debug: Print the values captured from awk
    echo "Values captured: $values"

    # Capture and update total_dN, total_dS, and total_NR from awk output
    total_dN=$(echo "$values" | awk '{print $1}')
    total_dS=$(echo "$values" | awk '{print $2}')
    total_NR=$(echo "$values" | awk '{print $3}')

    # Debug: Print total dN, dS, and number of records after processing current file
    echo "Total dN so far: $total_dN"
    echo "Total dS so far: $total_dS"
    echo "Total NR so far: $total_NR"
done

# Calculate average dN and dS values
if [ "$total_NR" -gt 0 ]; then
    # Calculate average dN and dS values
    avg_dN=$(echo "scale=6; $total_dN / $total_NR" | bc)
    avg_dS=$(echo "scale=6; $total_dS / $total_NR" | bc)

    # Debug: Print average dN and dS values
    echo "Average dN: $avg_dN"
    echo "Average dS: $avg_dS"
else
    echo "No records found to calculate averages."
    exit 1
fi

# Loop over all files to calculate smoothed dN/dS estimates
for file in "$folder"/*; do
    # Extract country names from the filename
    filename=$(basename "$file")
    country1=$(echo "$filename" | cut -d '-' -f 1 | cut -d '_' -f 1)
    country2=$(echo "$filename" | cut -d '-' -f 2 | cut -d '_' -f 1)

    # Calculate dN/dS using grep and awk with average smoothing
    result=$(awk 'NR>1 && NF>1 { sumdN += $8 + '"$avg_dN"'; sumdS += $11 + '"$avg_dS"' }    
    END {
        print sumdN/NR "\t" sumdS/NR "\t" (sumdN/(NR))/(sumdS/(NR))                 
    }' "$file")

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"
done

