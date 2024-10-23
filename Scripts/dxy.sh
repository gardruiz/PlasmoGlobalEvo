#!/bin/bash

if [ $# -ne 5 ]; then
    echo "Usage: $0 pop1_files.txt pop2_files.txt output_file gene_start gene_end"
    exit 1
fi

pop1_files_file="$1"
pop2_files_file="$2"
output_file="$3"
gene_start="$4"
gene_end="$5"

pop1_files=($(cat "$pop1_files_file"))
pop2_files=($(cat "$pop2_files_file"))

# Remove existing output file if it exists
if [ -e "$output_file" ]; then
   rm "$output_file"
   echo "Removed existing $output_file"
fi

run_dxy() {
    result=$(python Scripts/Dxy.py "$1" "$2" "$gene_start" "$gene_end") 
    
    # Check if result is None
    if [ "$result" == "None" ]; then
        echo "Error: Dxy.py returned None for $1 and $2. Stopping execution."
        return 1  # Indicate failure
    fi

    # Extract filenames without any extensions and parent directory
    country1=$(basename "$1" | sed 's/\.[^.]*$//' | sed 's/\.[^.]*$//' | sed 's/.*\///')
    country2=$(basename "$2" | sed 's/\.[^.]*$//' | sed 's/\.[^.]*$//' | sed 's/.*\///')
    echo -e "$country1\t$country2\t$result" >> "$output_file"
}

# Loop running dN/dS for unique pairs in population 1

process_all_files() {
    for ((i = 0; i < ${#pop1_files[@]}; i++)); do
            run_dxy "${pop1_files[i]}" "${pop1_files[i]}"
	    for ((j = i + 1; j < ${#pop1_files[@]}; j++)); do
            run_dxy "${pop1_files[i]}" "${pop1_files[j]}"
        done
    done

    # Loop running dN/dS for unique pairs in population 2
    for ((i = 0; i < ${#pop2_files[@]}; i++)); do
            run_dxy "${pop2_files[i]}" "${pop2_files[i]}"
            for ((j = i + 1; j < ${#pop2_files[@]}; j++)); do
            run_dxy "${pop2_files[i]}" "${pop2_files[j]}"
        done
    done

    # Loop running dN/dS for all pairs between population 1 and population 2
    for file1 in "${pop1_files[@]}"; do
        for file2 in "${pop2_files[@]}"; do
            run_dxy "$file1" "$file2"
        done
    done

    return 0  # Indicate success if all loops pass
}

# Call the function to process all files
process_all_files
if [ $? -ne 0 ]; then
    echo "Error occurred while processing files. Stopping execution."
fi

