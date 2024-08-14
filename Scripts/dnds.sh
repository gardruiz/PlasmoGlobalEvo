#!/bin/bash

if [ $# -ne 3 ]; then
        echo "Usage: $0 pop1_fasta_files_dir pop2_fast_files_dir output_file"
    exit 1
fi

pop1_directory="$1"
pop2_directory="$2"
output_file="$3"

pop1_files=("$pop1_directory"/*.fasta)
pop2_files=("$pop2_directory"/*.fasta)


# Create a folder to store results

gene_name=$(basename "$(dirname "$(dirname "$pop1_directory")")")
results_folder="Results/dnds/yn00/$gene_name"



mkdir -p "$results_folder"

#dN/dS function 
run_dnds() {
    file1="$1"
    file2="$2"

    # Check if tmp directory exists, create if not
    tmp_directory="tmp"
    if [ ! -d "$tmp_directory" ]; then
        mkdir "$tmp_directory"
    fi

    # Concatenate the fasta files
#    concatenated_file="tmp/concat.fasta"
    concatenated_file="tmp/align.fasta"
    cat "$file1" "$file2" > "$concatenated_file"

    # Convert the sequences to PHYLIP format
    python Scripts/fasta_to_phylip_sequential.py tmp/align.fasta tmp/concatenated.phy

    # Run CODEML using the control file, automatically press "Enter"
    yn00 Config/yn00_control_file.ctl

    # Calculate dN/dS using grep and awk
    result=$(awk '/seq\. seq./,/LWL85, LPB93 & LWLm methods/ {if (!/LWL85, LPB93 & LWLm methods/) print}' tmp/yn00_output.txt | awk '{sumdN += $8; sumdS += $11} END {print sumdN/NR "\t" sumdS/NR "\t" (sumdN/NR)/(sumdS/NR)}')
    
    country1=$(basename "$1" | cut -d '_' -f 1)
    country2=$(basename "$2" | cut -d '_' -f 1)

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"

    # Make a directory to store the results 

    # Construct the output filename
    filename="$(basename "$file1" .fasta)-$(basename "$file2" .fasta).out"

    # Store codeml output
    mv tmp/yn00_output.txt "$results_folder/$filename"

    # Clean up temporary files
    rm tmp/concat.fasta tmp/align.fasta  tmp/concatenated.phy 2ML.dN  2ML.dS  2ML.t 2NG.dN  2NG.dS  2NG.t   2YN.dN  2YN.dS  2YN.t rst rst1 rub

}

#Loop running dnds for files in population 1
for file1 in "${pop1_files[@]}"; do
    for file2 in "${pop1_files[@]}"; do
        run_dnds "$file1" "$file2"
    done
done

#Loop runnin dnds for files in population 2
for file1 in "${pop2_files[@]}"; do
    for file2 in "${pop2_files[@]}"; do
        run_dnds "$file1" "$file2"
        
    done
done

#Loop running dnds for files in population 1 vs files in population 2
for file1 in "${pop1_files[@]}"; do
    for file2 in "${pop2_files[@]}"; do
        run_dnds "$file1" "$file2"
    done
done
