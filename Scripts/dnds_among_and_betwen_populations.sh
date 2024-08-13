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
results_folder="Results/dnds/codeml/$gene_name"



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
#    concatenated_file="tmp/concatenated.fasta"
    concatenated_file="tmp/alignment.fasta"
    cat "$file1" "$file2" > "$concatenated_file"

    # Align the concatenated file using Muscle
#    muscle_output="tmp/alignment.fasta"
#    muscle -super5 "$concatenated_file" -output "$muscle_output"


    # Run CODEML using the control file, automatically press "Enter"
    codeml Config/codeml_control_file.ctl

    # Calculate dN/dS using grep and awk
    result=$(grep -i 'dN/dS=' "tmp/codeml.out" | awk '{sum11 += $11; sum14 += $14} END {print sum11/NR "\t" sum14/NR "\t" (sum11/NR)/(sum14/NR)}')
    
    country1=$(basename "$1" | cut -d '_' -f 1)
    country2=$(basename "$2" | cut -d '_' -f 1)

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"

    # Make a directory to store the results 

    # Construct the output filename
    filename="$(basename "$file1" .fasta)-$(basename "$file2" .fasta).out"

    # Store codeml output
    mv tmp/codeml.out "$results_folder/$filename"

    # Clean up temporary files
    rm tmp/concatenated.fasta tmp/alignment.fasta rst  rst1  rub 2ML.dN  2ML.dS  2ML.t  2NG.dN  2NG.dS  2NG.t  4fold.nuc

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




