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
results_folder="results/$gene_name"



mkdir -p "$results_folder"

#dN/dS function 
run_dnds() {
    file1="$1"
    file2="$2"

    # Creates unique id
    uuid=$(uuidgen)

    # Concatenate the fasta files
    concatenated_file="/tmp/align${uuid}.fasta"
    cat "$file1" "$file2" > "$concatenated_file"

    # Create YN00 input file name
    yn00_input_file="/tmp/aling${uuid}.phy"
    yn00_output_file="tmp/yn00_output_${uuid}.txt"

    # Convert the sequences to PHYLIP format
    python Scripts/fasta_to_phylip_sequential.py $concatenated_file $yn00_input_file
    

  

    # Create a unique control file for YN00
    control_file="Config/yn00_control_file_${uuid}.ctl"
    cat > "$control_file" <<EOL
seqfile = $yn00_input_file       * sequence data file name
outfile = $yn00_output_file      * main result file name

* Additional options (optional):
verbose = 1               * level of output detail (0: minimal, 1: detailed)
icode = 0                 * genetic code (0: universal code)
weighting = 0             * weighting for transition/transversion ratio estimation
commonf3x4 = 0            * whether to use common frequencies for 3x4 codon table
EOL

    # Run CODEML using the generated control file
    yn00 "$control_file"
  


    # Calculate dN/dS using grep and awk
    result=$(awk '/seq\. seq./,/LWL85, LPB93 & LWLm methods/ {if (!/LWL85, LPB93 & LWLm methods/) print}' "$yn00_output_file" | awk '{sumdN += $8; sumdS += $11} END {print sumdN/NR "\t" sumdS/NR "\t" (sumdN/NR)/(sumdS/NR)}')
    
    country1=$(basename "$1" | cut -d '_' -f 1)
    country2=$(basename "$2" | cut -d '_' -f 1)

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"

    # Make a directory to store the results 

    # Construct the output filename
    filename="$(basename "$file1" .fasta)-$(basename "$file2" .fasta).out"

    # Store codeml output
    mv "$yn00_output_file" "$results_folder/$filename"

    # Clean up temporary files
    rm  "$yn00_input_file" "$control_file" 

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
