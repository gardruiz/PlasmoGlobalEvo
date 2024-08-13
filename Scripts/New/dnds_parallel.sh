#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: $0 pop1_fasta_files_dir pop2_fast_files_dir output_file num_parallel_jobs"
    exit 1
fi

pop1_directory="$1"
pop2_directory="$2"
output_file="$3"
num_parallel_jobs="$4"

pop1_files=("$pop1_directory"/*.fasta)
pop2_files=("$pop2_directory"/*.fasta)

# dN/dS function
run_dnds() {
    file1="$1"
    file2="$2"

    # Create a unique temporary directory for each run
    tmp_directory=  "/media/Daniel/ADATA HM800/Malaria/tmp/$(date +%s)"
    mkdir -p "$tmp_directory"  # Use -p to avoid errors if directory already exists

    # Concatenate the fasta files
    concatenated_file="$tmp_directory/concatenated.fasta"
    cat "$file1" "$file2" > "$concatenated_file"

#    # Align the concatenated file using Muscle
 #   muscle_output="$tmp_directory/alignment.fasta"
#    muscle -super5 "$concatenated_file" -output "$muscle_output"

    # Generate a dynamic control file for CODEML
    dynamic_control_file="$tmp_directory/dynamic_codeml_control_file.ctl"
    cp /Config/codeml_dynamic_control_file.ctl "$dynamic_control_file"
    sed -i "s|<seqfile>|$concatenated_file|" "$dynamic_control_file"
    sed -i "s|<outfile>|$tmp_directory/codeml.out|" "$dynamic_control_file"

     # Run CODEML using the dynamically created control file
    cd "$tmp_directory" || exit 1  # Change to the temp directory
    codeml "$dynamic_control_file"
    cd - || exit 1  # Change back to the original directory

    # Calculate dN/dS using grep and awk
    result=$(grep -i 'dN/dS=' "$tmp_directory/codeml.out" | awk '{sum11 += $11; sum14 += $14} END {print (sum11/NR)/(sum14/NR)}')
    country1=$(basename "$1" | cut -d '_' -f 1)
    country2=$(basename "$2" | cut -d '_' -f 1)

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"

    # Clean up temporary directory
    rm -r "$tmp_directory"
}


export -f run_dnds  # Export the function for parallel execution

# Run dN/dS in parallel for files in population 1
parallel -j "$num_parallel_jobs" run_dnds ::: "${pop1_files[@]}" ::: "${pop1_files[@]}"

# Run dN/dS in parallel for files in population 2
parallel -j "$num_parallel_jobs" run_dnds ::: "${pop2_files[@]}" ::: "${pop2_files[@]}"

# Run dN/dS in parallel for files in population 1 vs files in population 2
parallel -j "$num_parallel_jobs" run_dnds ::: "${pop1_files[@]}" ::: "${pop2_files[@]}"

