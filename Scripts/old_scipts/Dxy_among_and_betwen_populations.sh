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

run_dxy() {
    result=$(python Scripts/Dxy.py "$1" "$2" "$gene_start" "$gene_end") 
    country1=$(basename "$1" | cut -d '_' -f 1)
    country2=$(basename "$2" | cut -d '_' -f 1)
    echo -e "$country1\t$country2\t$result" >> "$output_file"
}

for file1 in "${pop1_files[@]}"; do
    for file2 in "${pop1_files[@]}"; do
        if [ "$file1" != "$file2" ]; then
            run_dxy "$file1" "$file2"
        fi
    done
done

for file1 in "${pop2_files[@]}"; do
    for file2 in "${pop2_files[@]}"; do
        if [ "$file1" != "$file2" ]; then
            run_dxy "$file1" "$file2"
        fi
    done
done

for file1 in "${pop1_files[@]}"; do
    for file2 in "${pop2_files[@]}"; do
        run_dxy "$file1" "$file2"
    done
done

