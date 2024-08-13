#!/bin/bash

# Check for the correct number of arguments
if [ $# -ne 2 ]; then
    echo "This script calculates Dxy for 2 populations. Usage: $0 vcf_file sample_id_file.txt"
    exit 1
fi

# Assign arguments to variables
vcf_file="$1"
sampleid_file="$2"

# Run the pixy command
pixy --stats pi fst dxy \
--vcf "$vcf_file" \
--populations "$sampleid_file" \
--window_size 10000 \
--n_cores 8 \
--output_folder output \
--output_prefix pixy_output \
--bypass_invariant_check 'yes'\
