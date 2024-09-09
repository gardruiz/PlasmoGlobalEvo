#!/bin/bash

# Function to read the control file and export variables
parse_control_file() {
    local control_file=$1
    while IFS='=' read -r key value; do
        # Strip leading/trailing whitespace and export variable
        key=$(echo "$key" | xargs)
        value=$(echo "$value" | xargs)
        export "$key=$value"
    done < "$control_file"
}

# Function to run a command and handle errors
run_command() {
    local command="$1"
    echo "Running command: $command"
    eval "$command"
    if [ $? -ne 0 ]; then
        echo "Error: Command failed: $command"
        exit 1
    fi
}

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 control_file.txt"
    exit 1
fi

control_file="$1"

# Parse the control file
parse_control_file "$control_file"

# Run scripts with arguments from the control file
run_command "python3 Scripts/New/multysampler.py $Populations_1_file $Number_of_samples $Samples_ID_file $Years_file $Samples_1_output_dir"
run_command "python3 Scripts/New/multysampler.py $Populations_2_file $Number_of_samples $Samples_ID_file $Years_file $Samples_2_output_dir"
run_command "Scripts/New/vcf_cutter.sh $Gene_names_file $Samples_file_dir $VCF_file $GFF_file"

# Check which pipelines to run
IFS=',' read -r -a pipelines <<< "$Pipelines"

# Parse gene locations
gene_locations_file="${Gene_names_file%.*}_locations.txt"

for pipeline in "${pipelines[@]}"; do
    case "$pipeline" in
        dxy)
            echo "Calculating nucleotide pairwise diversity"

            # Loop through the gene list
            while IFS=$'\t' read -r gene_name chrom gene_start gene_end; do
                echo "Processing gene: $gene_name"

I		# Identify subdirectories inside the gene directory
                gene_dir="Variants/$gene_name"
                subdirs=($(find "$gene_dir" -mindepth 1 -maxdepth 1 -type d))

                if [ "${#subdirs[@]}" -lt 2 ]; then
                    echo "Error: Less than 2 subdirectories found in $gene_dir!"
                    exit 1
                fi

                # Dynamically assign Pop1_files and Pop2_files from the first two subdirectories
                Pop1_dir="${subdirs[0]}"
                Pop2_dir="${subdirs[1]}"

		# Generate temporary files for population lists
    		pop1_files=$(mktemp)
    		pop2_files=$(mktemp)

                # Get the list of files from each subdirectory
		find "$Pop1_dir" -type f > "$pop1_files"
		find "$Pop2_dir" -type f > "$pop2_files"

                output_gene_file="Results/$gene_name.dxy.txt"  # Create gene-specific output file
                
		run_command "Scripts/dxy.sh $pop1_files $pop2_files $output_gene_file $gene_start $gene_end"
		# Clean up temporary files
    		rm "$pop1_files" "$pop2_files"

            done < "$gene_locations_file"

            ;;
        *)
            echo "Warning: Unknown pipeline '$pipeline'. Skipping."
            ;;
    esac
done

echo "Pipeline executed successfully."

