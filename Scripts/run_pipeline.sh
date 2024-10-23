#!/bin/bash

# Function to read the control file and export variables
parse_control_file() {
    local control_file="Config/control_file.txt"
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


control_file="Config/control_file.txt"

# Parse the control file
parse_control_file "$control_file"

# Check if use_gene_name is set, if not default to using IDs
if [ -z "$use_gene_name" ]; then
    echo "use_gene_name is not set, defaulting to using IDs."
    use_gene_name="False"
fi

# Run scripts with arguments from the control file
run_command "python3 Scripts/multysampler.py $Populations_1_file $Number_of_samples $Samples_ID_file $Years_file $Samples_1_output_dir"
run_command "python3 Scripts/multysampler.py $Populations_2_file $Number_of_samples $Samples_ID_file $Years_file $Samples_2_output_dir"
run_command "Scripts/vcf_cutter.sh $Gene_names_file $Samples_file_dir $VCF_file $GFF_file $use_gene_name"

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
            # Check if any of the fields are empty (invalid line)
            if [ -z "$gene_name" ] || [ -z "$chrom" ] || [ -z "$gene_start" ] || [ -z "$gene_end" ]; then
                echo "Warning: Skipping invalid line (missing fields) - $gene_name $chrom $gene_start $gene_end"
                continue  # Skip to the next line
            fi

            echo "Processing gene: $gene_name"

		# Identify subdirectories inside the gene directory
                gene_dir="variants/$gene_name"
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

                output_gene_file="results/${gene_name}_dxy.txt"  # Create gene-specific output file
                
		# Run the dxy script and check for success
                if ! run_command "Scripts/dxy.sh $pop1_files $pop2_files $output_gene_file $gene_start $gene_end"; then
                    echo "Error: Failed to calculate dxy for $gene_name! Skipping this gene."
                    # Clean up temporary files
                    rm "$pop1_files" "$pop2_files"
                    continue  # Skip to the next gene
                fi
		
		# Plot data
                if ! run_command "python3 Scripts/heatmap2.py results/${gene_name}_dxy.txt"; then
                    echo "Error: Failed to plot smoothed heatmap for $gene_name! Skipping this gene."
                    continue  # Skip to the next gene
                fi


		# Clean up temporary files
    		rm "$pop1_files" "$pop2_files"

            done < "$gene_locations_file"

            ;;

	dnds)
            
            run_command() {
		command="$1"
    		echo "Running command: $command"
    		$command
    		if [ $? -ne 0 ]; then
        	    echo "Error: Command failed - $command"
                    exit 1
    		fi
	    }

	
	    echo "Calculating dN/dS"

            # Function to calculate dN/dS for each gene
	    calculate_dnds() {
		gene_name="$1"
                echo "Processing gene: $gene_name"

                # Identify subdirectories inside the gene directory
                gene_dir="variants/$gene_name"
                subdirs=($(find "$gene_dir" -mindepth 1 -maxdepth 1 -type d))

                if [ "${#subdirs[@]}" -lt 2 ]; then
                    echo "Error: Less than 2 subdirectories found in $gene_dir!"
                    exit 1
                fi

                # Dynamically assign Pop1_files and Pop2_files from the first two subdirectories
                Pop1_dir="${subdirs[0]}"
                Pop2_dir="${subdirs[1]}"

		# Check if the file exists
	

                output_dnds_file="results/${gene_name}_dnds.txt"  # Create gene-specific output file
		if [ -f "$output_dnds_file" ]; then
		    echo "File $output_dnds_file exists, replacing it."
		    rm "$output_dnds_file"  # Remove the existing file
		fi

		# Extract the ID from the GFF file
		if [ "$use_gene_name" = "True" ]; then
		    echo "Using gene names."
		    gene_info=$(grep "Name=$gene_name;" "$GFF_file")
		    # Extract the ID from the gene_info using awk
                    gene_ID=$(echo "$gene_info" | awk -F 'ID=' '{print $2}' | awk -F ';' '{print $1}')
                    echo "Gene ID: $gene_ID"

	        else 
		    echo "Using IDs names."
                    gene_info=$(grep "ID=$gene_name;" "$GFF_file")
		    gene_ID=$gene_name
		fi

		# Extract the genomic region (chromosome, start, and end) using awk
		Genomic_region=$(echo "$gene_info" | awk '{print $1":"$4"-"$5}')

		tmp_dir="/tmp/${gene_name}_$$"
                mkdir -p "$tmp_dir"
                
	        if ! run_command "Scripts/bed_extractor.sh $gene_ID $GFF_file $tmp_dir"; then
		    echo "Error: Failed to extract BED file for $gene_name! Skipping this gene."
                    rm -rf "$tmp_dir"  # Clean up temp directory
                    return  # Skip to the next gene
                fi	
                
		BED_file="$tmp_dir/$gene_ID.bed"

		# Extract CDS from populations dir
    		if ! run_command "Scripts/vcf_to_cds.sh $Pop1_dir $Genomic_region $BED_file"; then
        	    echo "Error: Failed to extract CDS for population 1 of $gene_name! Skipping this gene."
                    rm -rf "$tmp_dir"  # Clean up temp directory
                    return  # Skip to the next gene
                fi

    		if ! run_command "Scripts/vcf_to_cds.sh $Pop2_dir $Genomic_region $BED_file"; then
        	    echo "Error: Failed to extract CDS for population 2 of $gene_name! Skipping this gene."
        	    rm -rf "$tmp_dir"  # Clean up temp directory
		    return  # Skip to the next gene
    		fi

    		# Calculate dnds
    		echo "Calculating dN/dS"
    		if ! run_command "Scripts/dnds.sh $Pop1_dir/sequences $Pop2_dir/sequences $output_dnds_file"; then
        	    echo "Error: Failed to calculate dN/dS for $gene_name! Skipping this gene."
        	    rm -rf "$tmp_dir"  # Clean up temp directory
                    return  # Skip to the next gene
    		fi
    
    		# Clean up temporary files
    		rm "$BED_file"

    		# Plot heatmap
    		if ! run_command "python3 Scripts/heatmap.py $output_dnds_file"; then
        	    echo "Error: Failed to plot heatmap for $gene_name! Skipping this gene."
                    rm -rf "$tmp_dir"  # Clean up temp directory
                    return  # Skip to the next gene
    		fi

    		# Check if smoothing is enabled in the control file
    		if [ "$Smoothing" = "True" ]; then
       		    echo "Smoothing enabled. Calculating smoothed dnds estimates for $gene_name."
        
                    # Clean the PAML/YN00 output files
        	    if ! run_command "Scripts/clean_yn00_output.sh results/$gene_name"; then
            		echo "Error: Failed to clean YN00 output for $gene_name! Skipping this gene."
            		rm -rf "$tmp_dir"  # Clean up temp directory
            		return  # Skip to the next gene
        	    fi

        	    # Calculate smoothed estimates
        	    if ! run_command "python3 Scripts/smoothing.py results/$gene_name/clean results/${gene_name}_smoothed_dnds.txt"; then
            		echo "Error: Failed to calculate smoothed dN/dS for $gene_name! Skipping this gene."
            		rm -rf "$tmp_dir"  # Clean up temp directory
        		return  # Skip to the next gene
        	     fi

        	     # Plot data
        	     if ! run_command "python3 Scripts/heatmap.py results/${gene_name}_smoothed_dnds.txt"; then
                         echo "Error: Failed to plot smoothed heatmap for $gene_name! Skipping this gene."
            		 rm -rf "$tmp_dir"  # Clean up temp directory
            		 return  # Skip to the next gene
        	     fi
    	        fi

    		# Clean up temporary directory
    		rm -rf "$tmp_dir"
	    }

            
	    # Export dnds function
	    export -f calculate_dnds
		
	    export -f run_command

	    # Export necessary variables to be accessible within the parallel jobs
    	    export gene_ID GFF_file Smoothing gene_locations_file

            # Use GNU parallel to process each gene concurrently
            cat "$gene_locations_file" | awk '{print $1}' | parallel --line-buffer -j $(nproc) calculate_dnds
            ;;
    esac
done

if [ "$Use_average_ds" = "True" ]; then
    echo "Calculating dN/dS using average dS"

    # Calculate dNdS using the average dS
    run_command "python3 Scripts/get_dnds_using_average_ds.py $Gene_names_file"
fi

echo "Pipeline executed successfully."
