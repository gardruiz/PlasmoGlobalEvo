#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <folder> <output_file> <alpha>"
    exit 1
fi

# Specify the directory containing the files
folder="$1"

# Specify the output file
output_file="$2"

# Specify the Laplacian smoothing parameter
alpha="$3"

# Function to calculate probabilities after additive smoothing
calculate_probabilities() {
    local -n input_data=$1
    local smoothing_constant=$2

    total_count=0
    total_categories=${#input_data[@]}

    for category_count in "${input_data[@]}"; do
        total_count=$((total_count + category_count + smoothing_constant))
    done

    for category in "${!input_data[@]}"; do
        count=${input_data[$category]}
        smoothed_count=$((count + smoothing_constant))
        probability=$(bc -l <<< "scale=10; $smoothed_count / $total_count")
        probabilities[$category]=$probability
    done
}

# Loop over all files in the folder
for file in "$folder"/*; do
    # Extract country names from the filename
    filename=$(basename "$file")
    country1=$(echo "$filename" | cut -d '-' -f 1 | cut -d '_' -f 1)
    country2=$(echo "$filename" | cut -d '-' -f 2 | cut -d '_' -f 1)

    # Calculate dN/dS using grep and awk with Laplacian smoothing
    result=$(grep -i 'dN/dS=' "$file" | awk -v alpha="$alpha" '{
        sum11 += $11
        sum14 += $14
        if ($11 != 0) count11++
        if ($14 != 0) count14++
    }
    END {
        # Calculate probabilities after additive smoothing
        non_zero_count11 = (count11 == 0) ? 1 : count11
        non_zero_count14 = (count14 == 0) ? 1 : count14
        
        smoothed_sum11 = sum11 + (alpha * non_zero_count11)
        smoothed_sum14 = sum14 + (alpha * non_zero_count14)

        probability_sum11 = smoothed_sum11 / (NR + (alpha * non_zero_count11))
        probability_sum14 = smoothed_sum14 / (NR + (alpha * non_zero_count14))
        
        ratio = probability_sum11 / probability_sum14

        print probability_sum11 "\t" probability_sum14 "\t" ratio
    }')

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"
done

