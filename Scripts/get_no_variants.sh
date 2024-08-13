
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

# Loop over all files in the folder
for file in "$folder"/*; do
    # Extract country names from the filename
    filename=$(basename "$file")
    country1=$(echo "$filename" | cut -d '-' -f 1 | cut -d '_' -f 1)
    country2=$(echo "$filename" | cut -d '-' -f 2 | cut -d '_' -f 1)

    # Calculate dN/dS using grep and awk with Laplacian smoothing
    result=$(grep -i 'dN/dS=' "$file" | awk -v alpha="$alpha" '{
        sum11 += $11 + alpha
        sum14 += $14 + alpha
        if ($11 != 0) count11++
        if ($14 != 0) count14++
    }
    END {
        # Calculate non-zero counts
        non_zero_count11 = (count11 == 0) ? 1 : count11
        non_zero_count14 = (count14 == 0) ? 1 : count14
        
        # Calculate Laplacian smoothed ratios
        ratio = (sum11 / (NR + (alpha * non_zero_count11))) / (sum14 / (NR + (alpha * non_zero_count14)))
        
        print sum11 / (NR + (alpha * non_zero_count11)) "\t" sum14 / (NR + (alpha * non_zero_count14)) "\t" ratio
    }')

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"
done

