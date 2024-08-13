
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
     # Calculate dN/dS using grep and awk
    result=$(awk '/seq\. seq./,/LWL85, LPB93 & LWLm methods/ {if (!/LWL85, LPB93 & LWLm methods/) print}' $file | awk -v alpha="$alpha" '{
    sumdN += $8 + alpha
    sumdS += $11 + alpha
} 
END {
print sumdN/(NR - 2) "\t" sumdS/(NR - 2) "\t" (sumdN/(NR - 2))/(sumdS/(NR - 2))
}')

    # Append results to the output file
    echo -e "$country1\t$country2\t$result" >> "$output_file"
done


