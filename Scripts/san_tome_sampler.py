import pandas as pd
import sys


# Check for the right number of arguments and print usage
if len(sys.argv) != 4:
    print ("Usage: python script.py sample_file gene sample_size")
    exit(1)

# Asing arguments to variables
gene = sys.argv[2] 
n = int(sys.argv[3])

# Load your CSV file into a Pandas DataFrame
df = pd.read_csv(sys.argv[1])

# Get unique values in the "Year" column
unique_years = df["Year"].unique()
print(unique_years)

# Iterate through unique years
for year in unique_years:
    year_condition = (df["Year"] == year)
    
    # Iterate through unique gene values within each year
    for gene in df[year_condition]["Gene"].str.upper():
        gene_condition = (df["Gene"].str.upper == gene)
        
        # Apply both year and gene conditions
        filtered_df = df[year_condition & gene_condition]
        
        # Take a random sample of size 'n' from the "asv" column
        random_sample = filtered_df["asv"].sample(n)
        
        # Construct the file name based on Gene and Year
        file_name = f"{gene.upper}_{year}.csv"
        
        # Save the random sample to a CSV file
        random_sample.to_csv(file_name, header=True, index=False)

