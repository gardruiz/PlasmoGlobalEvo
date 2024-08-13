import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
# Directory containing CSV files
directory = sys.argv[1]
# Read the second file containing country information
country_df = pd.read_csv(sys.argv[2], sep='\t')
# Drop the 'snp_distance' column from the country DataFrame
country_df = country_df.drop(columns=['sequence'])
# List to store DataFrame from each CSV file
dfs = []
# Iterate over each CSV file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        filepath = os.path.join(directory, filename)
        # Read CSV file into a DataFrame
        df = pd.read_csv(filepath)
        # Append DataFrame to the list
        dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
combined_df = pd.concat(dfs)

# Calculate the number of entries to keep (top 10%)
top_entries_count = int(len(combined_df) * 0.1)

# Sort DataFrame by snp_distance and keep the top 10% entries
top_entries = combined_df.nsmallest(top_entries_count, 'snp_distance')

# Merge the DataFrames on 'sample' and 'sequence2' columns
result_df = pd.merge(top_entries, country_df, left_on=['sequence2'], right_on=['sample'], how='left')
result_df = result_df.drop(columns=['sample'])


# Plot the distribution of snp_distances
plt.figure(figsize=(12, 8))
sns.kdeplot(data=top_entries, x='snp_distance', fill=True, color='blue')
plt.xlabel('SNP Distance', fontsize=14)
plt.ylabel('Density', fontsize=14)
plt.title('Distribution of SNP Distances (Top 10% pairs)', fontsize=16)
plt.show()

# Save the filtered entries to a CSV file
result_df.to_csv('top_entries.csv', index=False)


