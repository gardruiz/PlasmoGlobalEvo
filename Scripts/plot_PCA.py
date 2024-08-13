
import sys
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle

if len(sys.argv) != 3:
    print("Usage: python plot_PCA.py input_file country_table")
    exit(1)

# Read data from a tab-delimited file
file_path = sys.argv[1]
df = pd.read_csv(file_path, delimiter='\t')
print(df.head())

# Read country information
country_table_path = sys.argv[2]
country_df = pd.read_csv(country_table_path, delimiter='\t')

# Merge dataframes on IID
merged_df = pd.merge(df, country_df[['Sample', 'Country']], left_on='#FID', right_on='Sample')

# Confirm data types
print("Data Types:")
print(merged_df.dtypes)

# Attempt to convert 'PC1' and 'PC2' to numeric, replace non-numeric values with NaN
merged_df['PC1'] = pd.to_numeric(merged_df['PC1'], errors='coerce')
merged_df['PC2'] = pd.to_numeric(merged_df['PC2'], errors='coerce')

# Get a unique list of countries and create a cycle of colors
unique_countries = merged_df['Country'].unique()
color_cycle = cycle(plt.cm.get_cmap('Set2', len(unique_countries)).colors)

# Plotting the scatter plot with individual colors for each country
for country in unique_countries:
    country_data = merged_df[merged_df['Country'] == country]
    color = next(color_cycle)
    plt.scatter(country_data["PC1"], country_data["PC2"], marker='o', c=color, label=country)

plt.title('PCA Plot')
plt.xlabel('PC1')
plt.ylabel('PC2')

# Move the legend to the right outside of the plot
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()



