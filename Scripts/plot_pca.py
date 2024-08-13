import sys
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle

# Define a dictionary mapping countries to continents
continent_mapping = {
    'South_America': ['Peru', 'Colombia', 'Venezuela'],
    'South_Asia': ['Bangladesh', 'India'],
    'West_Africa': ['Ghana', 'Mali', 'Gambia', 'Benin', 'Guinea', 'Senegal', 'Nigeria', 'Mauritania', "Côte d'Ivoire", 'Burkina Faso'],
    'Central_Africa': ['Cameroon','Gabon', 'Republic of the Congo'],
    'Southeast_Asia': ['Vietnam', 'Cambodia', 'Myanmar', 'Thailand', 'Laos'],
    'East_Africa': ['Kenya', 'Tanzania', 'Malawi', 'Sudan', 'Mozambique', 'Ethiopia', 'Madagascar', 'Uganda'],
    'Oceania': ['Papua New Guinea', 'Indonesia'],
    'Resistant' : ['Resistant'],
    'Sensitive' : ['Sensitive'],
    'Undetermined': ['Undertermined']
}

# Flatten the dictionary for a reverse mapping (country to continent)
country_to_continent = {country: continent for continent, countries in continent_mapping.items() for country in countries}

if len(sys.argv) != 4:
    print("Usage: python plot_PCA.py input_file metadata column")
    exit(1)

# Read data from a tab-delimited file
file_path = sys.argv[1]
df = pd.read_csv(file_path, delimiter='\t')
print(df.head())

# Read country information
country_table_path = sys.argv[2]
country_df = pd.read_csv(country_table_path, delimiter='\t')
column = sys.argv[3]

# Merge dataframes on IID
merged_df = pd.merge(df, country_df[['Sample', column]], left_on='#FID', right_on='Sample')

# Confirm data types
print("Data Types:")
print(merged_df.dtypes)

# Attempt to convert 'PC1' and 'PC2' to numeric, replace non-numeric values with NaN
merged_df['PC1'] = pd.to_numeric(merged_df['PC1'], errors='coerce')
merged_df['PC2'] = pd.to_numeric(merged_df['PC2'], errors='coerce')

# Drop rows with NaN values in 'PC1' or 'PC2'
merged_df = merged_df.dropna(subset=['PC1', 'PC2'])

# Map countries to continents using the dictionary
merged_df['Continent'] = merged_df[column].map(country_to_continent)

# Get a unique list of continents and create a cycle of colors
unique_continents = merged_df['Continent'].dropna().unique()
color_cycle = cycle(plt.cm.get_cmap('viridis', len(unique_continents)).colors)

# Plotting the scatter plot with individual colors for each continent
for continent in unique_continents:
    continent_data = merged_df[merged_df['Continent'] == continent]
    color = next(color_cycle)
    plt.scatter(continent_data["PC1"], continent_data["PC2"], marker='o', c=color, label=continent, s=10)

plt.xlabel('PC1', fontsize=14)
plt.ylabel('PC2', fontsize=14)

# Move the legend to the right outside of the plot
plt.legend(loc='upper right', frameon=False)
plt.tight_layout()
plt.show()
