import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print("Usage: ./heatmap.py Dxy_file Gene_name")
    sys.exit(1)

file = sys.argv[1]

# Load data from TSV file
data = pd.read_csv(file, sep='\t', header=None, names=['Country1', 'Country2', 'GeneticVariation'])

# Pivot data
pivot_data = data.pivot(index='Country1', columns='Country2', values='GeneticVariation')
pivot_data_filled = pivot_data.fillna(pivot_data.T)

# Custom sorting orders for Asia and Africa
asia_countries = ['Cambodia', 'Laos', 'Myanmar', 'Thailand', 'Vietnam']
west_africa_countries = ['Ghana', 'Mauritania', 'Senegal', 'Gabon', 'Cameroon', 'Guinea']
east_africa_countries = ['Kenya', 'Tanzania']

# Reorder rows and columns based on the custom sorting orders
data = pivot_data_filled.loc[asia_countries + east_africa_countries + west_africa_countries,
                             asia_countries + east_africa_countries + west_africa_countries]

# Create a clustermap using Seaborn
plt.figure(figsize=(10, 8))
cluster = sns.clustermap(data, annot=True, cmap='YlGnBu', linewidths=.5, col_cluster=False, row_cluster=False,
                         row_colors=[sns.color_palette("Set2", n_colors=len(asia_countries + east_africa_countries + west_africa_countries))],
                         col_colors=[sns.color_palette("Set2", n_colors=len(asia_countries + east_africa_countries + west_africa_countries))])

# Add labels for the groups
plt.text(0, -1, 'Asia', ha='center', va='center', fontsize=12)
plt.text(len(asia_countries) + len(east_africa_countries) - 1, -1, 'East Africa', ha='center', va='center', fontsize=12)
plt.text(len(asia_countries) + len(east_africa_countries) + len(west_africa_countries) - 1, -1, 'West Africa', ha='center', va='center', fontsize=12)

plt.title(sys.argv[2] + ' Genetic Variation')
plt.xlabel('')
plt.ylabel('')
plt.show()

