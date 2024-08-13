

import pandas as pd
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster import hierarchy

if len(sys.argv) != 3:
    print("Usage: ./heatmap.py Dxy_file Gene_name")
    sys.exit(1)


file = sys.argv[1]

# Load data from TSV file
data = pd.read_csv(file, sep='\t', header=None, names=['Country1', 'Country2', 'dN', 'dS', 'dN/dS',])
data.replace([np.inf, -np.inf], np.nan, inplace=True)
data.dropna()

# Drop duplicate entries with swapped order of countries
data['Country1'], data['Country2'] = np.sort(data[['Country1', 'Country2']], axis=1).T
print(data.columns)
print(data.head)

data.drop_duplicates(subset=['Country1', 'Country2'], keep='first', inplace=True)

pivot_dn = data.pivot(index='Country1', columns='Country2', values='dN')
pivot_ds = data.pivot(index='Country1', columns='Country2', values='dS')
pivot_dnds = data.pivot(index='Country1', columns='Country2', values='dN/dS')
dn_data = pivot_dn.fillna(pivot_dn.T)
ds_data = pivot_ds.fillna(pivot_ds.T)
dnds_data = pivot_dnds.fillna(pivot_dnds.T)

print('dN values')
print(dn_data.describe())
print(dn_data)
print('dS values')
print(ds_data.describe())
print(ds_data)
print('dN/dS values')
print(dnds_data.describe())
print(dnds_data)

# Create custom sorting orders for Asia and Africa
asia_countries = ['Cambodia', 'Thailand', 'Myanmar', 'Laos', 'Vietnam']
west_africa_countries = ['Ghana', 'Mali', 'Gambia', 'Mauritania', 'Senegal', 'Gabon', 'Cameroon', 'Guinea']
east_africa_countries =['Kenya', 'Tanzania']


# Remove countries that are not present in the data
asia_countries = [country for country in asia_countries if country in dn_data.index]
east_africa_countries = [country for country in east_africa_countries if country in dn_data.index]
west_africa_countries = [country for country in west_africa_countries if country in dn_data.index]

# Reorder rows and columns based on the custom sorting orders
dn = dn_data.loc[asia_countries + east_africa_countries + west_africa_countries, asia_countries + east_africa_countries + west_africa_countries ]
ds = ds_data.loc[asia_countries + east_africa_countries + west_africa_countries, asia_countries + east_africa_countries + west_africa_countries ]
dnds = dnds_data.loc[asia_countries + east_africa_countries + west_africa_countries, asia_countries + east_africa_countries + west_africa_countries ]
"""
# Create custom sorting orders for Asia and Africa
asia_countries = ['SA', 'SEA', 'IEA']                               
africa_countries = ['WA', 'CA', 'EA']                                                                 

# Remove countries that are not present in the data
asia_countries = [country for country in asia_countries if country in dn_data.index]
africa_countries = [country for country in africa_countries if country in dn_data.index]

# Reorder rows and columns based on the custom sorting orders
dn = dn_data.loc[africa_countries + asia_countries, africa_countries + asia_countries ]
ds = ds_data.loc[africa_countries + asia_countries, africa_countries + asia_countries  ]
dnds = dnds_data.loc[africa_countries + asia_countries, africa_countries + asia_countries ]
"""                                    

sns.heatmap(dn, annot=False, cmap='YlGnBu', annot_kws={"size": 10}, fmt='.2f', cbar_kws={'label': 'dN'})
plt.title(sys.argv[2] + ' dN')
plt.xlabel('')
plt.ylabel('')
plt.tight_layout()
plt.show()

sns.heatmap(ds, annot=False, cmap='YlGnBu', annot_kws={"size": 10}, fmt='.2f', cbar_kws={'label': 'dS'})
plt.title(sys.argv[2] + ' dS')
plt.xlabel('')
plt.ylabel('')
plt.tight_layout()
plt.show()

sns.heatmap(dnds, annot=False, cmap='YlGnBu', annot_kws={"size": 10}, fmt='.2f', cbar_kws={'label': 'dN/dS'})
plt.title(sys.argv[2] + ' dN/dS')
plt.xlabel('')
plt.ylabel('')
plt.tight_layout()
plt.show()
"""


import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 3:
    print("Usage: ./heatmap.py Dxy_file Gene_name")
    sys.exit(1)

file = sys.argv[1]

# Load data from TSV file
data = pd.read_csv(file, sep='\t', header=None, names=['Country1', 'Country2', 'dN', 'dS', 'dN/dS'])
data.replace([np.inf, -np.inf], np.nan, inplace=True)
data.dropna(inplace=True)

# Drop duplicate entries with swapped order of countries
data['Country1'], data['Country2'] = np.sort(data[['Country1', 'Country2']], axis=1).T
data.drop_duplicates(subset=['Country1', 'Country2'], keep='first', inplace=True)

# Pivot the data
pivot_dn = data.pivot(index='Country1', columns='Country2', values='dN')
pivot_ds = data.pivot(index='Country1', columns='Country2', values='dS')
pivot_dnds = data.pivot(index='Country1', columns='Country2', values='dN/dS')

# Fill missing values
dn_data = pivot_dn.fillna(pivot_dn.T)
ds_data = pivot_ds.fillna(pivot_ds.T)
dnds_data = pivot_dnds.fillna(pivot_dnds.T)

# Define custom sorting orders
asia_countries = ['Cambodia', 'Thailand', 'Myanmar', 'Laos', 'Vietnam']
west_africa_countries = ['Gambia', 'Mali', 'Ghana', 'Mauritania', 'Senegal', 'Gabon', 'Cameroon', 'Guinea']
east_africa_countries =['Kenya', 'Tanzania']

# Filter out countries not present in the data
asia_countries = [country for country in asia_countries if country in dn_data.index]
east_africa_countries = [country for country in east_africa_countries if country in dn_data.index]
west_africa_countries = [country for country in west_africa_countries if country in dn_data.index]


# Reorder rows and columns based on custom sorting orders if possible, otherwise use data as it is
if all(country in dn_data.index for country in asia_countries + east_africa_countries + west_africa_countries):
    dn = dn_data.loc[asia_countries + east_africa_countries + west_africa_countries, asia_countries + east_africa_countries + west_africa_countries]
    ds = ds_data.loc[asia_countries + east_africa_countries + west_africa_countries, asia_countries + east_africa_countries + west_africa_countries]
    dnds = dnds_data.loc[asia_countries + east_africa_countries + west_africa_countries, asia_countries + east_africa_countries + west_africa_countries]
else:
    print("Custom sorting not possible. Using data as is.")
    dn = dn_data
    ds = ds_data
    dnds = dnds_data


print('dN values')
print(dn_data.describe())
print(dn_data)
print('dS values')
print(ds_data.describe())
print(ds_data)
print('dN/dS values')
print(dnds_data.describe())
print(dnds_data)



# Generate heatmaps
sns.heatmap(dn_data, annot=False, cmap='YlGnBu', annot_kws={"size": 10}, fmt='.2f', cbar_kws={'label': 'dN'})
plt.title(sys.argv[2] + ' dN')
plt.xlabel('')
plt.ylabel('')
plt.tight_layout()
plt.show()

sns.heatmap(ds_data, annot=False, cmap='YlGnBu', annot_kws={"size": 10}, fmt='.2f', cbar_kws={'label': 'dS'})
plt.title(sys.argv[2] + ' dS')
plt.xlabel('')
plt.ylabel('')
plt.tight_layout()
plt.show()

sns.heatmap(dnds_data, annot=False, cmap='YlGnBu', annot_kws={"size": 10}, fmt='.2f', cbar_kws={'label': 'dN/dS'})
plt.title(sys.argv[2] + ' dN/dS')
plt.xlabel('')
plt.ylabel('')
plt.tight_layout()
plt.show()
"""
