import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

# Sample code for reading the CSV file
# Uncomment and use sys.argv if running from command line
df1 = pd.read_csv(sys.argv[1])
df2 = pd.read_csv(sys.argv[2])
#df3 = pd.read_csv(sys.argv[3])
country1 = sys.argv[3]
country2 = sys.argv[4]
#country3 = sys.argv[6]
df1['country']=country1
df2['country']=country2
#df3['country']=country3
#df= pd.concat([df1, df2, df3], ignore_index=True)
df= pd.concat([df1, df2], ignore_index=True)

print (df)
print(df.info())
print(df.describe())

# Replace the above data with your actual DataFrame
# df = pd.read_csv("path_to_your_data.csv")
#df['snp_distance'] = df['snp_distance'].round(5)

# Group by 'snp_distance' and 'country' and count the frequency
counts = df.groupby(['snp_distance', 'country']).size().reset_index(name='frequency')

# Get unique countries and count
unique_countries = df['country'].unique()
num_countries = len(unique_countries)


# Pivot the DataFrame
pivot_df = counts.pivot(index='snp_distance', columns='country', values='frequency').fillna(0)
print(pivot_df)

# Generate a custom color palette based on the number of unique countries
#custom_palette = sns.color_palette("husl", num_countries)

# Define custom pastel blue and red colors in hexadecimal format
custom_palette = ['green',"#1E90FF", "#FF0000"] 
#["#ADD8E6", "#FFA07A"]  # Pastel blue and red colors

# Set the custom palette
sns.set_palette(custom_palette)
# Get data for only two specific countries (e.g., 'Sao_Tome' and 'Cameroon')
countries_to_plot = [country1, country2]
df = df[df['country'].isin(countries_to_plot)]


# Set up the plot
plt.figure(figsize=(12, 8))
ax = plt.gca()
# Plot the distribution of snp_distances for each country using KDE plots
for i, (country, data) in enumerate(df.groupby('country')):
    sns.kdeplot(data=data['snp_distance'], label=country, fill=True, alpha=0.3, linewidth=0, color=custom_palette[i])



"""
# Use the same colors for the vertical lines
plt.axvline(x=0.005531127629702, linestyle='--', alpha = 0.4, color=custom_palette[0], label='MH-02-AC_HAM15-000778 Nigeria')
plt.axvline(x=0.0057985745677654, linestyle='--', alpha = 0.4, color=custom_palette[1], label='MH-02-AC_HAM15-000778 Sao Tome')
plt.axvline(x=0.0054549108438964, linestyle='-', alpha = 0.4, color=custom_palette[0], label='MH-02-AC_HAM15-000572 Nigeria')                       
plt.axvline(x=0.0055816080475626, linestyle='-', alpha = 0.4, color=custom_palette[1], label='MH-02-AC_HAM15-000572 Sao Tome')
"""

#for country, data in df.groupby('country'):
 #   sns.kdeplot(data=data['snp_distance'], label=country, fill=True, alpha=0.4, linewidth=0)

# Customize the plot
plt.xlabel('SNP Distance', fontsize=14)
plt.ylabel('Density', fontsize=14)
plt.legend(loc='upper left', fontsize=14, frameon=False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

