

import pandas as pd
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Usage: python plot_frequency.py <data_file>")
    exit(1)

# Read the input CSV file
df = pd.read_csv(sys.argv[1])

# Group by 'snp_distance' and 'country' and count the frequency
counts = df.groupby(['snp_distance', 'country']).size().reset_index(name='frequency')

# Filter for only 'Sao_Tome' and 'Nigeria'
filtered_counts = counts[counts['country'].isin(['Sao_Tome', 'Nigeria'])]

# Calculate the average SNP distance for each country
sao_tome_avg = df[df['country'] == 'Sao_Tome']['snp_distance'].mean()
nigeria_avg = df[df['country'] == 'Nigeria']['snp_distance'].mean()

# Plotting the data
plt.figure(figsize=(12, 8))

# Get unique countries (should be only 'Sao_Tome' and 'Nigeria' after filtering)
countries = filtered_counts['country'].unique()

# Plot each country's data separately
for country in countries:
    country_data = filtered_counts[filtered_counts['country'] == country]
    plt.bar(country_data['snp_distance'], country_data['frequency'], label=country, edgecolor='black')

# Add dotted lines for the average SNP distance
plt.axvline(x=sao_tome_avg, color='blue', linestyle='-', linewidth=1, label='Sao_Tome Average')
plt.axvline(x=nigeria_avg, color='orange', linestyle='-', linewidth=1, label='Nigeria Average')

# Specific samples and their line styles
specific_samples = {
    'MH-12-AC_HT10082010-1': 'dotted',
    'M11_HT04092010-01': 'dashdot'
}

# Add continuous lines for specific samples for both 'Sao_Tome' and 'Nigeria'
for sample, line_style in specific_samples.items():
    for country, color in zip(['Sao_Tome', 'Nigeria'], ['blue', 'orange']):
        sample_snp_distance = df[(df['sample'] == sample) & (df['country'] == country)]['snp_distance'].values
        if len(sample_snp_distance) > 0:
            plt.axvline(x=sample_snp_distance[0], color=color, linestyle=line_style, linewidth=1, label=f'{sample} ({country})')



plt.xlabel('SNP Distance')
plt.ylabel('Frequency')
plt.title('Frequency of SNP Distance by Country (Sao_Tome and Nigeria)')
plt.ylim(0, 20)  # Set y-axis limit to 20
plt.xticks(rotation='vertical')
plt.legend(title='Country')
plt.show()



"""
if len(sys.argv) != 2:
    print("Usage: python plot_frequency.py <data_file>")
    exit(1)

# Read the input CSV file
df = pd.read_csv(sys.argv[1])

# Group by 'snp_distance' and 'country' and count the frequency
counts = df.groupby(['snp_distance', 'country']).size().reset_index(name='frequency')

# Filter for only 'Sao_Tome' and 'Nigeria'
filtered_counts = counts[counts['country'].isin(['Sao_Tome', 'Nigeria'])]

# Calculate the average SNP distance for each country
sao_tome_avg = df[df['country'] == 'Sao_Tome']['snp_distance'].mean()
nigeria_avg = df[df['country'] == 'Nigeria']['snp_distance'].mean()

# Plotting the data
plt.figure(figsize=(12, 8))

# Get unique countries (should be only 'Sao_Tome' and 'Nigeria' after filtering)
countries = filtered_counts['country'].unique()

# Plot each country's data separately
for country in countries:
    country_data = filtered_counts[filtered_counts['country'] == country]
    plt.bar(country_data['snp_distance'], country_data['frequency'], label=country, edgecolor='black')

# Add dotted lines for the average SNP distance
plt.axvline(x=sao_tome_avg, color='blue', linestyle='--', linewidth=1, label='Sao_Tome Average')
plt.axvline(x=nigeria_avg, color='orange', linestyle='--', linewidth=1, label='Nigeria Average')


plt.xlabel('SNP Distance')
plt.ylabel('Frequency')
plt.title('Frequency of SNP Distance by Country (Sao_Tome and Nigeria)')
plt.xticks(rotation='vertical')
plt.legend(title='Country')
plt.show()
"""
