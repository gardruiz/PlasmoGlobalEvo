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

# Convert snp_distance to string to handle small values in the plot
counts['snp_distance'] = counts['snp_distance'].astype(float)

# Plotting the data
plt.figure(figsize=(12, 8))

# Get unique countries
countries = counts['country'].unique()

# Plot each country's data separately
for country in countries:
    country_data = counts[counts['country'] == country]
    plt.bar(country_data['snp_distance'], country_data['frequency'], label=country, edgecolor='black')

plt.xlabel('SNP Distance')
plt.ylabel('Frequency')
plt.title('Frequency of SNP Distance by Country')
plt.xticks(rotation='vertical')
# Set x-axis range

plt.legend(title='Country')
plt.show()

