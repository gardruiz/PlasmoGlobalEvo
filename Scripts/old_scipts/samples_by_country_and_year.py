import pandas as pd
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Usage: python plot_frequency.py <data_file>")
    exit(1)

# Assuming your data is in a CSV file named 'your_data.csv'
df = pd.read_csv(sys.argv[1], sep='\t')

# Group the data by both Country and Year, sort by frequency and then sort countries by total frequency in descending order
grouped_data = df.groupby(['Country', 'Year']).size().unstack(fill_value=0)
sorted_countries = grouped_data.sum(axis=1).sort_values(ascending=False).index
grouped_data = grouped_data.loc[sorted_countries]

# Plotting the data as a grouped bar chart with legend to the right
plt.figure(figsize=(12, 8))
ax = grouped_data.plot(kind='bar', stacked=True, colormap='viridis')
plt.xlabel('Country')
plt.ylabel('Frequency')
plt.title('Sample Frequency by Country and Year')

# Move the legend to the right
ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.show()






