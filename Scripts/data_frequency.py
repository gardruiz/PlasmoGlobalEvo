import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

if len(sys.argv) != 3:
    print("Usage: python path_to/data_frequency.py data_file column_name")
    exit(1)

file = sys.argv[1]
column = sys.argv[2]

# Load your CSV file into a Pandas DataFrame
df = pd.read_csv(file, sep='\t')

# Group by 'column' and count the frequency for each combination of the specified column and 'year'
grouped_data = df.groupby([column, 'year']).size().unstack(fill_value=0)

# Calculate the total frequency for each country
grouped_data['total'] = grouped_data.sum(axis=1)

# Sort the DataFrame by the total frequency in descending order
grouped_data = grouped_data.sort_values(by='total', ascending=False)

# Drop the 'total' column before plotting
grouped_data.drop(columns='total', inplace=True)

# Use Seaborn for a more visually appealing plot
palette = sns.color_palette("crest", as_cmap=True)  # You can choose a different color palette
ax = grouped_data.plot(kind='bar', stacked=True, width=0.8, cmap=palette)
ax.set_xlabel(column, fontsize=14)
ax.set_ylabel("Frequency", fontsize=14)
ax.set_title(f"Frequency of samples by {column}", fontsize=16)
ax.tick_params(axis='x', labelrotation=45)  # Rotate x-axis labels for better visibility
ax.legend(title='Year', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)  # Add a legend

# Remove top and right spines for a cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# Show the plot (you can also save it using plt.savefig() if needed)
plt.tight_layout()
plt.show()


