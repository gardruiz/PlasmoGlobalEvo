import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# Sample DataFrame (replace this with your actual data loading code)
df = pd.read_csv(sys.argv[1])

# Replace the above data with your actual DataFrame
# df = pd.read_csv("path_to_your_data.csv")

# Round snp_distance to 4 decimal places
df['snp_distance'] = df['snp_distance'].round(4)

# Get unique countries and count
unique_countries = df['country'].unique()
num_countries = len(unique_countries)

# Generate a custom color palette based on the number of unique countries
custom_palette = sns.color_palette("husl", num_countries)

# Get the tab10 palette colors
tab10_colors = sns.color_palette('tab10')

# Define colors for each group using tab10 colors
color_dict = {'Kenya': tab10_colors[0], 'Gabon': "#FFA500", 'Ghana': '#B2D63F', 'Sao_Tome': "#FF0000", 'Tanzania': tab10_colors[4], 'Cameroon': "green", 'Nigeria': "#00D000"}

# Create a dictionary to map each country to a specific color
countries = df['country'].unique()

# Set up the plot
plt.figure(figsize=(9, 4))
ax = plt.gca()
for country in countries:
    subset = df[df['country'] == country]
    sns.kdeplot(data=subset['snp_distance']*100, label=country, fill=True, alpha=0.3, linewidth=0, color=color_dict[country])


# Customize the plot
plt.xlabel('Average proportion of SNP difference(%)', fontsize=14, labelpad=10, fontweight='bold')
plt.ylabel('Density', fontsize=14, labelpad=10, fontweight='bold')
legend=plt.legend(loc='upper left', fontsize=14, frameon=False)
plt.gca().spines['bottom'].set_color('black')  # Adjust bottom axis line color
plt.gca().spines['bottom'].set_linewidth(2)    # Adjust bottom axis line width
plt.gca().spines['left'].set_color('black')    # Adjust left axis line color
plt.gca().spines['left'].set_linewidth(2)
plt.gca().tick_params(axis='both', which='major', width=2)
plt.xticks(fontweight='bold', fontsize=14)
plt.yticks(fontweight='bold', fontsize=14)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# Modify legend labels
for text in legend.get_texts():
    if text.get_text() == 'Sao_Tome':
        text.set_text('STP')
plt.tight_layout()
plt.show()

