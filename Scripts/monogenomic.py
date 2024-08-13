import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('metadata/Pf7_fws.txt', sep='\t')

# Calculate bin edges
bin_edges = np.histogram_bin_edges(df["Fws"], bins=100)

# Identify the values in the last bin
last_bin_min = bin_edges[-2]
last_bin_max = bin_edges[-1]
df_fws = df[(df["Fws"] >= last_bin_min) & (df["Fws"] <= last_bin_max)]


# Load the dataset with sample and country information
df_country = pd.read_csv('Samples/Pf7_samples.txt', sep='\t')
# Merge the datasets on the "Sample" column
merged_df = pd.merge(df_fws, df_country, on="Sample")
merged_df.drop(columns=['Fws'], inplace=True)
print(merged_df.columns)
merged_df['year'] = merged_df['year'].astype(int)
merged_df.to_csv("monogenomic_infections.csv", index = False, sep='\t')
# Group by "Country" and count the frequency of entries
country_counts = merged_df['Country'].value_counts().reset_index()
country_counts.columns = ['Country', 'Frequency']

# Plot the frequency using a bar plot
plt.figure(figsize=(10, 6))
sns.barplot(data=country_counts, x='Country', y='Frequency')
plt.title('Frequency of Entries by Country')
plt.xlabel('Country')
plt.ylabel('Frequency')
plt.xticks(rotation=90)  # Rotate the x labels for better readability
plt.show()

