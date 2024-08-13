import sys
import pandas as pd
import matplotlib.pyplot as plt

def plot_average_fws_by_country(file1, file2):
    # Read the data from the files
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')

    # Merge the dataframes on the 'Ids' column
    merged_df = pd.merge(df1, df2, on='Sample')

    # Define a dictionary mapping countries to continents
    continent_mapping = {
        'Oceania': ['Indonesia','Papua New Guinea'],
        'South_America': ['Peru', 'Colombia', 'Venezuela'],
        'South_Asia': ['Bangladesh', 'India'],
        'West_Africa': ['Ghana', 'Mali', 'Gambia', 'Benin', 'Guinea', 'Senegal', 'Nigeria', 'Mauritania', "Côte d'Ivoire", 'Burkina Faso'],
        'Central_Africa': ['Cameroon','Gabon', 'Democratic Republic of the Congo','Republic of the Congo'],
        'Southeast_Asia': ['Vietnam', 'Cambodia', 'Myanmar', 'Thailand', 'Laos'],
        'East_Africa': ['Kenya', 'Tanzania', 'Malawi', 'Sudan', 'Mozambique', 'Ethiopia', 'Madagascar', 'Uganda']
    }

    # Add continent column to merged dataframe
    continent_column = []
    for country in merged_df['Country']:
        found_continent = False
        for continent, countries in continent_mapping.items():
            if country in countries:
                continent_column.append(continent)
                found_continent = True
                break  # No need to continue checking other continents if the country is found
        if not found_continent:
            print(f"No continent found for country: {country}")


    merged_df['Continent'] = continent_column

    # Calculate the average Fws value for each country
    average_fws_by_country = merged_df.groupby(['Continent', 'Country'])['Fws'].mean().sort_values()

    # Plot average Fws value for each country
    colors = plt.cm.tab10.colors  # Choose a colormap with distinct colors
    ax = average_fws_by_country.plot(kind='bar', color=colors)
    plt.xlabel('Country')
    plt.title('Average Fws Value for Each Country')
    ax.set_ylim(0.5, average_fws_by_country.max() + 0.1)  # Setting y-axis limits
    plt.ylabel('Fws')
    plt.xticks(rotation=90)
    plt.show()

# Example usage
file1 = 'Samples/Pf7_samples.txt'
file2 = sys.argv[1]
plot_average_fws_by_country(file1, file2)

