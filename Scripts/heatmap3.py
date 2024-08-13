import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Define custom color palette
custom_palette = sns.color_palette("light:b", as_cmap=True)

# Function to read data from a file and extract gene name
def read_data(file_path):
    gene_name = os.path.splitext(os.path.basename(file_path))[0]
    data = pd.read_csv(file_path, sep='\t', header=None, names=['Country1', 'Country2', 'dN', 'dS', 'dN/dS'])
    
    
    # Combine 'Country1' and 'Country2' into a single column
    data['Continent'] = data['Country1'] + '-' + data['Country2']
    
    # Drop 'Country1' and 'Country2' columns
    data.drop(['Country1', 'Country2'], axis=1, inplace=True)
    

    return gene_name, data

# Function to create heatmaps
def create_heatmap(df, title):

    # Replace NaN values with zero
    df_filled = df.fillna(0)
    # Transpose the DataFrame
    df_transposed = df_filled.T
    
    plt.figure(figsize=(10, 6))
    sns.heatmap(df_transposed, cmap=custom_palette, fmt=".2f", linewidths=.5)
    plt.title(title)
    plt.xlabel('Continent')
    plt.ylabel('Genes')
    plt.tight_layout()
    plt.show()



# Main function to iterate through files and create heatmaps
def main(folder_path):
    dn_data = pd.DataFrame()
    ds_data = pd.DataFrame()
    dnds_data = pd.DataFrame()

    for file_name in os.listdir(folder_path):
        if file_name.endswith('.tab'):  # Assuming all files are tab-delimited text files
            file_path = os.path.join(folder_path, file_name)
            gene_name, df = read_data(file_path)
            print(df)

            # Set gene name as index
            df.set_index('Continent', inplace=True)

            dn_data[gene_name] = df['dN']
            ds_data[gene_name] = df['dS']
            dnds_data[gene_name] = df['dN/dS']



    # Create heatmaps
    print(dn_data)
    print(dnds_data)
    create_heatmap(dn_data, 'dN Values')
    create_heatmap(ds_data, 'dS Values')
    create_heatmap(dnds_data, 'dN/dS Values')






folder_path = sys.argv[1]
main(folder_path)





