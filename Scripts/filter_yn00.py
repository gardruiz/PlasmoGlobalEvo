import os
import pandas as pd
import sys 
import subprocess


def filter_and_save(file_path):
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep='\s+')


    # Print the values of dN*N and dS*S before filtering
    print("Before filtering:")
    print("dN*N:", df['dN'] * df['N'])
    print("dS*S:", df['dS'] * df['S'])


    # Filter rows where either dN * N or dS * S is not equal to 0 using pandas
#    df = df[~((abs(df['dN'] * df['N']) == 0) & (abs(df['dS'] * df['S']) == 0))]


    # Filter rows where either dN * N or dS * S is not equal to 0
    df = df[(df['dN'] * df['N'] != 0) & (df['dS'] * df['S'] != 0)]

    # Print the values of dN*N and dS*S before filtering
    print("After filtering:")
    print("dN*N:", df['dN'] * df['N'])
    print("dS*S:", df['dS'] * df['S'])
    



    # Create 'no_zero' directory if it doesn't exist
    no_zero_dir = os.path.join(os.path.dirname(file_path), 'no_zero')
    if not os.path.exists(no_zero_dir):
        os.makedirs(no_zero_dir)

    # Get the file name without extension
    filename_no_ext = os.path.splitext(os.path.basename(file_path))[0]

    # Save the filtered DataFrame to a new file inside the 'no_zero' directory
    output_file_path = os.path.join(no_zero_dir, f"{filename_no_ext}.tab")
    df.to_csv(output_file_path, sep='\t', index=False)

# Path to the folder containing the files
folder_path = sys.argv[1]

# Iterate over each file in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.tab'):
        file_path = os.path.join(folder_path, file_name)
        # Apply filter and save to a new file
        filter_and_save(file_path)

