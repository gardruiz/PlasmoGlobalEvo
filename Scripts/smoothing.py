import os
import numpy as np
import pandas as pd
import sys

def calculate_smoothed_averages(data, global_avg_dn, global_avg_ds):
    """
    Calculate the smoothed average of dN and dS values using additive smoothing.
    """

    print("data lenght",len(data), data.columns)
    smoothed_dn = data["dN"] + global_avg_dn
    smoothed_ds = data["dS"] + global_avg_ds

    
    # Print sample data for verification
    print("Sample data for smoothing:", len(data), data.head())

    # Calculate the mean using pandas
    smoothed_dn_mean = smoothed_dn.mean() if not smoothed_dn.empty else 0
    smoothed_ds_mean = smoothed_ds.mean() if not smoothed_ds.empty else 0
    
    print("Smoothed dN and dS:", smoothed_dn_mean, smoothed_ds_mean)
    
    return smoothed_dn_mean, smoothed_ds_mean

def process_file(file_path, global_avg_dn, global_avg_ds):
    """
    Process a single file to calculate smoothed dN, dS, and dN/dS values.
    """
    data = pd.read_csv(file_path, delimiter='\s+', engine='python', skip_blank_lines=True)
    
    print("process file data", len(data))
    smoothed_dn, smoothed_ds = calculate_smoothed_averages(data, global_avg_dn, global_avg_ds)

    dn_ds_ratio = smoothed_dn / smoothed_ds if smoothed_ds != 0 else np.nan
    
    return smoothed_dn, smoothed_ds, dn_ds_ratio

def extract_countries(file_name):
    """
    Extract country1 and country2 from the file name.
    """
    base_name = os.path.basename(file_name)
    parts = base_name.split('.')[0].split('-')
    country1 = parts[0].split('_')[0]
    country2 = parts[1].split('_')[0]
    return country1, country2

def main(folder, output_file):
    global_sum_dn = 0
    global_sum_ds = 0
    global_count_dn = 0
    global_count_ds = 0
    111
    # First pass: Calculate global averages
    for file_name in os.listdir(folder):
        file_path = os.path.join(folder, file_name)
        df = pd.read_csv(file_path, delimiter='\s+', engine='python', skip_blank_lines=True)
         
        # Calculate sums and counts for dN and dS
        global_sum_dn += df['dN'].sum()
        global_sum_ds += df['dS'].sum()
        global_count_dn += df["dN"].astype(bool).sum()  # Count non-zero values
        global_count_ds += df["dS"].astype(bool).sum()  # Count non-zero values

    
    global_avg_dn = global_sum_dn / global_count_dn if global_count_dn != 0 else 0
    global_avg_ds = global_sum_ds / global_count_ds if global_count_ds != 0 else 0
    
    # Second pass: Process each file and output results
    with open(output_file, 'w') as out_file:
        for file_name in os.listdir(folder):
            file_path = os.path.join(folder, file_name)
            country1, country2 = extract_countries(file_name)
            print("Global average dN and dS:",global_avg_dn, global_avg_ds)
            smoothed_dn, smoothed_ds, dn_ds_ratio = process_file(file_path, global_avg_dn, global_avg_ds)
            out_file.write(f"{country1}\t{country2}\t{smoothed_dn}\t{smoothed_ds}\t{dn_ds_ratio}\n")

if __name__ == "__main__":
    folder = sys.argv[1]
    output_file = sys.argv[2]
    main(folder, output_file)

"""

import os
import pandas as pd
import sys

# Folder containing the files

def calculate_global_averages(folder_path):
    global_sum_dN = 0.0
    global_count_dN = 0
    global_sum_dS = 0.0
    global_count_dS = 0

    # Iterate over all files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                lines = file.readlines()
                
                # Skip the header (assuming it's the first line)
                for line in lines[1:]:
                    columns = line.split()

                    # Ensure that the line has the expected number of columns
                    if len(columns) >= 9:
                        try:
                            dN = float(columns[7])
                            dS = float(columns[9])

                            # Only consider non-zero dN and dS values
                            if dN != 0:
                                global_sum_dN += dN
                                global_count_dN += 1

                            if dS != 0:
                                global_sum_dS += dS
                                global_count_dS += 1

                        except ValueError:
                            print(f"Skipping line with invalid float conversion in file: {filename}")
                            continue

    # Calculate global averages
    global_avg_dN = global_sum_dN / global_count_dN if global_count_dN > 0 else 0
    global_avg_dS = global_sum_dS / global_count_dS if global_count_dS > 0 else 0

    return global_avg_dN, global_avg_dS

# Example usage:
folder = sys.argv[1]
global_avg_dN, global_avg_dS = calculate_global_averages(folder)
print(f"Global Average dN: {global_avg_dN}")
print(f"Global Average dS: {global_avg_dS}")

"""
