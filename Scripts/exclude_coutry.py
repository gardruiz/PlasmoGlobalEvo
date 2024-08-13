"""
import sys
import os


def sampler(input_txt, country_list_file,  output):
    # Open the tab-delimited text file
    with open(input_txt, 'r') as f:
        # Split lines by tabs and create a list for each line
        lines = [line.strip().split('\t') for line in f]
    
    # Read the list of countries from the file
    with open(country_list_file, 'r') as countries_file:
        country_list = [line.strip() for line in countries_file]

    # Create a dictionary mapping sample IDs to country names for samples that passed QC
    data = {}
    # and the fifteenth column (index 14) contains QC information
    for line in lines[1:]:
        sample_id = line[0]
        country = line[2] 
        data[sample_id] = country

    if os.path.exists(output):
        os.remove(output) # Remove file if it exisits 
        print(f"file {output} was removed")
    
    # select samples from the specified country
    with open(output, "w") as file:
        for country_name in country_list:
            country_samples = {k: v for k, v in data.items() if v == country_name}
            for item in country_samples:
                # Separates the sample's name using a newline character
                file.write(item + "\t" + item  + "\n")


file = sys.argv[1]
country_list = sys.argv[2]
output = sys.argv[3]
sampler(file, country_list, output)
"""

import sys
import os

def sampler(input_txt, country_list_file, output):
    # Open the tab-delimited text file
    with open(input_txt, 'r') as f:
        # Split lines by tabs and create a list for each line
        lines = [line.strip().split('\t') for line in f]

    # Read the list of countries from the file
    with open(country_list_file, 'r') as countries_file:
        country_list = [line.strip() for line in countries_file]

    # Create a dictionary mapping sample IDs to country names for samples that passed QC
    data = {}
    # and the fifteenth column (index 14) contains QC information
    for line in lines[1:]:
        sample_id = line[0]
        country = line[2]
        data[sample_id] = country

    if os.path.exists(output):
        os.remove(output)  # Remove file if it exists
        print(f"File {output} was removed")

    # Open the output file once before the loop
    with open(output, "w") as file:
        # Iterate over each country
        for country_name in country_list:
            country_samples = {k: v for k, v in data.items() if v == country_name}
            print(country_name)
            print(country_samples)
            for item in country_samples:
                # Separates the sample's name using a newline character
                file.write(item + "\t" + item + "\n")

file = sys.argv[1]
country_list = sys.argv[2]
output = sys.argv[3]
sampler(file, country_list, output)

