import os
import sampler
import sys

if len(sys.argv) != 4:
        print("Usage: python multysampler.py <country_list_file> <num_samples> <samplesID_file>")
        sys.exit(1)

# Define the filename of the country list file
country_list_file = sys.argv[1]

# Define the number of samples for each country
num_samples_per_country = int(sys.argv[2])

# Define the samples ID file
samples_id_file = sys.argv[3]

with open(country_list_file, 'r') as file:
    for line in file:
        # Remove leading and trailing whitespace (e.g., newline characters)
        country_name = line.strip()
        
        # Call the sampler function for each country and save results in the Sa
        sampler.sampler(samples_id_file, country_name, num_samples_per_country)

