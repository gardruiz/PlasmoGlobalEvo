import sys
import random
import os


def sampler(input_txt, country_name, num_samples, output_directory):
    # Open the tab-delimited text file
    with open(input_txt, 'r') as f:
        # Split lines by tabs and create a list for each line
        lines = [line.strip().split('\t') for line in f]

    # Check if the file has at least 15 columns (0-based index)
    if len(lines[0]) < 15:
        print("The input file does not have at least 15 columns.")
        sys.exit(1)

    # Create a dictionary mapping sample IDs to country names for samples that passed QC
    data = {}
    
    # Assuming the first column (index 0) contains sample IDs, the third column (index 2) contains country names,
    # and the fifteenth column (index 14) contains QC information
    for line in lines[1:]:
        sample_id = line[0]
        country = line[2] 
        qc_flag = line[13] 

        if qc_flag == "True":
            data[sample_id] = country

    
    # Create a directory to store the results
    os.makedirs(output_directory, exist_ok=True)  # Create the directory if it doesn't exist
    
    
    # Make a name tag for the output file using 'output_directory'
    output_name = os.path.join(output_directory, "{}.txt".format(country_name))

    # Make a name tag for the output file
   # output_name = "Samples/SamplesID/{}.txt".format(country_name)
    
    # Check if file exist
    if os.path.exists(output_name):
        os.remove(output_name) # Remove file if it exisits 
        print(f"file {output_name} was removed")
    # select samples from the specified country
    country_samples = {k: v for k, v in data.items() if v == country_name}

    # Check that the country is exisit in the data file
    if len(country_samples) == 0:
        print("No samples available for", country_name)
        sys.exit(0)

    if num_samples <= 0 or num_samples > len(country_samples):
        print(f"Invalid number of samples requested for {country_name}.")
        return


    # Sample the specified number of samples
    samples = random.sample(list(country_samples.keys()), int(num_samples))

    # Write the samples to a text file
    with open(output_name, "w") as file:
        for item in samples:
            # Separates the sample's name using a newline character
            file.write(item + "\n")

#if __name__ == "__main__":
#    if len(sys.argv) != 4:
 #       print("Usage: python script_name.py samples_file country_name num_samples")
 #       sys.exit(1)

#input_txt = sys.argv[1]
#country_name = sys.argv[2]
#num_samples = sys.argv[3]

#sampler(input_txt, country_name, num_samples)

