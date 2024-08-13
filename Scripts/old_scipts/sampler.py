import csv
import random
import sys
import os

"""THIS SCRIPT TAKES A SAMPLES.CSV FILE FROM MalariaGENE, A COUNTRY NAME AND A NUMBER OF SAMPLES AS AN INPUT AND RETURS THE SPECIFED
NUMBER OF RANDOM SAMPLES AS TEXT FILE"""


# Open the file
with open(sys.argv[1], 'r') as f:
    # Create a CSV reader
    reader = csv.reader(f,)
    # Skip the header row
    next(reader)

    # Create an empty dictionary
    data = {}

    # Iterate over the rows
    for row in reader:
        # Get the key and value from the row
        key = row[0]
        value =row[5]

        # Add the key and value to the dictionary
        data[key] = value

#Makes a name tag for the output file
output_name = "{}.txt".format(sys.argv[2])

# Check if the output file already exists and remove it
if os.path.exists(output_name):
    os.remove(output_name)
    print(f"Removed existing {output_name}")

#Selects the samples from the specified country
sample = {k: v for k, v in data.items() if v == sys.argv[2]}

#Samples the specified number of samples
samples = random.sample(list(sample.keys()), int(sys.argv[3]))

#Writes the samples to a text file
with open(output_name, "w") as file:
    for item in samples:
        #Separates the samples 
        file.write(item + "\n")


