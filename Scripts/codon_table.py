#this scrips converts a codon table to a python dictionary 

import sys 

if len(sys.argv) != 3:
    print ("Usage: ",sys.argv[0] + " codon_table output_file")
    sys.exit(1)

codon_table = sys.argv[1]
output_file = sys.argv[2]

# Read the data from the input file
with open(codon_table, 'r') as file:
    next(file)  # Skip the header line
    codon_data = [line.strip().split() for line in file if len(line.strip().split()) == 2]

# Create the codons dictionary and replace '*' with 'STOP'
codons = {}
for codon, amino_acid in codon_data:
    if amino_acid == "*":
        codons[codon] = "STOP"
    else:
        codons[codon]=amino_acid

# Save the codons dictionary as a .py script
with open(output_file, 'w') as script_file:
    script_file.write('codons = {\n')
    for codon, amino_acid in codons.items():
        script_file.write(f'    "{codon}": "{amino_acid}",\n')
    script_file.write('}\n')

print('Codons dictionary saved as', output_file, 'with "*" replaced by "STOP"')

