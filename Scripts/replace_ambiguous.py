import sys

def replace_ambiguous_nucleotides(fasta_file):
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    # Iterate through the lines in the FASTA file
    for i in range(len(lines)):
        # Skip header lines
        if lines[i].startswith('>'):
            continue
        # Replace ambiguous nucleotides with 'N'
        else:
            lines[i] = lines[i].replace('r', 'n').replace('y', 'n').replace('w', 'n').replace('s', 'n').replace('m', 'n').replace('k', 'n').replace('b', 'n').replace('d', 'n').replace('h', 'n').replace('v', 'n').replace('R', 'N').replace('Y', 'N').replace('W', 'N').replace('S', 'N').replace('M', 'N').replace('K', 'N').replace('B', 'N').replace('D', 'N').replace('H', 'N').replace('V', 'N')

    # Write the modified lines back to the file
    with open(fasta_file, 'w') as f:
        f.writelines(lines)

if __name__ == "__main__":
    # Check if correct number of arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <fasta_file>")
        sys.exit(1)


fasta_file = sys.argv[1]

replace_ambiguous_nucleotides(fasta_file)

