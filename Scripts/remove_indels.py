import sys

def remove_insertions(consensus_file, temp_file):
    with open(consensus_file, 'r') as infile:
        lines = infile.readlines()

    processed_lines = []
    for line in lines:
        if line.startswith('>'):
            processed_lines.append(line)
        else:
            processed_line = ''.join([char for char in line if not char.islower()])
            processed_lines.append(processed_line)

    with open(temp_file, 'w') as outfile:
        outfile.writelines(processed_lines)


def add_deletions(temp_file, reference, header, output_file):
    ref_lines = reference.split('\n')
    
    with open(temp_file, 'r') as cons_file:
        cons_lines = cons_file.readlines()

    header
    ref_seq = ''
    cons_seq = ''
    
    # Extract headers and sequences
    for line in ref_lines:
        if not line.startswith('>'):
            ref_seq += line.strip()

    for line in cons_lines:
        if not line.startswith('>'):
            cons_seq += line.strip()
    
    print("Consensus sequence:\n", "length ", len(cons_seq), "\n", cons_seq)
    print("Reference sequence:\n",  "length ", len(ref_seq), "\n", ref_seq)
    # Process sequences
    final_seq = ''.join([ref_seq[i] if char == '*' else char for i, char in enumerate(cons_seq)])
    print("Final sequence:\n", "length ", len(final_seq), "\n", final_seq)
    # Format as FASTA string
    final_fasta = ">"+ header + '\n' + final_seq + '\n'

    with open(output_file, 'w') as outfile:
        outfile.write(final_fasta)


if __name__ == "__main__":
    consensus_file = sys.argv[1] # Input consensus sequence file
    reference = sys.argv[2]   # Reference sequence file
    
    output_file = sys.argv[3]      # Output final consensus sequence file
    temp_file = 'temp_no_insertions.fa'  # Temporary file to store intermediate result
    filename = consensus_file.split('/')[-1]
    filenae = filename.split('.')[0]
    header = filename.split('_')[0]
    remove_insertions(consensus_file, temp_file)
    add_deletions(temp_file, reference, header, output_file)
    # Clean up temporary file
    import os
  #  os.remove(temp_file)


