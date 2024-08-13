from Bio import SeqIO
import sys

def calculate_sequence_lengths(fasta_file):
    # Parse the FASTA file
    sequences = SeqIO.parse(fasta_file, "fasta")

    # Iterate through each sequence and print its length
    for record in sequences:
        print(f"Sequence: {record.id}, Length: {len(record.seq)}")


file=sys.argv[1]
calculate_sequence_lengths(file)

