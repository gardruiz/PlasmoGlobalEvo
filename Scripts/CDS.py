#!usr/bin/env python3

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys
import os

def reverse_complement(sequence):
    complement = {'*': '*', 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'R': 'Y', 'Y': 'R', 'W': 'W', 'S': 'S', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'r': 'y', 'y': 'r', 'w': 'w', 's': 's', 'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b'}
    reversed_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement[base] for base in reversed_sequence)
    return reverse_complement_sequence


def read_fasta(file):
   sequences=[]
   with open(file, "r") as fasta_file:
       for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record)
   return sequences

def read_bed(file):
    cds_positions =[]
    gene_start=''
    sence=''
    with open(file, "r") as bed_file:
        for line in bed_file:
            fields = line.strip().split(' ')
            if fields[1]=='gene':
                gene_start=int(fields[2])
                sence += fields[4]
            else:
                start=int(fields[2])-gene_start
                end=int(fields[3])-gene_start
                cds_positions.append([start,end])
    return cds_positions, sence


def get_cds(sequence, positions, sence):
    cds = ""
    for start, end in positions:
        cds += sequence[start:end + 1]
    if sence == "-":
        cds = reverse_complement(cds)
    # Check if the last codon is a stop codon
    last_codon = cds[-3:]  
    stop_codons = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}  # Set of stop codons
    if last_codon in stop_codons:
        cds = cds[:-3]  # Remove the last codon
    return cds

if len(sys.argv) != 4:
    print("Usage: python CDS.py fasta_file bed_file output_file")
    exit(1)

file=sys.argv[1]
sequences=read_fasta(file)
cds_positions, sence = read_bed(sys.argv[2])
print("CDS positions: ", cds_positions)

# Initialize an empty list to store SeqRecord objecs
cds_records = []

# Generate the sequence ID

for sequence in sequences:
    header = sequence.description  # Get the header from the original FASTA record
    seq_id = header.split()[0]  # Use the first part of the header as the seq_id
    cds = get_cds(sequence.seq, cds_positions, sence)
    # Extract the file name without extension
    file_name_without_extension = os.path.basename(file).split('.')[0]
    
    
    cds_record   = SeqRecord(Seq(cds), id=seq_id)
    
    # Add the SeqRecord to the list
    cds_records.append(cds_record)
    
output_file=sys.argv[3]
SeqIO.write(cds_records, output_file, "fasta")
print(f"CDS sequences saved to {output_file}")
