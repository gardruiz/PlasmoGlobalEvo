
from Bio import AlignIO

# Load the alignment file
alignment = AlignIO.read("alignment_file.fasta", "fasta")

# Define a function to identify indels in the alignment
def identify_indels(alignment):
    indels = []
    for record in alignment:
        seq = str(record.seq)
        ref_seq = str(alignment[0].seq)  # Reference sequence
        for i, (base_ref, base_seq) in enumerate(zip(ref_seq, seq)):
            if base_ref != base_seq:
                # Check if it's a gap in the reference sequence
                if base_ref == '-':
                    indels.append((i, "insertion", base_seq))
                # Check if it's a gap in the query sequence
                elif base_seq == '-':
                    indels.append((i, "deletion", base_ref))
    return indels

# Define function to adjust exon coordinates
def adjust_exon_coordinates(exons, indels):
    for exon in exons:
        for position, indel_type, base in indels:
            # Check if exon overlaps with indel
            if exon.start <= position <= exon.end:
                if indel_type == "insertion":
                    exon.start += 1
                    exon.end += 1
                elif indel_type == "deletion":
                    if exon.start > position:
                        exon.start -= 1
                    if exon.end > position:
                        exon.end -= 1
    return exons

# Load and process GTF file
exons = []  # Store exon coordinates from GTF file
# Read exon coordinates from GTF file and store in 'exons'

# Identify indels in the alignment
indels = identify_indels(alignment)

# Adjust exon coordinates based on identified indels
adjusted_exons = adjust_exon_coordinates(exons, indels)

# Update GTF file with adjusted exon coordinates
# Write adjusted_exons to the GTF file

