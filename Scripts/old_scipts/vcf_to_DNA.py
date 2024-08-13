""" This scrip converits VCF to DNA, it will convert all positions in the VCF file if no custom range is give"""

import sys
from Bio import SeqIO
import vcf

def extract_sequences_from_vcf(vcf_file, reference_genome_file):
    sequences = {}
    # Load the reference genome
    reference_genome = SeqIO.to_dict(SeqIO.parse(reference_genome_file, "fasta"))

    # Parse the VCF file to get the genomic range
    genomic_range, genomic_range_start, genomic_range_end, chrom = extract_genomic_range_from_vcf_header(vcf_file)
    print(chrom)  
    # Define the chromosome sequence
    chromosome_sequence = reference_genome[chrom].seq
    

    # Initialize sequences for each sample
    vcf_reader = vcf.Reader(filename=vcf_file)
        # Modify the reference sequence acording to the vcf file
    for record in vcf_reader:
        for sample in record.samples:
            sample_name = str(sample)
            sequences[sample_name]= chromosome_sequence[genomic_range_start:genomic_range_end] 
            if str(sample.gt_alt[0]) != '*':
                chrom_position = record.POS -1
                gene_position = chrom_position - genomic_range_start
                # If ALT is not an asterisk, include the alternate allele
                sequences[sample_name] = (
                sequences[sample_name][:gene_position]
                + str(sample.gt.alt[0])
                + sequences[sample_name][gene_position + 1 :]
                )

    sequence_list = [(sequence_name, sequence) for sequence_name, sequence in sequences.items()]
    return sequence_list

def extract_genomic_range_from_vcf_header(vcf_file):
    genomic_range_start, genomic_range_end = None, None
    with open(vcf_file, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('##bcftools_viewCommand=view -r'):
                parts = line.split('-r ')[1].split(' ')
                genomic_range = parts[0]
                # Extract the start and end positions from the genomic range
                positions = genomic_range.split(':')[1].split('-')
                genomic_range_start = int(positions[0]) -1
                genomic_range_end = int(positions[1]) -1
                chromosome = genomic_range.split(":")
                genomic_range= f"{chromosome[0]}:{genomic_range_start}-{genomic_range_end}"
                chrom = chromosome[0]
                break

    return genomic_range, genomic_range_start, genomic_range_end, chrom

# Example usage
vcf_file = sys.argv[1]
reference_genome_file = sys.argv[2]
result = extract_sequences_from_vcf(vcf_file, reference_genome_file)
output_fasta_file="test.fasta"
# Write the sequences to the FASTA file
with open(output_fasta_file, "w") as output_handle:
    for name, sequence in result:
        seq_record = SeqIO.SeqRecord(sequence, id=name)
        SeqIO.write(seq_record, output_handle, "fasta")
