from Bio import SeqIO
import gffutils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

if len(sys.argv) != 4:
    print("Usage: python CDS_extractor.py <gene_sequence_file.fasta> <anotation_file.gtf> <output_file>")
    exit(1)


# Asing files to variables
gene_sequence_file = sys.argv[1]
gtf_file = sys.argv[2]

# Load GTF annotation
db = gffutils.create_db(gtf_file, dbfn=':memory:', force=True, keep_order=True, disable_infer_transcripts=True)

# Function to extract CDS sequence
def extract_cds_sequence(db, gene_id, gene_sequence):
    exons = db.features_of_type('exon', order_by='start', limit=db[gene_id].end - db[gene_id].start)
    cds_sequence = Seq('')
    for exon in exons:
        cds_sequence += gene_sequence[exon.start - 1:exon.end]
    return cds_sequence

# Replace 'output_file.fasta' with the desired output file path
output_file = sys.argv[3]

# Open the output file for writing
with open(output_file, 'w') as output_handle:
    # Iterate over sequences in the input gene sequence file
    for gene_record in SeqIO.parse(gene_sequence_file, 'fasta'):
        # Extract CDS sequence for each sequence
        gene_id = gene_record.id  # Assuming the sequence IDs are the gene IDs
        cds_sequence = extract_cds_sequence(db, gene_id, gene_record.seq)
        
        # Create a SeqRecord for the CDS sequence
        cds_record = SeqRecord(cds_sequence, id=gene_id, description="CDS sequence")

        # Write the SeqRecord to the output file
        SeqIO.write([cds_record], output_handle, 'fasta')

print(f"CDS sequences have been written to {output_file}")

