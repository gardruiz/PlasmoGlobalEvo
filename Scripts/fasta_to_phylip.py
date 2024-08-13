import sys
from Bio import SeqIO

def truncate_name(seq_name):
    return seq_name[:6]

def fasta_to_phylip(input_fasta):
    output_phylip = input_fasta.rsplit('.', 1)[0] + ".phy"
    sequences = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    
    # Truncate and fix duplicate sequence names
    seen_names = {}
    for seq_record in sequences.values():
        seq_name = seq_record.id
        truncated_name = truncate_name(seq_name)
        if truncated_name in seen_names:
            count = seen_names[truncated_name]
            seq_record.id = f"{truncated_name}_{count} "
            seen_names[truncated_name] += 1
        else:
            seq_record.id = f"{truncated_name} "
            seen_names[truncated_name] = 1

    print("Sequence names after truncation and handling duplicates:")
    for seq_record in sequences.values():
        print(seq_record.id)

    with open(output_phylip, "w") as output_handle:
        SeqIO.write(sequences.values(), output_handle, "phylip")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.fasta")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    fasta_to_phylip(input_fasta)




