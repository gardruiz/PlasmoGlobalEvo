import pandas as pd
import subprocess
import sys
import tempfile
import os

if len(sys.argv) != 4:
    print("Usage: python script.py locus_file vcf_file reference_fasta_file")
    exit(1)


df = pd.read_csv(sys.argv[1])
region_column = 'locus'

# If bcftools is not in the system's PATH replace 'bcftools_path' with the actual path to the bcftools executable
bcftools_path = 'bcftools'
vcf_file = sys.argv[2]
reference_fasta = sys.argv[3]

# Get unique samples from the VCF using bcftools query
bcftools_query_command = [bcftools_path, 'query', '-l', vcf_file]
samples = subprocess.check_output(bcftools_query_command, text=True).splitlines()

# Create an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['sample', 'locus', 'sequence'])



# Iterate over unique genomic regions
for genomic_region in df[region_column].unique():
    processed_region = genomic_region.replace('-', ':', 1).rsplit('-', 1)[0]

    # Use samtools faidx to get the subset reference
    subset_reference_command = [
        'samtools',
        'faidx',
        reference_fasta,
        processed_region
    ]
    subset_reference = subprocess.check_output(subset_reference_command, text=True)

    # Write the subset reference to a temporary file
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_fasta_file:
        temp_fasta_file.write(subset_reference)

    # Cut the VCF for the genomic region and save it to a temporary file
    fd, temp_vcf_name = tempfile.mkstemp(suffix='.vcf.gz')
    bcftools_view_command = [
        bcftools_path,
        'view',
        '-r', processed_region,
        '-O', 'z',  # Output in compressed VCF format
        '-o', temp_vcf_name,
        vcf_file
    ]
    subprocess.run(bcftools_view_command)
    
    # Index the temporary VCF file using tabix
    tabix_command = ['tabix', '-p', 'vcf', temp_vcf_name]
    subprocess.run(tabix_command)


    # Iterate through each sample
    for sample in samples:
        # Construct the bcftools consensus command
        bcftools_consensus_command = [
            bcftools_path,
            'consensus',
            '-s', sample,
            '-f', temp_fasta_file.name,
            temp_vcf_name,  
        ]

        # Run the bcftools consensus command and capture the output
        consensus_sequence = subprocess.check_output(bcftools_consensus_command, text=True)

        # Append the result to the DataFrame
        result_df = result_df.append({'sample': sample, 'locus': genomic_region, 'sequence': consensus_sequence.strip()}, ignore_index=True)

    # Remove the temporary files
    os.remove(temp_fasta_file.name)
    os.remove(temp_vcf_file.name)
    os.remove(temp_vcf_file.name + '.tbi')


# Save the result DataFrame to a tab-delimited file in the specified output directory
result_df.to_csv('output_table.tsv', sep='\t', index=False)

print("Output table generated.")
