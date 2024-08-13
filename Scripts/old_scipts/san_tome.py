import pandas as pd
import subprocess
import sys
import tempfile
import os

if len(sys.argv) != 4:
    print("Usage: python script.py locus_file vcf_file samplesID_file")
    exit(1)

df = pd.read_csv(sys.argv[1])
region_column = 'locus'
# If bcftools is not in the system's PATH, replace 'bcftools_path' with the actual path to the bcftools executable
bcftools_path = 'bcftools'
vcf_file = sys.argv[2]
reference_fasta = 'Genome/Pfalciparum.genome.fasta'


# Get unique samples from the VCF using bcftools query
#bcftools_query_command = [bcftools_path, 'query', '-l', vcf_file]
#samples = subprocess.check_output(bcftools_query_command, text=True).splitlines()

# Get samples from the sample id file
sample_file = open(sys.argv[3], 'r')
data = sample_file.read()
samples = data.split('\n')
print(samples)


# Create an empty list to collect the results
result_rows = []

# Iterate over unique genomic regions
for genomic_region in df[region_column].unique():
    # Preprocess the region
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

    # Normalize the temporary VCF file
    '''
    bcftools_norm_command = [
        bcftools_path,
        'norm',
        '-a',
        '-f', reference_fasta,
        '-Oz',
        '-o', temp_vcf_name,
        temp_vcf_name
    ]
    subprocess.run(bcftools_norm_command)
'''
    # Index the normalized VCF file using tabix
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

        # Append the result to the list
        result_rows.append({'sample': sample, 'locus': processed_region, 'sequence': consensus_sequence.strip()})

    # Remove the temporary files
    os.remove(temp_fasta_file.name)
    os.remove(temp_vcf_name)
    os.remove(temp_vcf_name + '.tbi')

# Convert the list of rows to a DataFrame
result_df = pd.DataFrame(result_rows, columns=['sample', 'locus', 'sequence'])

# Save the result DataFrame to a tab-delimited file in the specified output directory
result_df.to_csv('output_table.tsv', sep='\t', index=False)

print("Output table generated.")

