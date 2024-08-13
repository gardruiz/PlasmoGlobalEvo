import sys
import subprocess
import matplotlib.pyplot as plt

def count_variants(vcf_dir):
    # Run shell command to count variants in each VCF file and save to a text file
    subprocess.run(['bash', '-c', f'for vcf_file in "{vcf_dir}"/*.vcf.gz; do filename=$(basename -- "$vcf_file"); filename_no_ext="${{filename%.vcf.gz}}"; variant_count=$(bcftools view -H "$vcf_file" | wc -l); echo "$filename_no_ext $variant_count" >> variant_counts.txt; done'])

    # Read variant counts from the text file
    with open('variant_counts.txt', 'r') as file:
        lines = file.readlines()

    # Extract filenames and variant counts
    filenames = [line.split()[0] for line in lines]
    variant_counts = [int(line.split()[1]) for line in lines]

    # Plot the number of variants per file
    plt.figure(figsize=(10, 6))
    plt.bar(filenames, variant_counts)
    plt.xlabel('VCF File')
    plt.ylabel('Number of Variants')
    plt.title('Number of Variants per VCF File')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <vcf_directory>")
        sys.exit(1)

    vcf_directory = sys.argv[1]
    count_variants(vcf_directory)

