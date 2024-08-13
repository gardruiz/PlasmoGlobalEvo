import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import allel
import sys

def plot_windowed_variant_density(pos, window_size, title=None):
    bins = np.arange(pos.min(), pos.max(), window_size)
    x = (bins[1:] + bins[:-1]) / 2
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size

    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)
    plt.axvline(221112, color='black', linestyle='--')
    plt.axvline(223145, color='black', linestyle='--')
    plt.show()

def process_vcf_chunk(file, region):
    vcf = allel.read_vcf(file, region=region)
    pos = vcf["variants/POS"]
    plot_windowed_variant_density(pos, window_size=1000, title=f"Chromosome 3 variant density ({region[0]} - {region[1]})")

def main(file):
    # Set the chunk size based on your needs
    chunk_size = 1000000000000

    # Read the total number of variants in the VCF file
    vcf = allel.read_vcf(file)
    total_variants = len(vcf['variants/POS'])
    regions = [(start, start + chunk_size) for start in range(0, total_variants, chunk_size)]

    # Process each chunk
    for region in regions:
        process_vcf_chunk(file, region)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py your_vcf_file.vcf")
        sys.exit(1)

    main(sys.argv[1])

