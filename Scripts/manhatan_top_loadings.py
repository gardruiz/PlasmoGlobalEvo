import matplotlib.pyplot as plt
import numpy as np
import sys

def parse_snp_loadings(file_path):
    snp_loadings = {}  # Dictionary to store SNP positions and their corresponding PC loadings

    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            parts = line.strip().split('\t')
            snp_id = parts[1]  # SNP ID in the second column

            # Extract chromosome number and position from SNP ID
            chrom_numeric = int(snp_id.split('_')[1])  # Extract the part between the first and second underscores
            position = int(snp_id.split(':')[1])  # Extract position from SNP ID

            # Capture PC loadings separately for each PC
            pcs = [float(pc_loading) for pc_loading in parts[5:9]]  # Adjusted to capture loadings from the 6th to the 9th column
            snp_loadings[(chrom_numeric, position)] = pcs  # Store PC loadings for the current SNP

    return snp_loadings


def extract_pc_loadings(snp_loadings, pc_index):
    """Extract only the loadings for the specified PC."""
    pc_specific_loadings = {}
    for pos, loadings in snp_loadings.items():
        pc_specific_loadings[pos] = loadings[pc_index]
    return pc_specific_loadings


def identify_top_negative_variants(pc_loadings, top_percentage):
    """Identify the top percentage of most negative loadings."""
    negative_loadings = [loading for loading in pc_loadings.values() if loading < 0]

    # Debug: Plot the distribution of negative loadings
    plt.hist(negative_loadings, bins=30, color='blue', alpha=0.7)
    plt.xlabel('PC Loadings')
    plt.ylabel('Frequency')
    plt.title('Distribution of Negative Loadings')
    
    # Calculate the threshold for the top percentage of most negative loadings
    if negative_loadings:  # Ensure there are negative loadings to process
        threshold = np.percentile(negative_loadings, 100 *  top_percentage)
        plt.axvline(x=threshold, color='red', linestyle='--', label=f'Threshold: {threshold:.2f} ({top_percentage:.2f}% most negative)')
    
    plt.legend()
    plt.show()

    # Identify variants with loadings less than or equal to the threshold
    top_negative_variants = {pos: loading for pos, loading in pc_loadings.items() if loading <= threshold}
    
    return top_negative_variants


def save_variants_to_file(variants, file_path):
    """Save the top negative variants to a file."""
    with open(file_path, 'w') as file:
        file.write('Chromosome\tPosition\tPC_Loading\n')
        for (chrom, pos), loading in variants.items():
            file.write(f'{chrom}\t{pos}\t{loading}\n')


def plot_manhattan(snp_loadings, pc_loadings, top_negative_variants):
    """Plot a Manhattan plot highlighting the top negative variants."""
    # Sort SNP positions based on chromosome number and position
    sorted_snp_positions = sorted(snp_loadings.keys())

    # Extract sorted SNP positions and corresponding PC loadings
    sorted_positions = [pos for chrom, pos in sorted_snp_positions]
    pc_loadings_sorted = [pc_loadings[pos] for pos in sorted_snp_positions]

    # Identify the positions of the top negative variants for highlighting
    top_negative_positions = set(top_negative_variants.keys())
    colors = ['red' if pos in top_negative_positions else 'blue' for pos in sorted_snp_positions]

    # Plot Manhattan plot
    plt.figure(figsize=(100, 3))  # Adjust the figsize to make the plot longer
    plt.scatter(sorted_positions, pc_loadings_sorted, c=colors, s=5)
    plt.xlabel('SNP position (ordered chromosomes)')
    plt.ylabel(f'Loadings on PC{pc_index + 1}')
    plt.title(f'PC{pc_index + 1}')
    
    # Remove box around plot
    plt.xticks([])
    plt.tight_layout()
    plt.show()


# Example usage
if __name__ == "__main__":
    snp_loadings = parse_snp_loadings(sys.argv[1])
    pc_index = 1  # Index for PC2 (0-based index, so 1 means PC2)

    # Extract only the loadings for the specified PC
    pc_loadings = extract_pc_loadings(snp_loadings, pc_index)
    
    top_negative_variants = identify_top_negative_variants(pc_loadings, top_percentage=0.05)
    save_variants_to_file(top_negative_variants, 'top_negative_variants_PC2.txt')

    plot_manhattan(snp_loadings, pc_loadings, top_negative_variants)
