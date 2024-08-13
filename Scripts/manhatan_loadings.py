import matplotlib.pyplot as plt
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
            pcs = [float(pc_loading) for pc_loading in parts[5:]]  # Adjusted to capture loadings from the 6th to the 9th column
            snp_loadings[(chrom_numeric, position)] = pcs  # Store PC loadings for the current SNP

    return snp_loadings




def plot_manhattan(snp_loadings, pc_index, filename):
    # Sort SNP positions based on chromosome number and position
    sorted_snp_positions = sorted(snp_loadings.keys())

    # Extract sorted SNP positions and corresponding PC loadings
    sorted_positions = [pos for chrom, pos in sorted_snp_positions]
    pc_loadings = [snp_loadings[pos][pc_index] for pos in sorted_snp_positions]

    # Plot Manhattan plot
    plt.figure(figsize=(50, 4))  # Adjust the figsize to make the plot longer
    plt.scatter(sorted_positions, pc_loadings, color='blue', s=5)
    plt.xlabel('SNP position (ordered chromosomes)')
    plt.ylabel(f'Loadings on PC{pc_index + 1}')
    plt.title(f'PC{pc_index + 1}')
    
    # Remove box around plot
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


# Example usage
snp_loadings = parse_snp_loadings(sys.argv[1])
pc_no = int(sys.argv[2])
for i in range(pc_no):  # Plot PCs
    filename = f'/home/Daniel/Pictures/results/PCA/loadings/PC{i + 1}.png'
    plot_manhattan(snp_loadings, i, filename)
