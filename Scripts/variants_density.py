import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import allel
import sys

def plot_windowed_variant_density(vcf_file, window_size, output_folder="."):
    """
    Plots windowed variant density.

    Parameters:
    - vcf_file: str, path to the VCF file
    - window_size: int, size of the genomic window
    - output_folder: str, path to the output folder (default is current directory)
    """
    
    # Read VCF file using scikit-allel
    vcf = allel.read_vcf(vcf_file, ['variants/CHROM', 'variants/POS', 'variants/DP', 'calldata/GT', 'ANN'])
    
    # Extract POS and INFO fields
    pos = vcf["variants/POS"]
    info_field = vcf["variants/ANN"]
    
    # Setup windows
    bins = np.arange(pos.min(), pos.max(), window_size)
    
    # Use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1]) / 2
    
    # Initialize arrays to store overall and missense variant densities
    y = np.zeros_like(x)
    y_missense = np.zeros_like(x)
    
    # Iterate over windows
    for i in range(len(bins) - 1):
        start, end = bins[i], bins[i + 1]
        
        # Filter missense variants within the current window
        missense_positions = [pos[j] for j in range(len(pos)) if start <= pos[j] < end and "missense_variant" in info_field[j]]
        
        # Compute histogram and density for overall variants within the window
        h, _ = np.histogram(pos, bins=[start, end])
        y[i] = h[0] / window_size
        
        # Compute histogram and density for missense variants within the window
        h_missense, _ = np.histogram(missense_positions, bins=[start, end])
        y_missense[i] = h_missense[0] / window_size
    
    # Plot overall variant density
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y, label='Overall Variants', color='blue')
    
    # Plot missense variant density
    ax.plot(x, y_missense, label='Missense Variants', color='red')

    # Plot SERA2 position
#    plt.axvline(x=317422 , color='black', linestyle='--', label='SERA2')
 #   plt.axvline(x=321552 , color='black', linestyle='--')

    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    ax.legend()

    # Save the plot with the file name as title
    plot_title = os.path.splitext(os.path.basename(vcf_file))[0] + "_variant_density"
    plot_filename = os.path.join(output_folder, f"{plot_title.replace(' ', '_')}.png")
    print(f"Saving plot to: {plot_filename}")
    plt.title(plot_title)
    plt.tight_layout()
    plt.savefig(plot_filename)

# Example usage:
# Assuming you have a VCF file 'your_file.vcf' and want to create a plot with a window size of 1000 bp

# plot_windowed_variant_density('your_file.vcf', window_size=1000, output_folder='your_output_folder')
plot_windowed_variant_density(sys.argv[1], window_size=100, output_folder='./')
