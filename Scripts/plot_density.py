#!/usr/bin/env python3

import numpy as np
import scipy
import pandas
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import allel
import sys
import os

def plot_windowed_variant_density(pos, window_size, output_folder="."):
    
    # setup windows 
    bins = np.arange(pos.min(), pos.max(), window_size)
    
    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2
    
    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size
    
    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')

    # Save the plot with the file name as title
    filename = os.path.splitext(os.path.basename(file))[0]
    plot_title = f"{filename}_variant_density"
    plot_filename = os.path.join(output_folder, f"{plot_title}.png")
    print(f"Saving plot to: {plot_filename}")
    plt.title(plot_title)
    print(plot_filename)
    print(output_folder)
    plt.savefig(plot_filename)

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script_name.py <vcf_file> <output_folder>")
    sys.exit(1)

file = sys.argv[1]
output_dir = sys.argv[2]

# Check if the output folder exists, create it if not
if not os.path.exists(output_dir):
    print(f"Creating output folder: {output_dir}")
    os.makedirs(output_dir)
else:
    print(f"Output folder already exists: {output_dir}")

vcf = allel.read_vcf(file)
pos = vcf["variants/POS"]

plot_windowed_variant_density(pos, window_size=1000, output_folder=output_dir)

