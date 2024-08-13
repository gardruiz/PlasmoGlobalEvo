import numpy as np
import scipy
import pandas
import matplotlib as mpl
import matplotlib.pyplot as plt
#%matplotlib inline
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import allel; print('scikit-allel', allel.__version__)
import sys
import h5py

def plot_variant_hist(f, bins=30):
    x = variants[f]
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.despine(ax=ax, offset=10)
    ax.hist(x, bins=bins)
    ax.set_xlabel(f)
    ax.set_ylabel('No. variants')
    ax.set_title('Variant %s distribution' % f)
    plt.show()

file=sys.argv[1]
allel.vcf_to_hdf5(file, 'file.h5', fields='*', overwrite=True)
variants=h5py.File('file.h5', mode='r')
print(variants)
