import pandas as pd
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Usage: python plink_missing.py input_file")
    exit(1)

file=sys.argv[1]
# Load sample missing data report
sample_missing_data = pd.read_csv(file, sep='\s+')

# Plot histogram of sample missingness rates
plt.hist(sample_missing_data["F_MISS"], bins=20, edgecolor='black')
plt.title("Histogram of Sample Missingness Rates")
plt.xlabel("Missingness Rate")
plt.ylabel("Frequency")
plt.show()

# Load variant missing data report
variant_missing_data = pd.read_csv("plink2.vmiss", sep='\s+')

# Plot histogram of variant missingness rates
plt.hist(variant_missing_data["F_MISS"], bins=20, edgecolor='black')
plt.title("Histogram of Variant Missingness Rates")
plt.xlabel("Missingness Rate")
plt.ylabel("Frequency")
plt.show()

