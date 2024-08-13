import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Usage: python plot_eigenvalues.py eigenvalues_file")
    exit(1)

# Read eigenvalues from the file
eigenvalues_file = sys.argv[1]
with open(eigenvalues_file, 'r') as file:
    eigenvalues = [float(line.strip()) for line in file]

# Plot the eigenvalues
plt.figure(figsize=(10, 6))
plt.bar(range(1, len(eigenvalues) + 1), eigenvalues)
plt.title('Variance explained')
plt.xlabel('Principal Component')
plt.ylabel('Eigenvalue')
plt.show()

