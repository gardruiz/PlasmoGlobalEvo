import sys

def calculate_variance_explained(eigenval_file):
    # Read eigenvalues from the .eigenval file
    with open(eigenval_file, 'r') as f:
        eigenvalues = [float(line.strip()) for line in f]
        print(eigenvalues)
    # Calculate total variance
    total_variance = sum(eigenvalues)
    print(total_variance)
    # Calculate variance explained by each PC
    variance_explained = [(eigval / total_variance) * 100 for eigval in eigenvalues]

    return variance_explained

# Example usage
eigenval_file = sys.argv[1]
variance_explained = calculate_variance_explained(eigenval_file)

# Print the percentage of variance explained by each PC
for i, pc_variance in enumerate(variance_explained, start=1):
    print(f"PC{i}: {pc_variance:.2f}%")

