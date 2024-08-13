import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Load the PCA output data
pca_df = pd.read_csv(sys.argv[1], delim_whitespace=True)

# Extract the principal components for clustering
pc_columns = [f'PC{i}' for i in range(1, 11)]
pc_data = pca_df[pc_columns]

# Perform k-means clustering
kmeans = KMeans(n_clusters=4, random_state=40)  # Change n_clusters to the desired number of clusters
pca_df['Cluster'] = kmeans.fit_predict(pc_data)

# Plot the PCA results with clusters
plt.figure(figsize=(10, 8))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', s=100, alpha=0.8)

# Add cluster centers to the plot
centers = kmeans.cluster_centers_
plt.scatter(centers[:, 0], centers[:, 1], c='black', s=200, alpha=0.5, marker='X')

# Add plot title and labels
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend(title='Cluster')
plt.show()

# Save the plot
plt.savefig('pca_with_clusters.png')

