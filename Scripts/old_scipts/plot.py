import matplotlib.pyplot as plt

# Given quantities
values = [4700.65, 1156.71, 595.87, 426.647, 369.779, 242.722, 230.937, 212.323, 204.681, 190.303]
labels = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']

# Create bar plot
plt.bar(labels, values, color='skyblue')
plt.xlabel('Labels')
plt.ylabel('Values')
plt.title('Variance explained')
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
plt.tight_layout()  # Adjust layout for better visualization

# Show the plot
plt.show()

