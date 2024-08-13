import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_frequency(data_file, column, chart_type):
    # Load data
    df = pd.read_csv(data_file, sep='\t')
    # Remove non-finite values
    
    # Count the frequency of entries by column
    counts = df[column].value_counts().sort_values(ascending=False)

    # Plotting the data based on the specified chart type
    plt.figure(figsize=(10, 6))

    if chart_type == 'bar':
        plt.bar(counts.index, counts.values, color='skyblue')
        plt.xlabel(column)
        plt.ylabel('Frequency', fontsize=14)
        plt.title(column)
        plt.xticks(rotation='vertical', fontsize=14)
        plt.show()
    elif chart_type == 'pie':
        # Define a threshold for "small" categories
        threshold = 0.05 * counts.sum()  # 1% of total counts

        # Separate large and small categories
        large_counts = counts[counts >= threshold]
        small_counts = counts[counts < threshold]

        # Aggregate small categories into "Others"
        if not small_counts.empty:
            other_count = small_counts.sum()
            large_counts = pd.concat([large_counts, pd.Series([other_count], index=['Others'])])


        # Plot the pie chart
        plt.figure(figsize=(10, 7))
        plt.pie(large_counts, labels=large_counts.index, autopct='%1.1f%%', startangle=140, colors=plt.cm.Paired.colors)
        plt.axis('equal')  # Equal aspect ratio ensures that the pie is drawn as a circle.
        plt.xlabel(column, fontsize=14)
        plt.xticks(rotation='vertical', fontsize=14)
        plt.show()
    else:
        print("Invalid chart type. Supported types: bar, pie")

# For testing purposes, replace sys.argv with function parameters
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python plot_frequency.py <data_file> <column> <chart_type>")
        print("Supported chart bar, pie")

    data_file = sys.argv[1]
    column = sys.argv[2]
    chart_type = sys.argv[3]

    plot_frequency(data_file, column, chart_type)
