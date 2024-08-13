import sys
import matplotlib.pyplot as plt

def count_zero_frequency(file_path):
    # Initialize a counter for zero frequency
    zero_count = 0
    # Open the file
    with open(file_path, 'r') as file:
        # Iterate through each line in the file
        for line in file:
            # Split the line by tab to get individual columns
            columns = line.strip().split('\t')
            # Check if the fourth column (index 3 since Python is 0-indexed) is '0'
            if columns[3] == '0':
                # Increment the zero count if the condition is met
                zero_count += 1
    return zero_count

def plot_zero_frequency(file_path):
    # Get the zero frequency count
    zero_frequency = count_zero_frequency(file_path)
    
    # Plotting
    plt.bar(['Zero Frequency'], [zero_frequency], color='blue')
    plt.xlabel('Frequency')
    plt.ylabel('Count')
    plt.title('Frequency of Zero Values in Fourth Column')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    plot_zero_frequency(file_path)
