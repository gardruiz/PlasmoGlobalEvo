import pandas as pd
import sys
import matplotlib.pyplot as plt
import os



def concatenate_and_check_files(folder_path):
    # Get list of files in the folder
    files = os.listdir(folder_path)
    # Remove any non-.tab files
    files = [file for file in files if file.endswith('.tab')]

    # Initialize an empty DataFrame to store the concatenated data
    concatenated_df = pd.DataFrame()

    # Flag to track if unusual values are found
    unusual_values_found = False

    # Iterate over each file in the folder
    for file in files:
        file_path = os.path.join(folder_path, file)
        # Read the current file into a DataFrame
        df = pd.read_csv(file_path, sep='\s+')
        # Calculate (dN * N) and (dS * S)
        df['dN*N'] = df['dN'] * df['N']
        df['dS*S'] = df['dS'] * df['S']
        # Check if any value of dN*N or dS*S is greater than 100
        if (df['dN*N'] > 100).any() or (df['dS*S'] > 100).any():
            print(f"Unusual dN*N or dS*S values found in file: {file}")
            unusual_values_found = True
            # Concatenate the current DataFrame to the concatenated DataFrame
            concatenated_df = pd.concat([concatenated_df, df])
        else:
            # Concatenate the current DataFrame to the concatenated DataFrame
            concatenated_df = pd.concat([concatenated_df, df])

    # If no unusual values were found, return the concatenated DataFrame
    return concatenated_df

# Usage
folder = sys.argv[1]
df = concatenate_and_check_files(folder)
print("Concatenated DataFrame:")
print(df)



print("dN", df['dN'])

print("N", df['N'])


# Calculate (dN * N) and (dS * S)
df['dN*N'] = df['dN'] * df['N']
df['dS*S'] = df['dS'] * df['S']
print("dN*N", df['dN*N'])
print("dS*S", df['dS*S'])

print(df.describe())

# Plotting
plt.figure(figsize=(10, 6))

plt.subplot(3, 3, 1)
plt.hist(df['dN'], bins=20, color='grey', alpha=0.7)
plt.title('Frequency of dN')
plt.xlabel('dN')
plt.ylabel('Frequency')
plt.xticks(rotation=45)

plt.subplot(3, 3, 2)
plt.hist(df['dS'], bins=20, color='grey', alpha=0.7)
plt.title('Frequency of dS')
plt.xlabel('dS')
plt.ylabel('Frequency')

plt.subplot(3, 3, 3)
plt.hist(df['omega'], bins=20, color='grey', alpha=0.7)
plt.title('Frequency of omega')
plt.xlabel('omega')
plt.ylabel('Frequency')

plt.subplot(3, 3, 4)
plt.hist(df['dN*N'], bins=20, color='grey', alpha=0.7)
plt.title('Frequency of dN*N')
plt.xlabel('dN*N')
plt.ylabel('Frequency')

plt.subplot(3, 3, 5)
plt.hist(df['dS*S'], bins=20, color='grey', alpha=0.7)
plt.title('Frequency of dS*S')
plt.xlabel('dS*S')
plt.ylabel('Frequency')

plt.subplot(3, 3, 6)
plt.hist(df['S'], bins=20, color='grey', alpha=0.7)
plt.title('Frequency of S')
plt.xlabel('S')
plt.ylabel('Frequency')

plt.subplot(3, 3, 7)
plt.hist(df['N'], bins=20, color='grey', alpha=0.7)
plt.title('Frequency of N')
plt.xlabel('N')
plt.ylabel('Frequency')


plt.tight_layout()
plt.show()

