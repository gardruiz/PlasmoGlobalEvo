import csv
import sys
import pandas as pd
# Check the number of arguments
if len(sys.argv) != 4:
    print("Usage: python scritp.py file1.csv file2.csv output_file")
    exit(1)

# Asing the CVS files to variables
file1 = sys.argv[1]
file2 = sys.argv[2]
output_file = sys.argv[3]

# Convert files into pandas dataframes
df1=pd.read_csv(file1)
df2=pd.read_csv(file2)

    
merged_df = pd.merge(df1, df2, left_on="sampleID", right_on= "SeqName",  how="inner")

merged_df.to_csv(output_file, index=False)

