import pandas as pd
import argparse

def filter_samples_by_year(input_file, output_file, year_range):
    # Read the tab-separated file into a DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Extract start and end years from the range
    start_year, end_year = map(int, year_range.split('-'))

    # Generate a list of years within the specified range
    years = list(range(start_year, end_year + 1))

    # Filter samples based on the specified years
    filtered_samples = df[df['year'].isin(years)]['Sample']

    result_df = pd.DataFrame({'Sample': filtered_samples, 'Sample_rep': filtered_samples})

    # Write the filtered samples to the output file
    result_df.to_csv(output_file, sep='\t', index=False)

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Filter samples by year from a tab-separated file')
    parser.add_argument('input_file', type=str, help='Path to the input tab-separated file')
    parser.add_argument('output_file', type=str, help='Path to the output tab-separated file')
    parser.add_argument('year_range', type=str, help='Year range (e.g., 2000-2010)')

    # Parse command line arguments
    args = parser.parse_args()

    # Call the function to filter samples by year range and write to the output file
    filter_samples_by_year(args.input_file, args.output_file, args.year_range)

    print(f"Filtered samples for years {args.year_range} written to {args.output_file}")

if __name__ == '__main__':
    main()


