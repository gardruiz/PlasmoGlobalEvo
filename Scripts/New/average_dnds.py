import sys

# Define lists of Asian and African countries
asia_countries = ['Cambodia', 'Thailand', 'Myanmar', 'Laos', 'Vietnam']
africa_countries = ['Ghana', 'Mauritania', 'Senegal', 'Gabon', 'Cameroon', 'Guinea', 'Kenya', 'Tanzania']

def calculate_averages(input_file):
    # Initialize sums and counts for dN and dS separately
    sum_dN_aa, sum_dS_aa, count_aa = 0, 0, 0
    sum_dN_af, sum_dS_af, count_af = 0, 0, 0
    sum_dN_asia, sum_dS_asia, count_asia = 0, 0, 0

    # Read the input file containing comparisons
    with open(input_file, "r") as file:
        for line in file:
            # Split the line into country1, country2, dN, and dS
            country1, country2, dN, dS = line.strip().split()
            
            # Convert dN and dS to floats
            dN, dS = float(dN), float(dS)
            
            # Check if both countries are from Asia
            if country1 in asia_countries and country2 in asia_countries:
                sum_dN_aa += dN
                sum_dS_aa += dS
                count_aa += 1
            # Check if both countries are from Africa
            elif country1 in africa_countries and country2 in africa_countries:
                sum_dN_af += dN
                sum_dS_af += dS
                count_af += 1
            # Check if one country is from Asia and the other from Africa
            elif (country1 in asia_countries and country2 in africa_countries) or (country1 in africa_countries and country2 in asia_countries):
                sum_dN_asia += dN
                sum_dS_asia += dS
                count_asia += 1

    # Calculate averages for dN and dS separately
    average_dN_aa = sum_dN_aa / count_aa if count_aa != 0 else 0
    average_dS_aa = sum_dS_aa / count_aa if count_aa != 0 else 0
    average_dN_af = sum_dN_af / count_af if count_af != 0 else 0
    average_dS_af = sum_dS_af / count_af if count_af != 0 else 0
    average_dN_asia = sum_dN_asia / count_asia if count_asia != 0 else 0
    average_dS_asia = sum_dS_asia / count_asia if count_asia != 0 else 0

    # Calculate average dN/dS
    average_dN_dS_aa = average_dN_aa / average_dS_aa if average_dS_aa != 0 else 0
    average_dN_dS_af = average_dN_af / average_dS_af if average_dS_af != 0 else 0
    average_dN_dS_asia = average_dN_asia / average_dS_asia if average_dS_asia != 0 else 0

    # Write averages to the output file
    with open("averages_output.txt", "w") as outfile:
        outfile.write("Africa Africa {:.6f} {:.6f} {:.6f}\n".format(average_dN_af, average_dS_af, average_dN_dS_af))
        outfile.write("Africa Asia {:.6f} {:.6f} {:.6f}\n".format(average_dN_asia, average_dS_asia, average_dN_dS_asia))
        outfile.write("Asia Asia {:.6f} {:.6f} {:.6f}\n".format(average_dN_aa, average_dS_aa, average_dN_dS_aa))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    calculate_averages(input_file)

