import pandas as pd
import glob

# Step 1: Find all CSV files in the directory
csv_files = glob.glob("/home/data/qs/data/reference_isoform/Homo_sapiens_GRCh38_cdna_all_modified_50mer_isoform_overlapping_count_otherisoform_2/*.csv")

# Step 2: Prepare an empty list to hold the result
max_values = []

# Step 3: Loop through each CSV file, calculate the maximum value of the 'overlapping_diff_frequency' column
for file in csv_files:
    df = pd.read_csv(file)

    # Find the maximum value in the 'overlapping_diff_frequency' column
    max_value = df['gene_set'].max()

    # Append the result as a tuple of file name and max value
    max_values.append((file, max_value))

# Step 4: Convert the list of tuples into a DataFrame
max_df = pd.DataFrame(max_values, columns=['File_Name', 'Max_Overlapping_Diff_Frequency'])

# Step 5: Output the result to a new CSV file
output_file = '/home/data/qs/data/reference_isoform/max_overlapping_diff_frequency_isoform_2.csv'
max_df.to_csv(output_file, index=False)

