import pandas as pd
import glob
import os

# Directory containing the CSV files
csv_dir = '/home/data/qs/Ndufs4_data/N2002873_ZJX_80-451570558/combined/trimmed_data/1P_kmer_count_merg_result/'

# List all CSV files in the directory
csv_files = glob.glob(os.path.join(csv_dir, 'results_file_GC_*.csv'))

# Sort the CSV files by their names
csv_files.sort()

# Initialize an empty DataFrame to store the combined results
combined_df = pd.DataFrame()

# Iterate over each sorted CSV file
for file in csv_files:
    # Read the CSV file
    df = pd.read_csv(file)

    # Extract the required columns
    df_subset = df[['File', 'Gene_Name', 'Transcript_ID', 'Global_Frequency', 'Present_in_Transcripts', 'Transcript_Length', 'Report']]

    # Extract the Sum or Fitted A (Abundance) for Count column
    abundance_column = df[['Sum or Fitted A (Abundance) for Normalized Count']]

    # Get the base file name without the directory and extension
    base_filename = os.path.basename(file).replace('.csv', '')

    # Rename the abundance column to the file name
    abundance_column.columns = [base_filename]

    # Concatenate the extracted columns and the abundance column
    df_combined = pd.concat([df_subset, abundance_column], axis=1)

    # Merge with the combined DataFrame on common columns
    if combined_df.empty:
        combined_df = df_combined
    else:
        combined_df = pd.merge(combined_df, df_combined, on=['File', 'Gene_Name', 'Transcript_ID', 'Global_Frequency', 'Present_in_Transcripts', 'Transcript_Length', 'Report'], how='outer')

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('Ndufs4_combined_GaussF_TPM.csv', index=False)
print('Combined results saved to dufs4_combined_GaussF_TPM.csv')
