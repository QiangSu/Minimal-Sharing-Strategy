import pandas as pd
import os

# Define the file names
files = [
    "/home/data/qs/HGC20231021002-0004/trimmed_data/RSEM_result/musHD_STR-1.isoforms.results.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/RSEM_result/musHD_STR-2.isoforms.results.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/RSEM_result/musHD_STR-3.isoforms.results.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/RSEM_result/musWT_STR-1.isoforms.results.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/RSEM_result/musWT_STR-2.isoforms.results.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/RSEM_result/musWT_STR-3.isoforms.results.csv"
]

# Initialize an empty DataFrame for the combined data
combined_df = pd.DataFrame()

# Define the data types for the columns
dtype_dict = {
    'transcript_id': str,
    'effective_length': float,
    'gene_id': str,
}

for file in files:
    # Read the CSV file with tab separation and specified data types
    df = pd.read_csv(file, delimiter='\t', dtype=dtype_dict)

    # Initialize the combined DataFrame with required columns only if it is empty
    if combined_df.empty:
        combined_df = df[['transcript_id', 'effective_length', 'gene_id']].copy()

    # Extract the sample name from the file name
    base_name = os.path.basename(file)

    # Remove ".isoforms.results.csv" from the filename to isolate the sample name
    sample_name = base_name.replace(".isoforms.results.csv", "")

    # Append the coverage column to the combined DataFrame
    combined_df[sample_name + '_TPM'] = df['TPM']

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('musWT_HD_RSEM_nobiascorrect_combined_est_TPM.csv', index=False)
