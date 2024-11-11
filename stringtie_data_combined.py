import pandas as pd
import os

# Define the file names
files = [
    "/home/data/qs/HGC20231021002-0004/trimmed_data/stringtie_result_musHD_STR-1/t_data_stringtie_result_musHD_STR-1.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/stringtie_result_musHD_STR-2/t_data_stringtie_result_musHD_STR-2.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/stringtie_result_musHD_STR-3/t_data_stringtie_result_musHD_STR-3.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/stringtie_result_musWT_STR-1/t_data_stringtie_result_musWT_STR-1.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/stringtie_result_musWT_STR-2/t_data_stringtie_result_musWT_STR-2.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/stringtie_result_musWT_STR-3/t_data_stringtie_result_musWT_STR-3.csv"
]

# Initialize an empty DataFrame for the combined data
combined_df = pd.DataFrame()

# Define the data types for the columns
dtype_dict = {
    't_name': str,
    'length': float,
    'gene_id': str,
    'cov': float
}

for file in files:
    # Read the CSV file with tab separation and specified data types
    df = pd.read_csv(file, delimiter='\t', dtype=dtype_dict)

    # Initialize the combined DataFrame with required columns only if it is empty
    if combined_df.empty:
        combined_df = df[['t_name', 'length', 'gene_id']].copy()

    # Extract the sample name from the file name
    base_name = os.path.basename(file)
    sample_name_parts = base_name.split('_')[2:6]
    sample_name = '_'.join(sample_name_parts).replace('.csv', '')

    # Append the coverage column to the combined DataFrame
    combined_df[sample_name + '_FPKM'] = df['FPKM']

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('musWT_HD_stringtie_nobiascorrect_combined_est_STR_FPKM.csv', index=False)
