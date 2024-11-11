import pandas as pd
import os

# Define the file names
files = [
    "/home/data/qs/HGC20231021002-0004/trimmed_data/musHD_STR-1_kallisto/abundance_musHD_STR-1_kallisto.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/musHD_STR-2_kallisto/abundance_musHD_STR-2_kallisto.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/musHD_STR-3_kallisto/abundance_musHD_STR-3_kallisto.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/musWT_STR-1_kallisto/abundance_musWT_STR-1_kallisto.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/musWT_STR-2_kallisto/abundance_musWT_STR-2_kallisto.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/musWT_STR-3_kallisto/abundance_musWT_STR-3_kallisto.csv"
]

# Initialize an empty DataFrame for the combined data
combined_df = pd.DataFrame()

for file in files:
    # Read the CSV file with tab separation
    df = pd.read_csv(file, delimiter='\t')

    # Initialize the combined DataFrame with required columns only if it is empty
    if combined_df.empty:
        combined_df = df[['target_id', 'length', 'eff_length']].copy()

    # Extract the sample name from the file name
    base_name = os.path.basename(file)
    sample_name_parts = base_name.split('_')[1:4]
    sample_name = '_'.join(sample_name_parts).replace('.csv', '')

    # Append the coverage column to the combined DataFrame
    combined_df[sample_name + '_tpm'] = df['tpm']

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('musWT_HD_kallisto_nobiascorrect_combined_STR_tpm.csv', index=False)
