import pandas as pd
import os

# Define the file names
files = [
    "musHD_STR-1_salmon_gcbiascorrect_quant.csv",
    "musHD_STR-2_salmon_gcbiascorrect_quant.csv",
    "musHD_STR-3_salmon_gcbiascorrect_quant.csv",
    "musWT_STR-1_salmon_gcbiascorrect_quant.csv",
    "musWT_STR-2_salmon_gcbiascorrect_quant.csv",
    "musWT_STR-3_salmon_gcbiascorrect_quant.csv"
]

# Initialize an empty DataFrame for the combined data
combined_df = pd.DataFrame()

for file in files:
    # Read the CSV file with tab separation
    df = pd.read_csv(file, delimiter='\t')

    # Initialize the combined DataFrame with required columns only if it is empty
    if combined_df.empty:
        combined_df = df[['Name', 'Length', 'EffectiveLength']].copy()

    # Extract the sample name from the file name, and append the read counts column to the combined DataFrame
    sample_name = os.path.basename(file).split('_salmon')[0]
    combined_df[sample_name + '_TPM'] = df['TPM']

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('/home/data/qs/HGC20231021002-0004/trimmed_data/musWT_HD_salmon_gcbiascorrect_combined_TPM.csv', index=False)
