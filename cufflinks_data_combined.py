import pandas as pd
import os

# Define the file names
files = [
    "/home/data/qs/HGC20231021002-0004/trimmed_data/cufflinks_result_HD/cufflinks_result_nonbiascorrect_musHD_STR-1/sorted_isoforms.fpkm_tracking_cufflinks_result_nonbiascorrect_musHD_STR-1.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/cufflinks_result_HD/cufflinks_result_nonbiascorrect_musHD_STR-2/sorted_isoforms.fpkm_tracking_cufflinks_result_nonbiascorrect_musHD_STR-2.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/cufflinks_result_HD/cufflinks_result_nonbiascorrect_musHD_STR-3/sorted_isoforms.fpkm_tracking_cufflinks_result_nonbiascorrect_musHD_STR-3.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/cufflinks_result_HD/cufflinks_result_nonbiascorrect_musWT_STR-1/sorted_isoforms.fpkm_tracking_cufflinks_result_nonbiascorrect_musWT_STR-1.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/cufflinks_result_HD/cufflinks_result_nonbiascorrect_musWT_STR-2/sorted_isoforms.fpkm_tracking_cufflinks_result_nonbiascorrect_musWT_STR-2.csv",
    "/home/data/qs/HGC20231021002-0004/trimmed_data/cufflinks_result_HD/cufflinks_result_nonbiascorrect_musWT_STR-3/sorted_isoforms.fpkm_tracking_cufflinks_result_nonbiascorrect_musWT_STR-3.csv"
]

# Initialize an empty DataFrame for the combined data
combined_df = pd.DataFrame()

for file in files:
    # Read the CSV file with tab separation
    df = pd.read_csv(file, delimiter='\t')

    # Initialize the combined DataFrame with required columns only if it is empty
    if combined_df.empty:
        combined_df = df[['tracking_id', 'gene_id', 'gene_short_name', 'length']].copy()

    # Extract the sample name from the file name
    base_name = os.path.basename(file)
    sample_name_parts = base_name.split('_')[6:8]
    sample_name = '_'.join(sample_name_parts).replace('.csv', '')

    # Append the coverage column to the combined DataFrame
    combined_df[sample_name + '_FPKM'] = df['FPKM']

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('/home/data/qs/HGC20231021002-0004/trimmed_data/musWT_HD_cufflinks_nobiascorrect_combined_STR_FPKM.csv', index=False)
