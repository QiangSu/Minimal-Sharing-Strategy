import os
import pandas as pd

# Define the base directory
base_dir = '/home/data/qs/Ndufs4_data/N2002874_ZJX_80-451580508/combined/trimmed_data/STAR_data/cufflinks_result/'

# Traverse the base directory and its subdirectories
for subdir, dirs, files in os.walk(base_dir):
    # Check if 'isoforms.fpkm_tracking' exists in the current subdirectory
    if 'isoforms.fpkm_tracking' in files:
        file_path = os.path.join(subdir, 'isoforms.fpkm_tracking')

        # Load the file into a pandas DataFrame
        df = pd.read_csv(file_path, sep='\t')

        # Sort the DataFrame based on the 'tracking_id' column
        df_sorted = df.sort_values(by='tracking_id')

        # Define a new file path by appending '_sorted' to the original file name
        new_file_path = os.path.join(subdir, 'isoforms_sorted.fpkm_tracking')

        # Output the sorted DataFrame to the new file
        df_sorted.to_csv(new_file_path, sep='\t', index=False)

        print(f"Sorted file saved as: {new_file_path}")

print("Processing complete.")
