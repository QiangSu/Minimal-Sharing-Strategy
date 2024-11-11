import csv
import sys
import re
import os
import shutil
from tqdm import tqdm  # Import tqdm for progress bar

# Increase the CSV field size limit
csv.field_size_limit(sys.maxsize)

# Set the threshold value
threshold_value = 1

# Specify the input and output directories
input_dir = '/home/data/qs/data/reference_isoform/gene_specific_output_500mer_wholefasta_allgene_geneID//'  # Change this to your actual input directory
output_dir = '/home/data/qs/data/reference_isoform/gene_specific_output_500mer_wholefasta_allgene_geneID_overlapping/'  # Change this to your actual output directory
new_output_dir = '/home/data/qs/data/reference_isoform/gene_specific_output_500mer_wholefasta_allgene_geneID_overlapping_count_othergene/'  # Directory for files exceeding threshold

# Create the output directories if they do not exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(new_output_dir):
    os.makedirs(new_output_dir)

# Regex pattern to match filenames and extract part of the name
filename_regex = re.compile(r'^(.*?)(?:_ENS[\w\.]+)?_kmers\.csv$')

# Iterate over all CSV files in the input directory
csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
for input_file in tqdm(csv_files, desc="Processing CSV files"):  # Add progress bar to file processing loop
    # Check if the filename matches the expected pattern
    term_to_count_match = filename_regex.match(input_file)
    if term_to_count_match:
        term_to_count = term_to_count_match.group(1)

        # Generate the output file paths
        base_name, ext = os.path.splitext(input_file)
        output_file = f"{base_name}_with_frequency{ext}"
        output_path = os.path.join(output_dir, output_file)
        new_output_path = os.path.join(new_output_dir, output_file)

        # Read the input CSV file
        input_path = os.path.join(input_dir, input_file)
        with open(input_path, mode='r', newline='') as infile:
            reader = csv.DictReader(infile)
            fieldnames = reader.fieldnames + ['Samegene_Frequency', 'overlapping_diff_frequency', 'gene_set']
            rows = list(reader)

        # Process each row
        for row in rows:
            present_in_transcripts = row['Present_in_Transcripts']
            transcripts = present_in_transcripts.split(', ')

            # Calculate the Samegene_Frequency
            term_frequency = sum(1 for transcript in transcripts if transcript.split('|')[0] == term_to_count)
            row['Samegene_Frequency'] = term_frequency

            # Calculate the number of isoforms in Present_in_Transcripts
            total_isoforms = len(transcripts)

            # Check and update Global_Frequency if necessary
            global_frequency = int(row['Global_Frequency'])
            if global_frequency != total_isoforms:
                print(f"Updating Global_Frequency for {base_name} kmer as it doesn't match total isoforms.")
                global_frequency = total_isoforms
                row['Global_Frequency'] = global_frequency  # Update the Global_Frequency in the row

            # Calculate the overlapping_diff_frequency
            overlapping_diff_frequency = global_frequency - term_frequency
            row['overlapping_diff_frequency'] = overlapping_diff_frequency

            # Extract unique gene names
            gene_names = {transcript.split('|')[0] for transcript in transcripts}  # Use a set to get unique names
            row['gene_set'] = len(gene_names)  # Store the number of unique gene names

        # Compute the sum of 'overlapping_diff_frequency' column
        overlapping_diff_total = sum(int(row['overlapping_diff_frequency']) for row in rows)

        # Write the updated data to the output CSV file
        with open(output_path, mode='w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        print(f"Updated CSV file has been saved as {output_path}")

        # Check if the total overlapping_diff_frequency is greater than threshold_value
        if overlapping_diff_total > threshold_value:
            # Copy the file to the new directory
            shutil.copy(output_path, new_output_path)
            print(f"File {output_file} has been copied to {new_output_dir} as the total overlapping_diff_frequency ({overlapping_diff_total}) exceeds the threshold ({threshold_value}).")
    else:
        print(f"Skipping file {input_file} as it does not match expected pattern.")
