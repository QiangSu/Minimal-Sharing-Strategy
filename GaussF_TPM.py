import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
import argparse
import re
from tqdm import tqdm  # Import tqdm for the progress bar

# Define a function to judge the suitability of the mean before fitting
def suitable_criteria_for_GC(xc_fitted_local, w_fitted_local):
    return xc_fitted_local > 0.5 * w_fitted_local

# Function to calculate total_normalized_kmer_count across all CSV files
def sum_normalized_kmer_counts(input_directory):
    total_normalized_kmer_count = 0

    # List all CSV files in the input directory
    csv_files = [f for f in os.listdir(input_directory) if f.endswith('_normalized.csv')]

    # Loop through each CSV file
    for csv_file in csv_files:
        file_path = os.path.join(input_directory, csv_file)

        # Read the CSV file into a DataFrame
        df = pd.read_csv(file_path)

        # Sum the 'Normalized_K-mer_Count' column
        file_normalized_kmer_sum = df['Normalized_K-mer_Count'].sum()

        # Add the file sum to the total sum
        total_normalized_kmer_count += file_normalized_kmer_sum

    return total_normalized_kmer_count

# Parse command line arguments
parser = argparse.ArgumentParser(description="Analyze GC content and fit Gaussian CDF.")
parser.add_argument('--input', type=str, required=True, help="Path to the input folder containing the CSV files.")
parser.add_argument('--output', type=str, required=True, help="Path and name of the output file to save the results.")
parser.add_argument('--threshold', type=int, default=10, help="Minimum number of k-mers required for fitting. Default is 10.")
args = parser.parse_args()

# GC content calculation function
def calculate_gc_content(kmer):
    gc_count = kmer.count('G') + kmer.count('C')
    total_bases = len(kmer)
    gc_content_percent = (gc_count / total_bases) * 100
    return round(gc_content_percent, 2)

# Gaussian CDF definition
def gaussian_cdf(x, A0, A, xc, w):
    return A0 + A * norm.cdf((x - xc) / w)

# Gaussian CDF with fixed parameters for Normalized K-mer Count
def gaussian_cdf_fixed(x, A0, A, xc_fixed, w_fixed):
    return A0 + A * norm.cdf((x - xc_fixed) / w_fixed)

# Gaussian CDF with fixed parameters for Count
def gaussian_cdf_fixed_count(x, A0, A, xc_fixed, w_fixed):
    return A0 + A * norm.cdf((x - xc_fixed) / w_fixed)

# Extract gene name and transcript ID from filename

def extract_gene_transcript_id(filename):
    # Match for mouse (Mus) data with transcript ID and hyphen in gene name
    match_mus = re.search(r'([\w-]+)_(ENSMUST[0-9]+)(\.\d+)?_kmers', filename)
    # Match for human (Homo sapiens) data with transcript ID and hyphen in gene name
    match_human = re.search(r'([\w-]+)_(ENST[0-9]+)(\.\d+)?_kmers', filename)
    # Match for files with no transcript ID (only gene name with hyphen)
    match_no_transcript = re.search(r'([\w-]+)_kmers', filename)

    if match_mus:
        gene_name = match_mus.group(1)
        transcript_id = match_mus.group(2)
        version = match_mus.group(3) if match_mus.group(3) else ""
        return gene_name, transcript_id + version
    elif match_human:
        gene_name = match_human.group(1)
        transcript_id = match_human.group(2)
        version = match_human.group(3) if match_human.group(3) else ""
        return gene_name, transcript_id + version
    elif match_no_transcript:
        gene_name = match_no_transcript.group(1)
        return gene_name, "-"
    else:
        return "Unknown", "Unknown"


# Calculate total_normalized_kmer_count across all CSV files in the input directory
total_normalized_kmer_count = sum_normalized_kmer_counts(args.input)

# List to store the results
results = []

# Get the list of files and sort them
files = sorted(os.listdir(args.input))

# Loop through each file in the directory with a progress bar
for filename in tqdm(files, desc="Processing files", unit="file"):
    if filename.endswith("merged_normalized.csv"):
        filepath = os.path.join(args.input, filename)

        # Extract gene name and transcript ID from the filename
        gene_name, transcript_id = extract_gene_transcript_id(filename)

        # Read the CSV file into a DataFrame
        df = pd.read_csv(filepath)

        # Ensure the 'kmer' column exists
        if 'kmer' not in df.columns:
            print(f"'kmer' column is missing in file {filename}. Skipping this file.")
            continue

        # Add the transcript length value extracted from the first row of the DataFrame
        transcript_length = df.at[0, 'Transcript_Length']

        # Calculate and add GC content to the DataFrame
        df['GC_Content'] = df['kmer'].apply(calculate_gc_content)

        # Group by GC content and sum frequencies
        gc_content_data = df.groupby('GC_Content').agg({
            'Local_Frequency': 'sum',
            'Normalized_K-mer_Count': 'sum',
            'Count': 'sum'
        }).reset_index()

        # Add the first row's Global_Frequency and Present_in_Transcripts to every result
        global_frequency = df.at[0, 'Global_Frequency']
        present_in_transcripts = df.at[0, 'Present_in_Transcripts']

        # Check if there are at least args.threshold distinct GC contents
        if len(gc_content_data['GC_Content'].unique()) < args.threshold:
            sum_normalized_kmer_count = gc_content_data['Normalized_K-mer_Count'].sum()
            # Normalize the sum by multiplying by 1000000/total_normalized_kmer_count
            normalized_sum = sum_normalized_kmer_count * 1000000 / total_normalized_kmer_count
            sum_count = gc_content_data['Count'].sum()
            results.append({
                'File': filename,
                'Gene_Name': gene_name,
                'Transcript_ID': transcript_id,
                'Global_Frequency': global_frequency,
                'Present_in_Transcripts': present_in_transcripts,
                'Transcript_Length': transcript_length,
                'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(normalized_sum),
                'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(sum_count),
                'Fixed Mean (xc)': 'N/A',
                'Fixed Standard Deviation (w)': 'N/A',
                'Report': 'Not enough distinct GC contents'
            })
            continue

        # Calculate cumulative sums
        gc_content_data_sorted = gc_content_data.sort_values(by='GC_Content')
        gc_content_data_sorted['Cumulative_Local_Frequency'] = gc_content_data_sorted['Local_Frequency'].cumsum()
        gc_content_data_sorted['Cumulative_Normalized_Count'] = gc_content_data_sorted['Normalized_K-mer_Count'].cumsum()
        gc_content_data_sorted['Cumulative_Count'] = gc_content_data_sorted['Count'].cumsum()

        # Get the data for fitting
        x_data = gc_content_data_sorted['GC_Content']
        y_data_local = gc_content_data_sorted['Cumulative_Local_Frequency']
        y_data_normalized = gc_content_data_sorted['Cumulative_Normalized_Count']
        y_data_count = gc_content_data_sorted['Cumulative_Count']

        # Fit the Gaussian CDF to Cumulative Local Frequency to get initial xc and w
        initial_guesses_local = [min(y_data_local), max(y_data_local) - min(y_data_local), x_data.mean(), x_data.std()]
        try:
            popt_local, pcov_local = curve_fit(gaussian_cdf, x_data, y_data_local, p0=initial_guesses_local)
            A0_fitted_local, A_fitted_local, xc_fitted_local, w_fitted_local = popt_local

            # Check if the fitted mean and standard deviation meet the suitability criteria
            if suitable_criteria_for_GC(xc_fitted_local, w_fitted_local):
                # Additional fitting for Normalized K-mer Count and Count using fixed xc and w
                initial_guesses_normalized = [min(y_data_normalized), max(y_data_normalized) - min(y_data_normalized)]
                popt_normalized, pcov_normalized = curve_fit(
                    lambda x, A0, A: gaussian_cdf_fixed(x, A0, A, xc_fitted_local, w_fitted_local),
                    x_data,
                    y_data_normalized,
                    p0=initial_guesses_normalized
                )
                A0_fitted_normalized, A_fitted_normalized = popt_normalized

                # Normalize A_fitted_normalized
                A_fitted_normalized = A_fitted_normalized * 1000000 / total_normalized_kmer_count

                initial_guesses_count = [min(y_data_count), max(y_data_count) - min(y_data_count)]
                popt_count, pcov_count = curve_fit(
                    lambda x, A0, A: gaussian_cdf_fixed_count(x, A0, A, xc_fitted_local, w_fitted_local),
                    x_data,
                    y_data_count,
                    p0=initial_guesses_count
                )
                A0_fitted_count, A_fitted_count = popt_count

                # Append successful fitting results
                results.append({
                    'File': filename,
                    'Gene_Name': gene_name,
                    'Transcript_ID': transcript_id,
                    'Global_Frequency': global_frequency,
                    'Present_in_Transcripts': present_in_transcripts,
                    'Transcript_Length': transcript_length,
                    'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(A_fitted_normalized),
                    'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(A_fitted_count),
                    'Fixed Mean (xc)': '{:.2f}'.format(xc_fitted_local),
                    'Fixed Standard Deviation (w)': '{:.2f}'.format(w_fitted_local),
                    'Report': 'OK'
                })
            else:
                raise ValueError("Unsuitable fit parameters: mean is not greater than standard deviation.")

        except (RuntimeError, ValueError) as e:
            error_message = str(e)
            # Handle fitting failures or unsuitable parameters. Potentially output sums as a fallback.
            sum_normalized_kmer_count = gc_content_data_sorted['Normalized_K-mer_Count'].sum()
            normalized_sum = sum_normalized_kmer_count * 1000000 / total_normalized_kmer_count
            sum_count = gc_content_data_sorted['Count'].sum()
            results.append({
                'File': filename,
                'Gene_Name': gene_name,
                'Transcript_ID': transcript_id,
                'Global_Frequency': global_frequency,
                'Present_in_Transcripts': present_in_transcripts,
                'Transcript_Length': transcript_length,
                'Sum or Fitted A (Abundance) for Normalized Count': '{:.2f}'.format(normalized_sum),
                'Sum or Fitted A (Abundance) for Count': '{:.2f}'.format(sum_count),
                'Fixed Mean (xc)': 'N/A',
                'Fixed Standard Deviation (w)': 'N/A',
                'Report': f'Unsuitable Fit - {error_message}'
            })

# Create results DataFrame
results_df = pd.DataFrame(results)

# Save to CSV file specified by the command-line argument
output_directory = os.path.dirname(args.output)  # Extract the directory part from the output path
os.makedirs(output_directory, exist_ok=True)  # Create the directory if it doesn't exist

results_df.to_csv(args.output, index=False)  # Save the dataframe
