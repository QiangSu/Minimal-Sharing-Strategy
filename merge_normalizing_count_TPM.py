import pandas as pd
import os
import glob
import argparse
from tqdm import tqdm

# Updated normalization function (removed total_read_count)
def normalize_frequency(kmer_count, total_csv_kmers, read_length, k):
    multiplier = 2 * 1000 / (total_csv_kmers * (read_length - k + 1))
    return kmer_count * multiplier

# Function to merge k-mer count files to original k-mer CSV and add normalized count
def merge_kmer_counts_to_kmers(kmers_filepath, kmer_counts_filepath, total_csv_kmers, read_length, k, output_filepath):
    # Load the original kmer data and kmer counts
    kmers_df = pd.read_csv(kmers_filepath)
    kmer_counts_df = pd.read_csv(kmer_counts_filepath)

    # Merge the dataframes on the 'kmer' column
    merged_df = kmers_df.merge(kmer_counts_df, left_on='kmer', right_on='K-mer')

    # Drop the redundant 'K-mer' column from the merged dataframe
    merged_df.drop('K-mer', axis=1, inplace=True)

    # Add normalized k-mer counts
    merged_df['Normalized_K-mer_Count'] = merged_df['Count'].apply(normalize_frequency, args=(total_csv_kmers, read_length, k))

    # Add the 'Mini_Shared_Length' column, which contains the total number of kmers for each row
    merged_df['Transcript_Length'] = total_csv_kmers

    # Write the merged data to a new CSV file
    merged_df.to_csv(output_filepath, index=False)

# Main function
def main(kmer_reference_directory, kmer_counts_directory, output_directory, read_length, k):
    # Ensure the output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Process all the files
    kmer_files = glob.glob(os.path.join(kmer_reference_directory, "*_kmers.csv"))
    for kmer_filepath in tqdm(kmer_files, desc="Processing kmer files"):
        # Calculate the number of k-mers in the current CSV file
        total_csv_kmers = pd.read_csv(kmer_filepath).shape[0]

        # Compute the filename of the corresponding kmer counts file
        basename = os.path.splitext(os.path.basename(kmer_filepath))[0]
        kmer_counts_filename = f"{basename}_counts.csv"
        kmer_counts_filepath = os.path.join(kmer_counts_directory, kmer_counts_filename)

        # Define the output file path in the output directory
        output_filename = f"{basename}_merged_normalized.csv"
        output_filepath = os.path.join(output_directory, output_filename)

        # Merge the kmer counts into the kmer file and add normalized counts
        merge_kmer_counts_to_kmers(kmer_filepath, kmer_counts_filepath, total_csv_kmers, read_length, k, output_filepath)

# Run script with command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge k-mer counts into original k-mer CSV files and add normalized counts.')
    parser.add_argument('--kmer_reference_directory', type=str, required=True, help='Directory containing the k-mer reference CSV files.')
    parser.add_argument('--kmer_counts_directory', type=str, required=True, help='Directory containing the k-mer counts CSV files.')
    parser.add_argument('--output_directory', type=str, required=True, help='Output directory for storing merged CSV files.')
    parser.add_argument('--read_length', type=int, default=150, help='Read length of FASTQ sequences.')
    parser.add_argument('--k', type=int, default=50, help='K-mer length.')
    args = parser.parse_args()
    main(args.kmer_reference_directory, args.kmer_counts_directory, args.output_directory, args.read_length, args.k)
