from collections import defaultdict
import csv
import os
from tqdm import tqdm  # Importing tqdm for progress bar

def read_fasta(file_path):
    headers = []
    sequences = []
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                sequence = ""
                headers.append(line[1:])  # Remove the '>' character
            else:
                sequence += line
        if sequence:
            sequences.append(sequence)
    return headers, sequences

def sanitize_filename(header):
    return header.replace('|', '_')

# Assumes that 'transcript_headers' contains the associated names for the transcripts
transcript_headers, transcripts = read_fasta('./Homo_sapiens.GRCh38.cdna.all.modified.geneID.fasta')

output_directory = './Homo_sapiens_GRCh38_cdna_all_modified_500mer_isoform-geneID'
os.makedirs(output_directory, exist_ok=True)

kmer_length = 500
global_kmer_counts = defaultdict(int)
kmer_transcript_sets = defaultdict(set)  # Stores which transcripts contain the kmer

# First pass: build global kmer statistics across all transcripts
print("Counting global k-mer frequencies across all transcripts...")
for isoform_index, sequence in tqdm(enumerate(transcripts), total=len(transcripts), desc="Global K-mer Counting"):
    for i in range(len(sequence) - kmer_length + 1):
        kmer = sequence[i:i + kmer_length]
        global_kmer_counts[kmer] += 1
        kmer_transcript_sets[kmer].add(isoform_index)  # Ensure that all isoforms with this kmer are added

# Second pass: write the kmers for each individual transcript
print("Writing k-mers to CSV files for each transcript...")
for isoform_index, header in tqdm(enumerate(transcript_headers), total=len(transcript_headers), desc="Writing CSVs"):
    output_csv_path = os.path.join(output_directory, sanitize_filename(header) + '_kmers.csv')

    # Local kmer occurrences for this transcript
    local_kmer_counts = defaultdict(int)

    for i in range(len(transcripts[isoform_index]) - kmer_length + 1):
        kmer = transcripts[isoform_index][i:i + kmer_length]
        local_kmer_counts[kmer] += 1

    # Skip creating the file if there are no k-mers
    if not local_kmer_counts:
        continue

    with open(output_csv_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['kmer', 'Local_Frequency', 'Global_Frequency', 'Present_in_Transcripts'])

        # Write the kmers, frequencies, and sets of transcripts where the kmer is present
        for kmer, local_freq in local_kmer_counts.items():
            global_freq = global_kmer_counts[kmer]
            # Filter out transcripts where the k-mer is not actually present
            transcripts_containing_kmer = [
                transcript_headers[i] for i in kmer_transcript_sets[kmer]
                if kmer in transcripts[i]
            ]
            if transcripts_containing_kmer:
                # Join the headers for transcripts containing the kmer
                transcripts_containing_kmer_str = ', '.join(transcripts_containing_kmer)

                # Count the number of transcripts by counting commas + 1 or splitting by ','
                num_transcripts_containing_kmer = len(transcripts_containing_kmer_str.split(', '))

                # Check if the global frequency matches the number of transcripts
                if global_freq != num_transcripts_containing_kmer:
                    global_freq = num_transcripts_containing_kmer  # Correct the global frequency

                csv_writer.writerow([kmer, local_freq, global_freq, transcripts_containing_kmer_str])

print(f"Kmers for each transcript have been saved as CSV files in the directory: {output_directory}")
