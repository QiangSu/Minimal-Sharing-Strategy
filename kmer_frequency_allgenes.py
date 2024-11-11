from collections import defaultdict
import csv
import os
from tqdm import tqdm  # Import tqdm for progress bars

def read_fasta(file_path):
    headers = []
    sequences = []
    sequence = ""
    with open(file_path, 'r') as file:
        for line in tqdm(file, desc="Reading FASTA file"):
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                sequence = ""
                headers.append(line[1:])  # Remove the '>' character
            else:
                sequence += line
        sequences.append(sequence)  # Don't forget the last one
    return headers, sequences

def parse_gene_name(header):
    return header.split('|')[0]  # Assumes the gene name is the first part before the '|'

output_directory = './gene_specific_output_500mer_wholefasta_allgene_geneID'
os.makedirs(output_directory, exist_ok=True)

kmer_length = 500

# Global data structures to hold kmer counts and mapping
gene_kmer_data = defaultdict(lambda: defaultdict(int))  # Gene -> Kmer -> Global Frequency (per gene)
gene_transcript_mapping = defaultdict(lambda: defaultdict(set))  # Gene -> Kmer -> Set of Transcripts
transcript_kmer_data = defaultdict(lambda: defaultdict(int))  # Transcript -> Kmer -> Local Frequency

# New structure to track global frequency across all transcripts
global_kmer_data = defaultdict(int)  # Kmer -> Global Frequency across all genes

# New structure to track kmer to all transcripts mapping
kmer_transcripts = defaultdict(set)  # Kmer -> Set of all Transcripts containing it

file_path = './Homo_sapiens.GRCh38.cdna.all.modified.geneID.fasta'
transcript_headers, transcripts = read_fasta(file_path)

# First, populate kmer data for global and local frequency
for header, sequence in tqdm(zip(transcript_headers, transcripts), desc="Processing transcripts", total=len(transcript_headers)):
    gene_name = parse_gene_name(header)
    seen_kmers = set()  # Track which kmers have been seen for this transcript
    for j in range(len(sequence) - kmer_length + 1):
        kmer = sequence[j:j + kmer_length]

        # Global kmer count for the gene
        gene_kmer_data[gene_name][kmer] += 1

        # Transcript kmer count (Local Frequency)
        transcript_kmer_data[header][kmer] += 1

        # Track which transcripts of the gene have this kmer
        gene_transcript_mapping[gene_name][kmer].add(header)

        # Track which transcripts (from all genes) have this kmer
        kmer_transcripts[kmer].add(header)

        # If the kmer hasn't been seen yet in this transcript, increment global frequency
        if kmer not in seen_kmers:
            global_kmer_data[kmer] += 1
            seen_kmers.add(kmer)  # Mark kmer as seen for this transcript

# Now, write the output CSV files with Global and Local Frequency based on all transcripts
for gene, kmers in tqdm(gene_kmer_data.items(), desc="Writing CSV files"):
    with open(os.path.join(output_directory, f'{gene}_kmers.csv'), 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['kmer', 'Local_Frequency', 'Global_Frequency', 'Present_in_Transcripts'])

        for kmer, gene_specific_global_frequency in kmers.items():
            # Calculate the sum of local frequencies across transcripts of the same gene
            local_frequencies = sum(
                transcript_kmer_data[transcript][kmer]
                for transcript in gene_transcript_mapping[gene][kmer]
            )

            # Global frequency from the global_kmer_data (number of distinct transcripts containing the kmer)
            global_frequency = global_kmer_data[kmer]

            # List of all transcripts sharing this kmer (from all genes)
            transcripts = ', '.join(kmer_transcripts[kmer])

            # Write kmer data with summed Local Frequency and Global Frequency (across all genes)
            csv_writer.writerow([kmer, local_frequencies, global_frequency, transcripts])

print(f"K-mer mapping with global frequency across all transcripts has been saved in specific CSV files in the directory: {output_directory}")
