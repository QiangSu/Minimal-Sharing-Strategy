from collections import defaultdict
import os
import csv
from tqdm import tqdm
from pyfaidx import Fasta
from multiprocessing import Pool

# Define constants
KMER_LENGTH = 50
OUTPUT_DIRECTORY = './parsed_exons_sequence_Homo_sapiens.GRCh38.112_50mer_exonID'
os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

# Step 1: Parse GTF file to retrieve exon coordinates and associated gene/transcript information
def parse_gtf(file_path):
    exon_info = []

    # Count total lines for progress bar
    with open(file_path, 'r') as file:
        total_lines = sum(1 for line in file if not line.startswith('#'))

    with open(file_path, 'r') as file:
        for line in tqdm(file, total=total_lines, desc="Parsing GTF"):
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'exon':
                chrom = fields[0]
                start = int(fields[3]) - 1  # Convert to 0-based indexing
                end = int(fields[4])
                strand = fields[6]
                # Parse attributes
                attributes = {}
                for attribute in fields[8].split(';'):
                    if attribute.strip():
                        key, value = attribute.strip().split(' ', 1)
                        attributes[key] = value.strip('"')
                gene_id = attributes.get('gene_id', '')
                transcript_id = attributes.get('transcript_id', '')
                exon_number = attributes.get('exon_number', '')
                exon_id = attributes.get('exon_id', '')
                exon_info.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_id': gene_id,
                    'transcript_id': transcript_id,
                    'exon_number': exon_number,
                    'exon_id': exon_id
                })
    return exon_info

# Step 2: Generate the reverse complement of a DNA sequence
def reverse_complement(seq):
    complement = str.maketrans('ATCGatcgNn', 'TAGCtagcNn')
    return seq.translate(complement)[::-1]

# Step 3: Extract exon sequence and format for FASTA output
def extract_exon_sequence(args):
    exon, fasta_path = args
    fasta = Fasta(fasta_path)
    chrom = exon['chrom']
    start = exon['start']
    end = exon['end']
    strand = exon['strand']
    sequence = fasta[chrom][start:end].seq
    if strand == '-':
        sequence = reverse_complement(sequence)
    header = f"{exon['gene_id']}|{exon['transcript_id']}|{exon['exon_id']}|{chrom}:{start+1}-{end}|{strand}"
    return header, sequence

# Step 4: Extract and save exon sequences with k-mers to CSV
def extract_and_save_exon_kmers(exon_info, fasta_path):
    # Prepare global k-mer statistics
    global_kmer_counts = defaultdict(int)
    kmer_transcript_sets = defaultdict(set)  # Stores which exons contain the kmer

    with Pool() as pool:
        args = [(exon, fasta_path) for exon in exon_info]
        results = pool.imap_unordered(extract_exon_sequence, args)

        # First pass: Build global k-mer statistics across all exons
        print("Counting global k-mer frequencies across all exons...")
        exon_sequences = []  # Store sequences temporarily for second pass
        for header, sequence in tqdm(results, total=len(exon_info), desc="Global K-mer Counting"):
            exon_sequences.append((header, sequence))
            exon_index = len(exon_sequences) - 1  # Track index for current exon
            for i in range(len(sequence) - KMER_LENGTH + 1):
                kmer = sequence[i:i + KMER_LENGTH]
                global_kmer_counts[kmer] += 1
                kmer_transcript_sets[kmer].add(exon_index)

    # Second pass: Write the k-mers for each individual exon
    print("Writing k-mers to CSV files for each exon...")
    for exon_index, (header, sequence) in tqdm(enumerate(exon_sequences), total=len(exon_sequences), desc="Writing CSVs"):
        output_csv_path = os.path.join(OUTPUT_DIRECTORY, sanitize_filename(header) + '_kmers.csv')

        # Local k-mer occurrences for this exon
        local_kmer_counts = defaultdict(int)
        for i in range(len(sequence) - KMER_LENGTH + 1):
            kmer = sequence[i:i + KMER_LENGTH]
            local_kmer_counts[kmer] += 1

        # Skip creating the file if there are no k-mers
        if not local_kmer_counts:
            continue

        with open(output_csv_path, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(['kmer', 'Local_Frequency', 'Global_Frequency', 'Present_in_Transcripts'])

            # Write the kmers, frequencies, and sets of exons where the kmer is present
            for kmer, local_freq in local_kmer_counts.items():
                global_freq = global_kmer_counts[kmer]
                # Filter out exons where the k-mer is not actually present
                transcripts_containing_kmer = [
                    exon_sequences[i][0] for i in kmer_transcript_sets[kmer]
                    if kmer in exon_sequences[i][1]
                ]
                if transcripts_containing_kmer:
                    transcripts_containing_kmer_str = ', '.join(transcripts_containing_kmer)
                    num_transcripts_containing_kmer = len(transcripts_containing_kmer_str.split(', '))

                    # Check if the global frequency matches the number of exons
                    if global_freq != num_transcripts_containing_kmer:
                        global_freq = num_transcripts_containing_kmer  # Correct the global frequency

                    csv_writer.writerow([kmer, local_freq, global_freq, transcripts_containing_kmer_str])

    print(f"K-mers for each exon have been saved as CSV files in the directory: {OUTPUT_DIRECTORY}")

# Utility function to sanitize file names
def sanitize_filename(header):
    return header.replace('|', '_')

# Main script execution
if __name__ == "__main__":
    # Paths to your GTF and FASTA files
    gtf_path = './Homo_sapiens.GRCh38.112.gtf'
    fasta_path = './Homo_sapiens.GRCh38.dna.primary_assembly.fa'

    # Parse GTF file to get exon information
    exon_info = parse_gtf(gtf_path)

    # Extract exon sequences and save k-mer information
    extract_and_save_exon_kmers(exon_info, fasta_path)
