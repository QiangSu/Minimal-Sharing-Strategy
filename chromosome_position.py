from Bio import SeqIO
import re

def find_all_sequences_in_fasta(fasta_file, sequence):
    matches = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        for match in re.finditer(sequence, str(record.seq)):
            matches.append((record.id, match.start(), match.end()))
    return matches

def find_gene_in_gtf(gtf_file, chrom, start, end):
    with open(gtf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[0] == chrom and fields[2] == 'gene':
                gene_start = int(fields[3])
                gene_end = int(fields[4])
                if gene_start <= start <= gene_end or gene_start <= end <= gene_end:
                    gene_info = fields[8]
                    gene_name_match = re.search('gene_name "([^"]+)"', gene_info)
                    if gene_name_match:
                        return gene_name_match.group(1)
    return None

fasta_file = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
gtf_file = 'Homo_sapiens.GRCh38.112.gtf'
sequence = 'TACAGGTCTTTGCGGATGTCCACGTCACACTTCATGATGGAGTTGAAGGT'

matches = find_all_sequences_in_fasta(fasta_file, sequence)
if matches:
    for chrom, start, end in matches:
        gene_name = find_gene_in_gtf(gtf_file, chrom, start, end)
        if gene_name:
            print(f"Sequence found in {chrom}:{start}-{end}, Gene: {gene_name}")
        else:
            print(f"Sequence found in {chrom}:{start}-{end}, but no gene found.")
else:
    print("Sequence not found in the FASTA file.")
