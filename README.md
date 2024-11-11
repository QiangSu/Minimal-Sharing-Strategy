# Minimal-Sharing-Strategy
Transcript-Specific Minimal-Shared Regions Refining Transcriptome Profiles in RNA-Seq
## 1.1 <span style="color:#4CAF50">exon_kmer_frequency_distribution.py</span> 
Exon k-mer Overlapping Frequency Analysis: The `<span style="color:#4CAF50">exon_kmer_frequency_distribution.py</span>` script performs a comprehensive analysis of k-mer occurrences at the exon level within genomic data. This Python tool extracts exon sequences, examines overlapping k-mer substrings, and outputs the results in CSV format, streamlining k-mer-based studies on exon structure. Initially, it parses a GTF file to gather details about each exon’s location, strand orientation, and identifiers—crucial elements for accurately pinpointing sequences in the genome. For exons located on the reverse strand, it calculates the reverse complement to maintain strand specificity. Leveraging a FASTA genome file, the script retrieves the exact nucleotide sequence for each exon, then divides it into overlapping k-mer segments. It computes both local k-mer frequencies (unique to each exon) and global frequencies (across all exons) and saves this data in CSV files, each documenting k-mer sequences, their frequencies, and associated transcripts for each exon. This systematic approach supports in-depth analyses of exon structures in genomic studies. 
kmer_frequency_overlapping_exons.py
The Python script ```kmer_frequency_overlapping_exons.py``` processes a set of CSV files containing exon and gene information, calculating and updating specific frequency values for each entry. It iterates through CSV files in a specified input directory, extracts exon and gene details, and performs calculations to determine the frequency differences between global and local occurrences. It also counts unique gene occurrences based on the data within each file.
## 1.2 kmer_frequency_allgenes.py 
Transcript k-mer Frequency Analysis: The script kmer_frequency_allgenes.py is designed to analyze k-mer occurrences within a FASTA file of transcript sequences, providing a comprehensive k-mer frequency analysis. It systematically extracts k-mers from each transcript, calculating both their local and global frequencies, and identifies the transcripts in which they appear. For each gene, the script generates CSV files that detail the specific k-mer sequence, its local frequency (the total count of the k-mer across all transcripts of the same gene), its global frequency (the number of distinct transcripts across all genes that contain the k-mer at least once), and a list of all transcripts from all genes that contain the k-mer. 
## 1.3 kmer_frequency_overlapping.py 
K-mer overlapping Analysis: The script, ```kmer_frequency_overlapping.py```, processes CSV files containing k-mer data from gene-specific outputs to analyze and update frequency information related to gene isoforms. It identifies files using a specific pattern and calculates the homologous gene overlapping frequency (Samegene_Frequency), which measures the presence of k-mers in transcripts of the same gene. Additionally, it computes the heterologous gene overlapping frequency (overlapping_diff_frequency) by finding the difference between the global frequency and the same-gene frequency. The script also identifies the gene_set, encompassing all related gene names involved in the analysis.
## 1.4 Max_Overlapping_Diff_Frequency.py 
Max Heterologous Gene Frequency: The script ```Max_Overlapping_Diff_Frequency.py``` outputs the maximum number of heterologous genes corresponding to each gene’s k-mers.
## 2. Isoform Level k-mer overlapping Analysis
## 2.1 kmer_frequency_distribution.py 
Isoform K-mer Frequency Distribution: The ```script kmer_frequency_distribution.py``` is designed to analyze and catalog k-mer sequences in transcript sequences from a FASTA file containing human cDNA data. It begins by reading the transcripts and storing them, then processes each transcript to calculate global k-mer frequencies across all transcripts and stores which transcripts contain each k-mer. In the next phase, the script generates a CSV file for each transcript, documenting each k-mer’s local frequency (within that transcript), global frequency (across all transcripts), and listing the specific transcripts where the k-mer is found. This allows for detailed comparison of k-mer presence and frequency across isoforms, aiding in genomic analysis and comparative studies. 
## 2.2 kmer_frequency_overlapping_isoforms.py 
Comprehensive Isoform K-mer Frequency Analysis: ```kmer_frequency_overlapping_isoforms.py``` processes a set of CSV files located in the specified directory to identify the maximum value of the overlapping_diff_frequency column in each file. It begins by searching for all CSV files within the target directory, then iterates through each file to read its contents and calculate the maximum value in the gene_set column (assumed to be a representation of overlapping differences). These maximum values are stored in a list, along with the respective file names.
## 2.3 Max_Overlapping_Diff_Frequency_isoform.py 
Max Heterologous Gene Frequency for Isoforms: ```Max_Overlapping_Diff_Frequency_isoform.py``` outputs the maximum number of heterologous isoforms corresponding to each isoform.
## 2.4 kmer_frequency_distribution_mini_shared.py 
The minimal-shared k-mer set filtering: ```kmer_frequency_distribution_mini_shared.py``` script processes a FASTA file containing nucleotide sequences, breaking each entire transcript isoform sequence down into k-mers (subsequences of a specified length). It calculates both local frequency (how often each k-mer appears within the given transcript isoform) and global frequency (how often each k-mer appears across all other transcript isoforms) for these k-mers. The script identifies k-mers with the minimum global frequency for each transcript, tracks which transcripts contain each k-mer, and writes this information into individual CSV files for each transcript. The script employs command-line arguments to specify input and output paths and uses tqdm for progress tracking.
## 2.5 kmer_count_v2_20240803.py 
k-mer counting from sequencing data: This Python script ```kmer_count_v2_20240803.py``` is designed for processing large FASTQ files, specifically for counting occurrences of k-mers (subsequences of length k) in sequencing data and matching those counts against a list of k-mers from CSV files. It operates by leveraging multithreading for efficient parallel processing, where a producer function reads the FASTQ file in chunks and extracts sequences, while multiple consumer threads process these sequences into k-mers and count their frequencies. The script uses a Queue to manage communication between the producer and consumers, ensuring memory efficiency. After counting k-mers, it matches the counts with a set of specified k-mers from input CSV files, generating output CSVs that contain the frequency of each k-mer, making it suitable for large-scale genomic data analysis.
## 2.6 merge_normalizing_count_TPM.py 
Merging and normalizing k-mer counting data: This Python script ```merge_normalizing_count_TPM.py``` is designed to merge k-mer count data with a reference k-mer dataset and calculate normalized k-mer counts for each entry. Given a directory of k-mer reference CSV files and a corresponding directory of k-mer count files, it processes each pair of files by merging them on the k-mer column, calculates a normalized frequency for each k-mer based on a specified read length and k-mer size, and saves the resulting merged data into an output directory. The script uses command-line arguments to specify input directories, output directory, read length, and k-mer length, making it adaptable for various data sets and sequencing configurations.
## 2.7 GaussF_TPM.py 
The transcript isoform quantification: This Python script ```GaussF_TPM.py``` is designed to analyze GC content in k-mer data across multiple gene transcript files and fit a Gaussian cumulative distribution function (CDF) to model the data. It processes multiple CSV files in a specified input directory, where each file contains k-mer data for different genes or transcripts. The script calculates the GC content for each k-mer, then fits a Gaussian CDF to various cumulative distributions (such as local frequency, normalized count, and absolute count) within the data. It extracts essential information like gene names and transcript IDs from filenames, applies fitting criteria to ensure meaningful results, and summarizes each file’s fitting parameters in an output CSV file. This streamlined approach enables comprehensive GC content analysis with robust error handling and suitability checks to ensure the reliability of Gaussian fitting for normalized k-mer data.
## 2.8 pipeline_isoform_quantification.sh
Integrated pipeline for merging steps 2.5–2.7 to output isoform tpm quantification: This ```script pipeline_isoform_quantification.sh``` automates the process of k-mer counting, merging, and quantification for isoform-specific analysis on a set of trimmed FASTQ files specified in a sample list. The script performs three main tasks: (1) counts isoform-specific k-mers in each sample’s FASTQ file using a specified k-mer counting script, (2) merges the k-mer counts with reference data for normalization, and (3) applies Gaussian filtering for GC-content-based isoform quantification. It supports parallel processing with a defined maximum of concurrent jobs, enhancing efficiency for large datasets. Intermediate files are managed carefully, with directories created as needed and temporary data optionally removed after processing. This setup provides a streamlined approach for high-throughput isoform quantification.
