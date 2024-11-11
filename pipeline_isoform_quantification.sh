#!/bin/bash

# Define the path to the sample list file
sample_list="/home/data/qs/Ndufs4_data/N2002873_ZJX_80-451570558/combined/trimmed_data/Ndufs4_sample_list.txt"

# Define the paths to the scripts and directories
kmer_count_script="/home/data/qs/mus_str.sn/HGC20231021002-0003-2/kmer_count_v2_20240803.py"
merge_script="/home/data/qs/colorec_data/merge_normalizing_count_TPM.py"
gaussf_script="/home/data/qs/colorec_data/GaussF_TPM.py"
reference_dir="/home/data/qs/data/reference_isoform/mus_ref/50mer_Mus_musculus_GRCc39_isoform_mini_shared_filter"
input_fastq_dir="/home/data/qs/Ndufs4_data/N2002873_ZJX_80-451570558/combined/trimmed_data"
output_dir="/home/data/qs/Ndufs4_data/N2002873_ZJX_80-451570558/combined/trimmed_data/1P_kmer_count"
merged_output_dir="/home/data/qs/Ndufs4_data/N2002873_ZJX_80-451570558/combined/trimmed_data/1P_kmer_count_merg"
result_output_dir="/home/data/qs/Ndufs4_data/N2002873_ZJX_80-451570558/combined/trimmed_data/1P_kmer_count_merg_result"

# Ensure output directories exist
mkdir -p "$output_dir"ll

mkdir -p "$merged_output_dir"
mkdir -p "$result_output_dir"

# Initialize a count to keep track of concurrently running jobs
max_jobs=2
count=0

# Read the sample list and process each sample
while IFS= read -r sample; do
  echo "Processing sample: $sample"

  # Function to process each sample
  process_sample() {
    # Step 1: Counting isoform specific kmer
    python "$kmer_count_script" --fastq_path "$input_fastq_dir/${sample}_trim_1P.fastq" --num_threads 20 --chunk_size 10000 --csv_input_dir "$reference_dir" --csv_output_dir "$output_dir/output_csv_${sample}"

    # Step 2: Merging the isoform specific reference csv and corresponding counting csv file
    python "$merge_script" --kmer_reference_directory "$reference_dir" --kmer_counts_directory "$output_dir/output_csv_${sample}/" --output_directory "$merged_output_dir/${sample}_merged_data" --read_length 150 --k 50

    # Step 3: GC based GaussF for isoform quantification
    python "$gaussf_script" --threshold 20 --input "$merged_output_dir/${sample}_merged_data" --output "$result_output_dir/results_file_GC_${sample}.csv"

    # Clean up intermediate directories
    #echo "Cleaning up intermediate directories for sample: $sample"
    #rm -rf "$output_dir/output_csv_${sample}"
    #rm -rf "$merged_output_dir/${sample}_merged_data"
  }

  # Run processing in background
  process_sample "$sample" &
  let count+=1

  # Wait if we have reached maximum jobs
  if [[ "$count" -ge "$max_jobs" ]]; then
    wait
    count=0
  fi

done < "$sample_list"

wait # Wait for all jobs to finish if any remains

echo "All samples processed and intermediate files cleaned up."
