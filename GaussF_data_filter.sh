#!/bin/bash

# Define the input and output files
transcript_file="/home/data/qs/HGC20231021002-0004/trimmed_data/protein_folding_isoform_list.txt"
csv_file="combined_results_GaussF_GC_mus_RPM.csv"
output_file="/home/data/qs/HGC20231021002-0004/trimmed_data/protein_folding_isoform_list_GaussF_RPKM.txt"

# Empty the output file to start fresh
> "$output_file"

# Loop through each transcript ID in the list
while IFS= read -r transcript_id; do
    # Use awk to search for the transcript ID in the third column
    awk -F, -v id="$transcript_id" '$3 ~ id' "$csv_file" >> "$output_file"
done < "$transcript_file"

echo "Search results saved to $output_file"

