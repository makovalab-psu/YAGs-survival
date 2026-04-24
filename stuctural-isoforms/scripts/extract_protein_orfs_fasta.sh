#!/bin/bash

id_file="$1"
fasta_file="$2"
output_file="$3"

# Initialize or clear the output file
> "$output_file"

# Loop through each ID in the ID file
while read -r id; do
# Extract the matching sequence using awk
awk -v target_id="$id" 'BEGIN {RS=">"; ORS=""} NR > 1 && $1 == target_id {print ">" $0}' "$fasta_file" >> "$output_file"
done < "$id_file"

echo "Extraction complete. Subset saved to $output_file"
