#!/bin/bash

# Assign arguments to variables
id_file=$1
fasta_file=$2
output_file=$3

# Loop through the IDs in the id_file
while IFS= read -r id; do
    # Use awk to extract the sequence with the given ID from the fasta file
    awk -v id="$id" 'BEGIN {RS=">"; ORS=""} NR > 1 && $1 == id {print ">" $0}' "$fasta_file" >> "$output_file"
done < "$id_file"


