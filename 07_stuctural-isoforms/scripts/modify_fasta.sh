#!/bin/bash

# Input FASTA file from the argument
input_fasta="$1"
output_fasta="$2"

# Extract the base name without the .fasta extension
base_name=$(basename "$input_fasta" .fasta)

# Replace the names in the FASTA headers
awk -v suffix="_$base_name" '
/^>/ {
    # Remove the _1 suffix and append the custom suffix
    sub(/_[0-9]+$/, "", $1)
    $1 = $1 suffix
}
1' "$input_fasta" > "$output_fasta"

echo "Modified FASTA file saved as $output_fasta"

