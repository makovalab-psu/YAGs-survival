#!/usr/bin/env python3
"""
01_convert_trf_to_bed.py
========================
Convert Tandem Repeat Finder (TRF) .dat output to BED format.

Usage:
    python3 01_convert_trf_to_bed.py <input.dat> <output.bed>

TRF .dat column layout (data lines only):
    1   Start position (1-based)
    2   End position
    3   Period size
    4   Copy number
    5   Consensus size
    6   Percent matches
    7   Percent indels
    8   Score
    9   A content (%)
   10   C content (%)
   11   G content (%)
   12   T content (%)
   13   Entropy
   14   Consensus sequence

Output BED fields: chrom, start (0-based), end, name, score
    name = "Period_{period}_Copies_{copies}"
"""

import sys


def convert_trf_to_bed(input_file: str, output_file: str) -> None:
    chromosome = None
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            # Header lines identify the sequence being described
            if line.startswith("Sequence:"):
                chromosome = line.split()[1]
                continue
            if line.startswith("Parameters:") or line.strip() == "" or line.startswith("@"):
                continue

            fields = line.strip().split()
            if len(fields) < 14:
                continue

            start  = int(fields[0]) - 1   # convert 1-based TRF → 0-based BED
            end    = int(fields[1])
            period = fields[2]
            copies = fields[3]
            score  = fields[7]

            name = f"Period_{period}_Copies_{copies}"
            fout.write(f"{chromosome}\t{start}\t{end}\t{name}\t{score}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.dat> <output.bed>")
        sys.exit(1)
    convert_trf_to_bed(sys.argv[1], sys.argv[2])
