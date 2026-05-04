#!/usr/bin/env python3

import re
from Bio import SeqIO
import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np


def plot_length_distribution(lengths, output_prefix, motif):

    if not lengths:
        print("no sequences found")
        return

    lengths_array = np.array(lengths)
    mean_len = np.mean(lengths_array)
    median_len = np.median(lengths_array)
    std_len = np.std(lengths_array)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Histogram
    ax1.hist(lengths, bins=min(50, len(set(lengths))), alpha=0.7, edgecolor='black')
    ax1.axvline(mean_len, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_len:.1f}')
    ax1.axvline(median_len, color='orange', linestyle='--', linewidth=2, label=f'Median: {median_len:.1f}')
    ax1.set_xlabel('Sequence Length (nucleotides)')
    ax1.set_ylabel('Frequency')
    ax1.set_yscale('log')
    ax1.set_title(f'Distribution of Extracted Sequence Lengths\nMotif: {motif}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Box plot
    ax2.boxplot(lengths, vert=True)
    ax2.set_ylabel('Sequence Length (nucleotides)')
    ax2.set_yscale('log')
    ax2.set_title(f'Box Plot of Sequence Lengths\n(n={len(lengths)} sequences)')
    ax2.grid(True, alpha=0.3)
    
    # Add statistics text
    stats_text = f'Statistics:\nMean: {mean_len:.1f}\nMedian: {median_len:.1f}\nStd: {std_len:.1f}\nMin: {min(lengths)}\nMax: {max(lengths)}'
    ax2.text(1.1, max(lengths), stats_text, transform=ax2.transData, 
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8),
             verticalalignment='top')
    
    plt.tight_layout()
    
    # Save the plot
    plot_filename = f"{output_prefix}_length_distribution.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Length distribution plot saved as: {plot_filename}")
    
    # Also save as PDF for publications
    pdf_filename = f"{output_prefix}_length_distribution.pdf"
    plt.savefig(pdf_filename, bbox_inches='tight')
    print(f"Length distribution plot saved as: {pdf_filename}")
    
    plt.show()
    
    # Print summary statistics
    print(f"\nSequence Length Statistics:")
    print(f"Total sequences: {len(lengths)}")
    print(f"Mean length: {mean_len:.2f} nucleotides")
    print(f"Median length: {median_len:.2f} nucleotides")
    print(f"Standard deviation: {std_len:.2f}")
    print(f"Range: {min(lengths)} - {max(lengths)} nucleotides")
    
    return lengths_array


def find_rightmost_motif_and_extract(sequence, motif, case_sensitive=False):
    """Find rightmost occurence of motif in read sequence."""
     
    search_seq = sequence if case_sensitive else sequence.upper()
    search_motif = motif if case_sensitive else motif.upper()

    iupac_codes = {
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }

    regex_motif = search_motif

    for code, pattern in iupac_codes.items():
        regex_motif = regex_motif.replace(code, pattern)

    matches = list(re.finditer(regex_motif,search_seq))

    if not matches:
        return False, -1 , ""
    
    rightmost_match = matches[-1]
    motif_end_pos = rightmost_match.end()
    motif_start_pos = rightmost_match.start()

    extracted_seq = sequence[motif_start_pos:]

    return True, motif_end_pos, extracted_seq


def process_fasta_file(input_files, output_file, output_full_sequences, motif, case_sensitive=False, min_length=0, plot_distribution=True):
    """Process fasta file"""
   
    sum_length = 0
    extracted_lengths = []

    
    with open(output_file, 'w') as out_handle, open(output_full_sequences, 'w') as out_full_handle:
        for input_file in input_files:
            print(f"Processing {input_file}")
            found_count = 0
            total_count = 0
            for record in SeqIO.parse(input_file, "fastq"):
                total_count += 1
                sequence = str(record.seq)
                sum_length += len(str(record.seq))

                found, position, extracted_seq = find_rightmost_motif_and_extract(
                    sequence, motif, case_sensitive
                )

                if found and len(extracted_seq) >= min_length:
                    found_count += 1
                    extracted_lengths.append(len(extracted_seq))

                    new_description = f"{record.description} | RBD_pos:{position} | extracted_length:{len(extracted_seq)}"
                    out_handle.write(f">{new_description}\n")

                    for i in range(0, len(extracted_seq), 80):
                        out_handle.write(f"{extracted_seq[i:i+80]}\n")

                    out_full_handle.write(f">{record.id}")
                    out_full_handle.write(f"{sequence}")


            print(f"Processed {total_count} sequences")
            print(f"Found motif in {found_count} sequences\n")

    print(f"Results written to {output_file}")
    print(f"Average length = {sum_length/total_count}\n")

    if plot_distribution and extracted_lengths:
        output_prefix = output_file.rsplit('.', 1)[0]  # Remove file extension
        plot_length_distribution(extracted_lengths, output_prefix, motif)


def main():
    parser = argparse.ArgumentParser(
        description="Extract DAZ \"tails\" after righmost RBD motif occurence"
    )
    parser.add_argument("input_files", help="Input FASTQ file", nargs="*")
    parser.add_argument("-m", "--motif", default = "MHMAYCACGVYGARKCCTRTAACTCAGYAYGTTCAG", 
                        help="Motif to search for (supports IUPAC codes)")
    parser.add_argument("-o", "--output", default="extracted_sequences.fasta",
                        help="Output FASTA file (default: extracted_sequences.fasta)")
    parser.add_argument("-s", "--full_sequences", default = "extracted_sequences_full.fasta",
                        help="Output FASTA file (default: extracted_sequences_full.fasta)")
    parser.add_argument("-c", "--case-sensitive", action="store_true",
                        help="Case sensitive search")
    parser.add_argument("-l", "--min_length", type=int, default=0,
                        help="Minimum length of extracted sequences (default:0)")
    
    args = parser.parse_args()

    # PonAbe Motif for RBD = GCCTGCAGCCCCTAATCACGCTGAATCCTRTAACTCAGTACGTT

    try:
        process_fasta_file(
            args.input_files,
            args.output,
            args.full_sequences,
            args.motif, 
            args.case_sensitive,
            args.min_length
        )
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_files}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
