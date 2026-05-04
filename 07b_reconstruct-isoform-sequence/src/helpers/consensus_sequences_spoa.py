#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict
from Bio.Seq import Seq
import subprocess
import re
import pandas as pd
import matplotlib.pyplot as plt
import math


def parse_arguments():
    parser = argparse.ArgumentParser(description='Create consensus sequences and find longest ORFs')
    parser.add_argument('reads_file', help='Tab-delimited file with read ID and sequence')
    parser.add_argument('mapping_file', help='Tab-delimited file with isoform, signature, and read ID')
    parser.add_argument("-m", '--min-orf-length', type=int, default=30,
                       help='Minimum ORF length in amino acids (default: 30)')
    parser.add_argument("-c", "--counts-table",help="counts table for isoform/signature combination", required=True)
    parser.add_argument("-o", "--output", help="name of output file")
    parser.add_argument("-p", "--plot", help="file name of PNG plot with counts per signature per isoform", required=True)
    parser.add_argument("-s", "--species", help="name of species being processed (for catalogueing purposes)", required=True)

    return parser.parse_args()

def cut_trailing_as(sequence, min_length, threshold = 0.94):
    """ Cut trailing A's from a string sequence if there are more than min_length A' (a is threshod*100 percentage)."""
    if not sequence:
        return sequence

    for i in range(len(sequence)-min_length, -1, -1):
        window = sequence[i:i+min_length]
        a_content = window.count('A') / len(window)
        if a_content >= threshold:
            continue
        else:
            #add 10% of the sequence back in case the window is cutting beyond the poly A tail
            to =i + (len(sequence) -i)/10
            return sequence[:math.ceil(to)]

    return sequence

def read_sequences(reads_file):
    """Read sequences into dictionary with read_id as key"""
    sequences = {}
    try:
        with open(reads_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 4:
                    print(f"Warning: Line {line_num} in reads file doesn't have 4 columns: {line}", file=sys.stderr)
                    # this essentially filters out seqquences extracted from the reference
                    continue
                read_id, sequence, signature, quality = parts
                sequences[read_id] = {}
                sequences[read_id]['seq'] = sequence
                sequences[read_id]['qual'] = quality
                avg_perbase_qual = sum([ord(x) - 33 for x in quality])
                sequences[read_id]["avg_qual"] = avg_perbase_qual /len(quality)

    except FileNotFoundError:
        print(f"Error: Could not find reads file: {reads_file}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(sequences)} sequences", file=sys.stderr)
    return sequences

def read_mapping(mapping_file):
    """Read mapping file and group read IDs by (isoform, signature)"""
    groups = defaultdict(list)
    try:
        with open(mapping_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 3:
                    print(f"Warning: Line {line_num} in mapping file doesn't have 3 columns: {line}", file=sys.stderr)
                    continue
                signature, read_id, isoform = parts
                groups[(isoform, signature)].append(read_id)
    except FileNotFoundError:
        print(f"Error: Could not find mapping file: {mapping_file}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(groups)} unique isoform-signature combinations", file=sys.stderr)
    return groups

def create_consensus_msa(sequences, output, isoform, signature):
    """Create consensus using multiple sequence alignment approach"""
    if not sequences:
        return ""

    if len(sequences) == 1:
        return sequences[0]

    ordered_sequences = {k:v for k, v in sorted(sequences.items(), key=lambda item: item[1]["avg_qual"])}
    print("extracting fasta")

    try:
        # Create temporary FASTA file
        # with tempfile.NamedTemporaryFil(mode='w', suffix='.fasta', delete=False) as temp_fasta:
        out = f"{output}_{isoform}_{signature}.fasta"
        with open(out, 'w') as out_fasta:
            for key, seq in ordered_sequences.items():
                # out_fasta.write(f"@seq{i}\n{seq['seq']}\n+\n{seq['qual']}\n")
                out_fasta.write(f">{key}\n{seq['seq']}\n")
            temp_fasta_path = out_fasta.name

        # Use SPOA to extract consensus
        try:
            # Try to run muscle
            # with tempfile.NamedTemporaryFile(mode='w', suffix='.aln', delete=False) as temp_aln:
            print("running subprocess")
            result = subprocess.run(['spoa', temp_fasta_path],
                                  capture_output=True, text=True, timeout=300)
            consensus = result.stdout.split("\n")[1] # remove fasta header >Consensus...

        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
            print("Muscle not available or failed", file=sys.stderr)
            sys.exit(1)

        # os.unlink(temp_fasta_path)
        return consensus

    except Exception as e:
        print(f"Warning: MSA failed ({e})", file=sys.stderr)
        sys.exit(1)


def find_orfs_biopython(sequence, min_length=30):
    """Find ORFs using BioPython"""
    seq_obj = Seq(sequence)
    orfs = []

    # Check all 6 reading frames
    for strand, seq in [(1, seq_obj), (-1, seq_obj.reverse_complement())]:
        for frame in range(3):
            # Translate in this frame
            subseq = seq[frame:]
            if len(subseq) < 3:
                continue

            # Translate to protein
            protein = subseq.translate()
            protein_str = str(protein)

            # Find ORFs (start with M, end with *)
            start = 0
            while True:
                # Find next start codon
                start_pos = protein_str.find('M', start)
                if start_pos == -1:
                    break

                # Find next stop codon
                stop_pos = protein_str.find('*', start_pos)
                if stop_pos == -1:
                    # No stop codon found, ORF goes to end
                    orf_protein = protein_str[start_pos:]
                    orf_length = len(orf_protein)
                else:
                    # Stop codon found
                    orf_protein = protein_str[start_pos:stop_pos]
                    orf_length = len(orf_protein)

                # Check if ORF meets minimum length
                if orf_length >= min_length:
                    # Calculate nucleotide positions
                    if strand == 1:
                        nt_start = frame + start_pos * 3
                        nt_end = frame + (start_pos + orf_length) * 3
                        frame_name = f"+{frame + 1}"
                    else:
                        nt_start = len(sequence) - (frame + (start_pos + orf_length) * 3)
                        nt_end = len(sequence) - (frame + start_pos * 3)
                        frame_name = f"-{frame + 1}"

                    # Get the DNA sequence
                    if strand == 1:
                        orf_dna = sequence[nt_start:nt_end]
                    else:
                        orf_dna = str(seq_obj[nt_start:nt_end].reverse_complement())

                    orfs.append({
                        'start': nt_start,
                        'end': nt_end,
                        'length': orf_length,
                        'protein': orf_protein,
                        'dna': orf_dna,
                        'frame': frame_name,
                        'strand': strand
                    })

                # Move to next position
                if stop_pos == -1:
                    break
                start = stop_pos + 1

    return orfs


def parse_cell_format(cell_value):
    pattern = r'^(\d+(?:\.\d+)?)\((\d+(?:\.\d+)?),(\d+(?:\.\d+)?),(\d+(?:\.\d+)?)\)$'
    match = re.match(pattern, str(cell_value).strip())

    if not match:
        pattern = r'^(\d+(?:\.\d+)?)\((\d+(?:\.\d+)?),(\d+(?:\.\d+)?),(\d+(?:\.\d+)?)\)\+(\d+(?:\.\d+)?)$'
        match = re.match(pattern, str(cell_value).strip())
        if not match:
            return None
        n, x, y, z, r = map(int, match.groups())
    else:

        n, x, y, z = map(int, match.groups())
    at_least_two_positive = sum(1 for val in [x, y, z] if val > 0) >= 2

    return {'n': n, 'x': x, 'y': y, 'z': z, 'at_least_two_positive': at_least_two_positive}


def main():
    args = parse_arguments()

    # Read input files
    sequences = read_sequences(args.reads_file)
    groups = read_mapping(args.mapping_file)
    print(groups)

    
    number_of_used_reads = sum([len(value) for key,value in groups.items()])
    limit = math.ceil(number_of_used_reads/100*2) # to process a group for consensus it has to have at least "limit" reads

    try:
        table = pd.read_csv(args.counts_table, sep='\t', index_col=0)

        if 'Column_Total' in table.index:
            table = table.drop('Column_Total')

        table_numeric = table.copy()
        plt.figure(figsize=(10, 6))


        for column in table.columns:
            if column != 'Row_Total':
                table_numeric[column] = table[column].astype(str).str.extract(r'^(\d+)')[0].fillna(0).astype(int)

                plt.plot(table.index, table_numeric[column], label=column, linewidth=2)


        plt.axhline(y=limit, color='orange', linestyle='-.', alpha = 0.7, linewidth=1.5)
        plt.axhline(y=5, color='red', linestyle='--', alpha=0.7, linewidth=1.5)

        plt.ylim(bottom=0)

        y_min, y_max = plt.ylim()
        actual_min = min(y_min, 0)  # Make sure we start from 0 or the actual minimum
        plt.axhspan(0, 5, alpha=0.1, color='pink', zorder=0)

        plt.xlabel('Signatures')
        plt.ylabel('Read Counts')
        plt.title('Isoform Expression Across Signatures')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  # Place legend outside plot
        plt.grid(True, alpha=0.3)

        plt.xticks(rotation=90, ha='right')
        plt.tight_layout()

        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        plt.savefig(f"{args.plot}")

    except Exception as e:
        print(f"Error: failed to load counts table, {e}", file = sys.stderr)
        sys.exit(1)
    # Process each group

    print(table)

    print("\nProcessing groups...")
    with open(args.output, "w") as outfile:
        for (isoform, signature), read_ids in groups.items():
            print(f"\nIsoform: {isoform}, Signature: {signature}")
            try:
                cell = table.loc[signature][isoform]
                cell_value = parse_cell_format(cell)
            except Exception as e:
                print(f"An error occured while extracting counts from the counts table for sig:{signature}, iso:{isoform}")
                print(e)
                sys.exit(1)

            try:
                if limit < 5:
                    limit = 5
                if (cell_value['n'] < limit or not cell_value['at_least_two_positive']):
                    continue
            except Exception as e:
                print(f"Issue with cell_value: {cell_value}", file = sys.stderr)
                print(f"Original cell value: {cell}")
                print(f"Error: {e}")
                sys.exit(1)
            # Collect sequences for this group
            group_sequences = {}
            missing_reads = []

            for read_id in read_ids:
                if read_id in sequences:
                    group_sequences[read_id] = sequences[read_id]
                else:
                    missing_reads.append(read_id)

            if missing_reads:
                print(f"  Warning: {len(missing_reads)} read IDs not found in sequences file")

            if not group_sequences:
                print(f"  No sequences found for this group")
                continue

            print(f"  Found {len(group_sequences)} sequences")

            # Create consensus
            consensus = create_consensus_msa(group_sequences, args.output, isoform, signature)
            print(f"  Consensus length: {len(consensus)} bp")

            if not consensus:
                print(f"  Empty consensus sequence")
                continue
            consensus = cut_trailing_as(consensus, 50)
            # Find ORFs using BioPython
            orfs = find_orfs_biopython(consensus, args.min_orf_length)

            if orfs:
                # Sort by length (longest first)
                orfs.sort(key=lambda x: x['length'], reverse=True)
                longest_orf = orfs[0]

                # print(f"  Found {len(orfs)} ORFs (≥{args.min_orf_length} aa)", file = outfile)
                # print(f"  Longest ORF:", file = outfile)
                # print(f"    Frame: {longest_orf['frame']}", file = outfile)
                # print(f"    Position: {longest_orf['start']}-{longest_orf['end']}", file = outfile)
                # print(f"    Length: {longest_orf['length']} amino acids", file = outfile)
                # print(f"    Protein: {longest_orf['protein'][:50]}{'...' if len(longest_orf['protein']) > 50 else ''}", file = outfile)

                print(f">{args.species}::{isoform}::{signature}::{len(group_sequences)}_reads", file = outfile)
                print(f"{longest_orf['protein']}", file = outfile)

            else:
                print(f"  No ORFs found (≥{args.min_orf_length} amino acids)")

if __name__ == "__main__":
    main()
