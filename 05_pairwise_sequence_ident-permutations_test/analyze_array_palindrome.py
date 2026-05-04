#!/usr/bin/env python3
"""
Analyze distance vs similarity for YAG gene copies in gene arrays and palindromes.
Creates separate plots for:
1. Gene copies in gene arrays (from S14 table)
2. Gene copies in palindrome regions (from S13 table) with Q-number based coloring
"""

import os
import sys
import gzip
import subprocess
import tempfile
import re
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from urllib.parse import quote

import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Constants
SHEET_ID = "1Fak_7kYxBVyFmQ0HPcEGn9U2z4kbqwGtZSLt2CSHw3k"
GID_S1 = "0"
GID_S13_PALINDROMES = "766397806"
GID_S14_ARRAY = "987600569"
GFF_DOWNLOAD_DIR = "gff_files"
OUTPUT_DIR = "./output"
SEQUENCES_DIR = "./sequences"
ALIGNMENTS_DIR = "./alignments"


## Fore each reference genome the Y chromosome contig was extracted with samtools faidx
## e.g.: samtools faidx ./GCF_028878055.2_NHGRI_mSymSyn1-v2.0_pri_genomic.fna NC_072448.2 > NC_072448.2.fa
REFERENCE_GENOMES = {
    "GorGor": "/Users/kxp5629/proj/Y/data/references/GorGor/ncbi_dataset_v2.0/data/GCF_029281585.2/NC_073248.2.fa",
    "PonAbe": "/Users/kxp5629/proj/Y/data/references/PonAbe/ncbi_dataset_v2.0/data/GCF_028885655.2/NC_072009.2.fa",
    "PonPyg": "/Users/kxp5629/proj/Y/data/references/PonPyg/ncbi_dataset_v2.0/data/GCF_028885625.2/NC_072397.2.fa",
    "PanPan": "/Users/kxp5629/proj/Y/data/references/PanPan/ncbi_dataset_v2.0/data/GCF_029289425.2/NC_073273.2.fa",
    "PanTro": "/Users/kxp5629/proj/Y/data/references/PanTro/ncbi_dataset_v2.0/data/GCF_028858775.2/NC_072422.2.fa",
    "SymSyn": "/Users/kxp5629/proj/Y/data/references/SymSyn/ncbi_dataset_v2.0/data/GCF_028878055.2/NC_072448.2.fa",
    "HomSap": "/Users/kxp5629/proj/Y/data/references/HomSap/ncbi_dataset/data/GCF_009914755.1/NC_060948.1.fa",
}

CHROM_MAPPING = {
    "NC_073248.2": "chrY",
    "NC_060948.1": "chrY",
    "NC_072422.2": "chrY",
    "NC_073273.2": "chrY",
    "NC_072009.2": "chrY",
    "NC_072397.2": "chrY",
    "NC_072448.2": "chrY",
}


def fetch_google_sheet_by_gid(gid: str, sheet_name: str = "", skip: int = 1) -> pd.DataFrame:
    """Fetch data from Google Sheet by gid."""
    print(f"Fetching Google Sheet: {sheet_name} (gid={gid})")
    url = f"https://docs.google.com/spreadsheets/d/{SHEET_ID}/export?format=csv&gid={gid}"
    df = pd.read_csv(url, skiprows=skip)
    # print(df)
    print(f"Found {len(df)} entries")
    return df


def parse_gff_for_gene(gff_path: str, gene_id: str) -> Optional[Dict]:
    """Parse GFF file to extract all CDS coordinates for a specific gene."""
    if not os.path.exists(gff_path):
        return None

    direct_exons = []
    transcript_ids: set[str] = set()
    chrom = None
    strand = None

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            feature_type = parts[2]
            attrs = parts[8] if len(parts) >= 9 else parts[7]

            if feature_type == 'CDS':
                if (f'gene={gene_id};' in attrs or f'gene={gene_id},' in attrs
                        or attrs.endswith(f'gene={gene_id}')
                        or f'ID=CDS-{gene_id}' in attrs or f'ID=CDS_{gene_id}' in attrs
                        or f'ID=cds-{gene_id}' in attrs):
                    chrom = parts[0]
                    strand = parts[6]
                    direct_exons.append({'start': int(parts[3]), 'end': int(parts[4])})

            elif feature_type in ('mRNA', 'transcript'):
                if f'Parent={gene_id}' in attrs or f'Parent={gene_id};' in attrs:
                    for attr in attrs.split(';'):
                        if attr.startswith('ID='):
                            transcript_ids.add(attr[3:])
                            break

    if direct_exons:
        direct_exons.sort(key=lambda x: x['start'])
        return {'chrom': chrom, 'exons': direct_exons, 'strand': strand, 'gene_id': gene_id}

    # Fallback: follow Parent chain for Liftoff-style annotations
    if not transcript_ids:
        return None

    fallback_exons = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attrs = parts[8]
            for tid in transcript_ids:
                if f'Parent={tid}' in attrs or f'Parent={tid};' in attrs:
                    chrom = parts[0]
                    strand = parts[6]
                    fallback_exons.append({'start': int(parts[3]), 'end': int(parts[4])})
                    break

    if not fallback_exons:
        return None

    fallback_exons.sort(key=lambda x: x['start'])
    return {'chrom': chrom, 'exons': fallback_exons, 'strand': strand, 'gene_id': gene_id}


def extract_sequence(genome_path: str, chrom: str, start: int, end: int, strand: str) -> Optional[str]:
    """Extract sequence from reference genome."""
    genome_path = genome_path.strip()  # Remove trailing whitespace
    if not os.path.exists(genome_path):
        return None

    if genome_path.endswith('.gz'):
        handle = gzip.open(genome_path, 'rt')
    else:
        handle = open(genome_path, 'r')

    try:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == chrom or chrom in record.description:
                seq = str(record.seq[start-1:end])
                if strand == '-':
                    seq = str(Seq(seq).reverse_complement())
                return seq
    finally:
        handle.close()
    return None


def extract_gene_sequence(genome_path: str, gene_info: Dict) -> Optional[str]:
    """Extract multi-exon gene sequence from reference genome."""
    genome_path = genome_path.strip()  # Remove trailing whitespace
    if not os.path.exists(genome_path):
        return None

    chrom = gene_info['chrom']
    exons = gene_info['exons']
    strand = gene_info['strand']

    if genome_path.endswith('.gz'):
        handle = gzip.open(genome_path, 'rt')
    else:
        handle = open(genome_path, 'r')

    try:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == chrom or chrom in record.description:
                # Extract all exon sequences
                exon_seqs = []
                for exon in exons:
                    seq = str(record.seq[exon['start']-1:exon['end']])
                    exon_seqs.append(seq)

                # Concatenate exons
                full_seq = ''.join(exon_seqs)

                # Reverse complement if on minus strand
                if strand == '-':
                    full_seq = str(Seq(full_seq).reverse_complement())

                return full_seq
    finally:
        handle.close()
    return None


def save_sequence(gene_id: str, sequence: str, species: str, gene_family: str, analysis_type: str):
    """Save sequence to file for review and reuse."""
    os.makedirs(SEQUENCES_DIR, exist_ok=True)
    seq_dir = os.path.join(SEQUENCES_DIR, f"{species}_{gene_family}_{analysis_type}")
    os.makedirs(seq_dir, exist_ok=True)

    seq_file = os.path.join(seq_dir, f"{gene_id}.fasta")
    with open(seq_file, 'w') as f:
        f.write(f">{gene_id}\n{sequence}\n")
    return seq_file


def load_sequence(gene_id: str, species: str, gene_family: str, analysis_type: str) -> Optional[str]:
    """Load previously saved sequence."""
    seq_file = os.path.join(SEQUENCES_DIR, f"{species}_{gene_family}_{analysis_type}", f"{gene_id}.fasta")
    if os.path.exists(seq_file):
        with open(seq_file, 'r') as f:
            lines = f.readlines()
            if len(lines) >= 2:
                return ''.join(line.strip() for line in lines[1:])
    return None


def get_alignment_filename(gene1: str, gene2: str, species: str, gene_family: str, analysis_type: str) -> str:
    """Generate alignment filename."""
    os.makedirs(ALIGNMENTS_DIR, exist_ok=True)
    aln_dir = os.path.join(ALIGNMENTS_DIR, f"{species}_{gene_family}_{analysis_type}")
    os.makedirs(aln_dir, exist_ok=True)

    genes_sorted = tuple(sorted([gene1, gene2]))
    aln_file = os.path.join(aln_dir, f"{genes_sorted[0]}__vs__{genes_sorted[1]}.fasta")
    return aln_file


def align_sequences_mafft(sequences: List[Tuple[str, str]],
                          species: str = "", gene_family: str = "",
                          analysis_type: str = "") -> Optional[AlignIO.MultipleSeqAlignment]:
    """Align sequences using MAFFT, with caching support."""
    if len(sequences) < 2:
        return None

    if species and gene_family and analysis_type:
        aln_file = get_alignment_filename(sequences[0][0], sequences[1][0],
                                         species, gene_family, analysis_type)

        if os.path.exists(aln_file):
            try:
                alignment = AlignIO.read(aln_file, 'fasta')
                return alignment
            except Exception:
                pass

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        temp_input = f.name
        for name, seq in sequences:
            f.write(f">{name}\n{seq}\n")

    try:
        temp_output = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False).name
        cmd = ['mafft', '--globalpair','--maxiterate', '1000', '--quiet', temp_input]

        with open(temp_output, 'w') as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            raise Exception(f"Alignment error:\n{result}")

        alignment = AlignIO.read(temp_output, 'fasta')

        if species and gene_family and analysis_type:
            aln_file = get_alignment_filename(sequences[0][0], sequences[1][0],
                                             species, gene_family, analysis_type)
            AlignIO.write(alignment, aln_file, 'fasta')

        return alignment

    finally:
        os.unlink(temp_input)
        if os.path.exists(temp_output):
            os.unlink(temp_output)


def calculate_identity(seq1: str, seq2: str) -> dict:
    """Calculate percent identity between two aligned sequences."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be aligned")

    matches        = 0
    aligned_cols   = 0  # neither is a gap
    gap_in_seq1    = 0  # insertion in seq2 relative to seq1
    gap_in_seq2    = 0  # insertion in seq1 relative to seq2
    double_gap     = 0  # shouldn't happen in pairwise but just in case
    diff_col_indices: List[int] = []  # ungapped column index of each mismatch

    for a, b in zip(seq1, seq2):
        if a == "-" and b == "-":
            double_gap += 1 # should not be
        elif a == '-':
            gap_in_seq1 += 1
        elif b == '-':
            gap_in_seq2 += 1
        else:
            aligned_cols += 1
            if a == b:
                matches += 1
            else:
                diff_col_indices.append(aligned_cols)

    diff_positions = [round(i / aligned_cols * 100, 1) for i in diff_col_indices] if aligned_cols else []

    aln_len   = len(seq1) - double_gap
    len_seq1  = aligned_cols + gap_in_seq2   # original seq1 length
    len_seq2  = aligned_cols + gap_in_seq1   # original seq2 length
    len_short = min(len_seq1, len_seq2)
    len_long  = max(len_seq1, len_seq2)

    if aligned_cols == 0:
        return 0.0

    # return (matches / aligned_cols) * 100
    return {
        'matches':         matches,
        'len_seq1':        len_seq1,
        'len_seq2':        len_seq2,
        'pct_id_ungapped': matches/aligned_cols * 100 if aligned_cols else 0, #BLAST style
        'pct_id_aln':      matches/aln_len      * 100 if aln_len      else 0, #including gaps
        'rel_len_diff':    abs(len_seq1 - len_seq2) / len_long,
        'abs_len_diff':    abs(len_seq1 - len_seq2),
        'diff_positions':  diff_positions,

    }

def extract_q_number(palindrome_name: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract Q number and arm from palindrome name.

    Example: NC_073248.2.Q1.A -> ('Q1', 'A')
             NC_073248.2.Q5.1.B -> ('Q5.1', 'B')
    """
    match = re.search(r'\.Q([\d\.]+)\.([AB])$', palindrome_name)
    if match:
        return f"Q{match.group(1)}", match.group(2)
    return None, None


def find_palindrome_for_position(gene_start: int, gene_end: int, chrom: str, palindromes_df: pd.DataFrame) -> Optional[Dict]:
    """Find palindrome region overlapping with gene coordinates.

    Returns information about complete or partial overlap.
    """
    for _, pal in palindromes_df.iterrows():
        if pal['Chromosome'] == chrom or CHROM_MAPPING.get(pal['Chromosome']) == chrom:
            pal_start = pal['Start']
            pal_end = pal['End']

            # Gene must be fully contained in palindrome
            fully_contained = (pal_start <= gene_start and gene_end <= pal_end)
            if fully_contained:
                partial_overlap = False

                q_num, arm = extract_q_number(pal['Palindrome name'])
                return {
                    'palindrome_name': pal['Palindrome name'],
                    'q_number': q_num,
                    'arm': arm,
                    'start': pal_start,
                    'end': pal_end,
                    'fully_contained': fully_contained,
                    'partial_overlap': partial_overlap,
                    'gene_start': gene_start,
                    'gene_end': gene_end
                }


def find_array_for_position(pos: int, chrom: str, gene_family: str, arrays_df: pd.DataFrame) -> Optional[Dict]:
    """Find array region containing the given position and matching gene family."""
    for _, arr in arrays_df.iterrows():
        if arr['chrom'] == chrom or CHROM_MAPPING.get(arr['chrom']) == chrom:
            if arr['gene_families'] == gene_family:
                if arr['start'] <= pos <= arr['end']:
                    return {
                        'array_name': arr['gene_families'],
                        'start': arr['start'],
                        'end': arr['end']
                    }
    return None


def process_species_gene_combined(all_genes_df: pd.DataFrame, arrays_df: pd.DataFrame,
                                   palindromes_df: pd.DataFrame, species: str, gene: str) -> Optional[Dict]:
    """Process all gene copies for a species-gene combination (arrays + palindromes + others)."""
    print(f"\nProcessing COMBINED: {species} - {gene}")

    mask = (all_genes_df['Species'] == species) & (all_genes_df['Gene Family'] == gene) & \
           (all_genes_df['Chromosome'] == 'chrY')
    gene_data = all_genes_df[mask].copy()

    if len(gene_data) < 2:
        print(f"  Skipping: only {len(gene_data)} copies found")
        return None

    species_arrays = arrays_df[arrays_df['Species'] == species]
    species_palindromes = palindromes_df[palindromes_df['Species'] == species]

    gene_copies = []
    for _, row in gene_data.iterrows():
        gene_id = row['Gene ID']
        gff_filename = row['File']

        if pd.isna(gff_filename) or gff_filename == "":
            continue

        gff_path = os.path.join(GFF_DOWNLOAD_DIR, gff_filename)
        gene_info = parse_gff_for_gene(gff_path, gene_id)

        if gene_info is None:
            print(f"missing {gene_id}")
            continue

        gene_start = min(exon['start'] for exon in gene_info['exons'])
        gene_end = max(exon['end'] for exon in gene_info['exons'])

        # Check if in array
        array_info = find_array_for_position(gene_start, gene_info['chrom'], gene, species_arrays)

        # Check if in palindrome
        palindrome_info = find_palindrome_for_position(gene_start, gene_end, gene_info['chrom'], species_palindromes)

        # Load or extract sequence
        sequence = load_sequence(gene_id, species, gene, 'combined')
        if sequence is None:
            genome_path = REFERENCE_GENOMES.get(species)
            if genome_path is None:
                continue

            sequence = extract_gene_sequence(genome_path, gene_info)
            if sequence is None:
                print('no sequence')
                continue

            save_sequence(gene_id, sequence, species, gene, 'combined')

        gene_copies.append({
            'gene_id': gene_id,
            'start': gene_start,
            'end': gene_end,
            'strand': gene_info['strand'],
            'sequence': sequence,
            'in_array': array_info is not None,
            'in_palindrome': palindrome_info is not None,
            'array_start': array_info['start'] if array_info else None,
            'array_end': array_info['end'] if array_info else None,
            'q_number': palindrome_info['q_number'] if palindrome_info else None,
            'arm': palindrome_info['arm'] if palindrome_info else None,
            'palindrome_fully_contained': palindrome_info['fully_contained'] if palindrome_info else None,
            'palindrome_partial_overlap': palindrome_info['partial_overlap'] if palindrome_info else None
        })

    if len(gene_copies) < 2:
        print(f"  Skipping: only {len(gene_copies)} valid copies")
        return None

    print(f"  Successfully extracted {len(gene_copies)} gene copies")

    comparisons = []
    for i in range(len(gene_copies)):
        for j in range(i+1, len(gene_copies)):
            copy1 = gene_copies[i]
            copy2 = gene_copies[j]

            palindrome_tag = ""

            # If both genes are in palindromes, only compare if same palindrome and opposite arms
            if copy1['in_palindrome'] and copy2['in_palindrome']:
                palindrome_tag = "palindrome"
                main_q1 = copy1['q_number'].split('.')[0] if copy1['q_number'] else None
                main_q2 = copy2['q_number'].split('.')[0] if copy2['q_number'] else None
                if copy1['q_number'] == copy2['q_number']:
                    palindrome_tag += "_sameQ"
                    if copy1['arm'] != copy2['arm']: 
                        if copy1['strand'] != copy2['strand']:
                            palindrome_tag += "_trans"
                    else:
                        # if copy1['strand'] == copy2['strand']:
                        palindrome_tag += "_cis"
                elif main_q1 != main_q2:
                    palindrome_tag += "_diffQ"
                # else: same main Q but different sub-Q — leave as plain "palindrome" (excluded from stats)


            # array_tag = ""
            # If both genes are in arrays, only compare if same array region
            # this comparison does not make sense
            # if copy1['in_array'] and copy2['in_array']:
            #     if copy1['array_start'] != copy2['array_start'] or copy1['array_end'] != copy2['array_end']:
            #         array_tag = "same_array"
            #     else:
            #         array_tag = "different_array"

            distance = abs(copy1['start'] - copy2['start'])

            # Try to load existing alignment from any analysis type
            alignment = None
            for analysis_type in ['combined', 'array', 'palindrome']:
                alignment = align_sequences_mafft([
                    (copy1['gene_id'], copy1['sequence']),
                    (copy2['gene_id'], copy2['sequence'])
                ], species, gene, analysis_type)
                if alignment is not None:
                    break

            if alignment is None:
                continue

            seq1_aligned = str(alignment[0].seq)
            seq2_aligned = str(alignment[1].seq)
            res = calculate_identity(seq1_aligned, seq2_aligned)
            identity = res['pct_id_ungapped']

            # Skip if either gene is in both array and palindrome (ambiguous category)
            if (copy1['in_array'] and copy1['in_palindrome']) or \
               (copy2['in_array'] and copy2['in_palindrome']):
                continue

            # Skip mixed pairs (one in array/palindrome, the other not)
            if copy1['in_array'] != copy2['in_array']:
                continue
            if copy1['in_palindrome'] != copy2['in_palindrome']:
                continue

            # # Skip cross-array comparisons (genes from different array regions)
            # if copy1['in_array'] and copy2['in_array']:
            #     if copy1['array_start'] != copy2['array_start'] or copy1['array_end'] != copy2['array_end']:
            #         continue

            # Determine plot category
            in_array = copy1['in_array'] and copy2['in_array']
            in_palindrome = copy1['in_palindrome'] and copy2['in_palindrome']

            same_orientation = (copy1['strand'] == copy2['strand'])
            same_q_diff_arm = False
            if in_palindrome and not same_orientation:
                same_q_diff_arm = (copy1['q_number'] == copy2['q_number'] and copy1['arm'] != copy2['arm'])

            palindrome_fully_contained = copy1['palindrome_fully_contained']  and copy2['palindrome_fully_contained']
            # print(f"{copy1['palindrome_fully_contained']} {copy2['palindrome_fully_contained']} {palindrome_fully_contained}")

            same_array_region = (
                in_array and
                copy1['array_start'] == copy2['array_start'] and
                copy1['array_end'] == copy2['array_end']
            )

            comparisons.append({
                'species': species,
                'gene_family': gene,
                'gene1': copy1['gene_id'],
                'gene1_len': res['len_seq1'],
                'gene2': copy2['gene_id'],
                'gene2_len': res['len_seq2'],
                'distance': distance,
                'identity': identity,
                'relative_len_diff': res['rel_len_diff'],
                'absolute_len_diff': res['abs_len_diff'],
                'in_array': in_array,
                'in_palindrome': in_palindrome,
                'same_orientation': same_orientation,
                'same_q_diff_arm': same_q_diff_arm,
                'palindrome_tag': palindrome_tag,
                'fully_contained_both': palindrome_fully_contained,
                'gene1_palindrome_fully_contained': copy1['palindrome_fully_contained'],
                'gene1_palindrome_partial_overlap': copy1['palindrome_partial_overlap'],
                'gene2_palindrome_fully_contained': copy2['palindrome_fully_contained'],
                'gene2_palindrome_partial_overlap': copy2['palindrome_partial_overlap'],
                'same_array_region': same_array_region,
                'diff_positions': res['diff_positions'],
            })

    print(f"  Completed {len(comparisons)} pairwise comparisons")

    return {
        'species': species,
        'gene': gene,
        'comparisons': comparisons,
        'plot_type': 'combined'
    }

def millions(x, pos):
    return f'{x*1e-6:.0f}'

def hundred_K(x, pos):
    return f'{x*1e-6:.2f}'

# Apply the formatter
def plot_combined_results(result: Dict):
    """Create combined scatterplot with arrays and palindromes."""
    species = result['species']
    gene = result['gene']
    comparisons = result['comparisons']

    if len(comparisons) == 0:
        return

    # Separate by category
    arrays = [c for c in comparisons if c['in_array']]
    pal_cis = [c for c in comparisons if c['in_palindrome'] and (c['palindrome_tag'] == "palindrome_sameQ_cis")]
    pal_trans = [c for c in comparisons if c['in_palindrome'] and (c['palindrome_tag'] == "palindrome_sameQ_trans")]
    pal_other = [c for c in comparisons if c['in_palindrome'] and (c['palindrome_tag'] == 'palindrome')]
    neither = [c for c in comparisons if not c['in_array'] and not c['in_palindrome']]

    # Function to plot data on an axis
    def plot_data(ax):
        # Plot neither (gray circles)
        # Plot palindromes - orange triangles for opposite arms same Q
        if pal_trans:
            x = [c['distance'] for c in pal_trans]
            y = [c['identity'] for c in pal_trans]
            ax.scatter(x, y, marker='^', c='orange', label='Palindrome (opposite arms, same Q)',
                      s=50, alpha=0.5, edgecolors='none')
        if pal_cis:
            x = [c['distance'] for c in pal_cis]
            y = [c['identity'] for c in pal_cis]
            ax.scatter(x, y, marker='^', c='yellowgreen', label='Palindrome (same arm, same Q)',
                      s=50, alpha=0.5, edgecolors='none')
                # Plot palindromes - gray triangles for opposite orientation
        if pal_other:
            x = [c['distance'] for c in pal_other]
            y = [c['identity'] for c in pal_other]
            ax.scatter(x, y, marker='^', c='gray', label='Palindrome other',
                      s=50, alpha=0.5, edgecolors='none')
                # Plot arrays - light blue circles
        if arrays:
            x = [c['distance'] for c in arrays]
            y = [c['identity'] for c in arrays]
            ax.scatter(x, y, marker='o', c='lightblue', label='Array region',
                      s=50, alpha=0.5, edgecolors='none')
        if neither:
            x = [c['distance'] for c in neither]
            y = [c['identity'] for c in neither]
            ax.scatter(x, y, marker='o', c='gray', label='Neither',
                      s=50, alpha=0.5, edgecolors='none')


        ax.legend()
        ax.spines[['right', 'top']].set_visible(False)
        ax.grid(True, alpha=0.3)

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # Plot 1: Full view
    fig1, ax1 = plt.subplots(figsize=(10, 6))

    plot_data(ax1)
    ax1.set_xlabel('Genomic Distance (Mb)', fontsize=12)
    ax1.set_ylabel('Sequence Identity (%)', fontsize=12)
    # ax1.set_title(f'{species} - {gene} (Combined)\nDistance vs Sequence Identity',
    #               fontsize=14, fontweight='bold')
    ax1.set_title(f"Pairwise sequence identity in relationship with genomic distance.")

    # Apply the formatter
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(millions))

    plt.grid(axis='y', color='0.9', linewidth=1)
    plt.grid(axis='x', color='0.9', linewidth=1)   

    filename1 = f"{species}_{gene}_COMBINED_distance_similarity_{timestamp}.pdf"
    filepath1 = os.path.join(OUTPUT_DIR, filename1)
    plt.tight_layout()
    plt.savefig(filepath1, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {filepath1}")

    # Plot 2: Zoomed view (95-100% identity, 0-2000000 distance)
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    plot_data(ax2)
    ax2.set_xlim(0, 2000000)
    ax2.set_ylim(95, 100)
    ax2.set_xlabel('Genomic Distance (Mb)', fontsize=12)
    ax2.set_ylabel('Sequence Identity (%)', fontsize=12)
    ax2.set_title(f'{species} - {gene} (Zoomed)\n95-100% Identity, 0-2M bp Distance',
                  fontsize=14, fontweight='bold')
    
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(hundred_K))

    filename2 = f"{species}_{gene}_COMBINED_distance_similarity_ZOOMED_{timestamp}.pdf"
    filepath2 = os.path.join(OUTPUT_DIR, filename2)
    plt.tight_layout()
    plt.savefig(filepath2, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {filepath2}")

    # Save comparison table
    csv_filename = f"{species}_{gene}_COMBINED_comparisons_{timestamp}.csv"
    csv_filepath = os.path.join(OUTPUT_DIR, csv_filename)
    df = pd.DataFrame(comparisons)



    # Check if species and gene_family are in the comparisons (for All_All combined)
    df = df[['species', 'gene_family', 'gene1', 'gene1_len', 'gene2','gene2_len', 'distance', 'identity', 'relative_len_diff', 'absolute_len_diff', 'in_array', 'in_palindrome',
             'palindrome_tag', 'fully_contained_both', 'gene1_palindrome_fully_contained', 'gene1_palindrome_partial_overlap',
             'gene2_palindrome_fully_contained', 'gene2_palindrome_partial_overlap', 'same_array_region', 'diff_positions']].rename(columns={
        'species': 'Species', 'gene_family': 'GeneFamily',
        'gene1': 'GeneA', 'gene1_len': 'GeneA_length','gene2': 'GeneB', 'gene2_len': 'GeneB_length',
        'distance': 'Distance', 'identity': 'Similarity',
        'relative_len_diff':'Relative_lenght_difference',
        'absolute_len_diff':'Absolude_length_difference',
        'palindrome_tag': "PalindromeTag",
        'fully_contained_both' : 'BothGenesFullyInPalindrome',
        'gene1_palindrome_fully_contained': 'GeneA_Fully_In_Palindrome',
        'gene1_palindrome_partial_overlap': 'GeneA_Partial_Palindrome_Overlap',
        'gene2_palindrome_fully_contained': 'GeneB_Fully_In_Palindrome',
        'gene2_palindrome_partial_overlap': 'GeneB_Partial_Palindrome_Overlap',
        'same_array_region': 'SameArrayRegion',
        'diff_positions':'Diff_positions'
    })

    # Add Type column based on in_array and in_palindrome
    df['Type'] = df.apply(lambda row: 'ARRAY' if row['in_array'] else ('PALINDROME' if row['in_palindrome'] else 'OTHER'), axis=1)
    df = df.drop(columns=['in_array', 'in_palindrome'])
    df = df.sort_values('Distance')
    df.to_csv(csv_filepath, index=False)
    print(f"  Table saved: {csv_filepath}")


def main():
    """Main execution function."""
    # Parse command-line arguments
    filter_species = None
    filter_gene = None

    if len(sys.argv) > 1:
        filter_species = sys.argv[1]
    if len(sys.argv) > 2:
        filter_gene = sys.argv[2]

    print("=" * 60)
    print("Distance vs Similarity Analysis")
    print("Gene Arrays and Palindromes")
    if filter_species or filter_gene:
        print(f"Filtering: Species={filter_species or 'ALL'}, Gene={filter_gene or 'ALL'}")
    print("=" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(GFF_DOWNLOAD_DIR, exist_ok=True)


    print("\n" + "=" * 60)
    print("COMBINED (ALL GENES)")
    print("=" * 60)

    if 'all_genes_df' not in locals():
        all_genes_df = fetch_google_sheet_by_gid(GID_S1, "S1. YAG genes accessions")
        all_genes_df = all_genes_df[all_genes_df['Use'] == 'yes']
        # assert False

    if 'arrays_df' not in locals():
        arrays_df = fetch_google_sheet_by_gid(GID_S14_ARRAY, "S14. Array", 1)
        arrays_df = arrays_df[arrays_df["species"].notna()]
        arrays_df.rename(columns={'species': 'Species'}, inplace=True)

    if 'palindromes_df' not in locals():
        palindromes_df = fetch_google_sheet_by_gid(GID_S13_PALINDROMES, "S13. Palindromes")

    species_list = all_genes_df['Species'].unique()
    gene_list = all_genes_df['Gene Family'].unique()

    # Apply filters
    if filter_species:
        species_list = [s for s in species_list if s == filter_species]
    if filter_gene:
        gene_list = [g for g in gene_list if g == filter_gene]

    print(f"\nSpecies: {list(species_list)}")
    print(f"Genes: {list(gene_list)}")

    combined_results = []
    for species in species_list:
        for gene in gene_list:
            result = process_species_gene_combined(all_genes_df, arrays_df, palindromes_df, species, gene)
            if result:
                combined_results.append(result)

    print(f"\nCombined analysis: {len(combined_results)} species-gene combinations")

    # Create one combined plot with all data
    if combined_results:
        print("\nCreating combined plot with all species-gene data...")
        all_comparisons = []
        for result in combined_results:
            all_comparisons.extend(result['comparisons'])

        combined_all = {
            'species': 'All',
            'gene': 'All',
            'comparisons': all_comparisons,
            'plot_type': 'combined'
        }
        # print(combined_all)
        # assert False
        plot_combined_results(combined_all)

    print(f"\n{'=' * 60}")
    print(f"Analysis complete!")
    print(f"Plots: {len(combined_results)}")
    print(f"Plots saved to: {OUTPUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
