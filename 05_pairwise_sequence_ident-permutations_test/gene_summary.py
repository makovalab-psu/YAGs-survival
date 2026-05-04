#!/usr/bin/env python3
"""
Per-gene summary table: one row per gene with palindrome and array annotations.
"""

import os
from datetime import datetime

import pandas as pd

from analyze_array_palindrome import (
    fetch_google_sheet_by_gid,
    parse_gff_for_gene,
    find_array_for_position,
    extract_q_number,
    CHROM_MAPPING,
    GID_S1, GID_S13_PALINDROMES, GID_S14_ARRAY,
    GFF_DOWNLOAD_DIR,
)


def find_palindrome_any_overlap(gene_start: int, gene_end: int, chrom: str,
                                palindromes_df: pd.DataFrame):
    for _, pal in palindromes_df.iterrows():
        if pal['Chromosome'] == chrom or CHROM_MAPPING.get(pal['Chromosome']) == chrom:
            pal_start, pal_end = pal['Start'], pal['End']
            if gene_start < pal_end and gene_end > pal_start:
                fully_contained = pal_start <= gene_start and gene_end <= pal_end
                q_num, arm = extract_q_number(pal['Palindrome name'])
                return {
                    'palindrome_name': pal['Palindrome name'],
                    'q_number': q_num,
                    'arm': arm,
                    'start': pal_start,
                    'end': pal_end,
                    'fully_contained': fully_contained,
                    'partial_overlap': not fully_contained,
                }
    return None


def build_gene_table(all_genes_df: pd.DataFrame, arrays_df: pd.DataFrame,
                     palindromes_df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for _, row in all_genes_df.iterrows():
        species = row['Species']
        gene_family = row['Gene Family']
        gene_id = row['Gene ID']
        gff_filename = row['File']

        if pd.isna(gff_filename) or gff_filename == "" or row.get('Chromosome') != 'chrY':
            continue

        gene_info = parse_gff_for_gene(os.path.join(GFF_DOWNLOAD_DIR, gff_filename), gene_id)
        if gene_info is None:
            continue

        gene_start = min(e['start'] for e in gene_info['exons'])
        gene_end = max(e['end'] for e in gene_info['exons'])

        sp_arrays = arrays_df[arrays_df['Species'] == species]
        sp_palindromes = palindromes_df[palindromes_df['Species'] == species]

        array_info = find_array_for_position(gene_start, gene_info['chrom'], gene_family, sp_arrays)
        palindrome_info = find_palindrome_any_overlap(gene_start, gene_end, gene_info['chrom'], sp_palindromes)

        if gene_family == "DAZ":
            print(f"  DAZ debug | {species} {gene_id} | chrom={gene_info['chrom']} start={gene_start} end={gene_end} | "
                  f"sp_palindromes={len(sp_palindromes)} | palindrome_info={palindrome_info is not None}")
            if len(sp_palindromes) > 0:
                print(f"    palindrome chroms: {sp_palindromes['Chromosome'].unique()}")

        in_array = array_info is not None
        in_palindrome = palindrome_info is not None

        if in_palindrome:
            category = "both" if in_array else "palindrome"
        elif in_array:
            category = "array"
        else:
            category = "neither"

        records.append({
            'Species': species,
            'GeneFamily': gene_family,
            'GeneID': gene_id,
            'Category': category,
            'InPalindrome': in_palindrome,
            'PalindromeName': palindrome_info['palindrome_name'] if palindrome_info else None,
            'PalindromeFullyContained': palindrome_info['fully_contained'] if palindrome_info else None,
            'PalindromePartialOverlap': palindrome_info['partial_overlap'] if palindrome_info else None,
            'InArray': in_array,
            'ArrayName': array_info['array_name'] if array_info else None,
            'ArrayStart': array_info['start'] if array_info else None,
            'ArrayEnd': array_info['end'] if array_info else None,
        })

    return pd.DataFrame(records)


def main():
    all_genes_df = fetch_google_sheet_by_gid(GID_S1, "S1. YAG genes accessions")
    all_genes_df = all_genes_df[all_genes_df['Use'] == 'yes']

    arrays_df = fetch_google_sheet_by_gid(GID_S14_ARRAY, "S14. Array", 1)
    arrays_df = arrays_df[arrays_df["species"].notna()].rename(columns={'species': 'Species'})

    palindromes_df = fetch_google_sheet_by_gid(GID_S13_PALINDROMES, "S13. Palindromes")

    genes_df = build_gene_table(all_genes_df, arrays_df, palindromes_df)

    print(genes_df.to_string(index=False))
    print(f"\nTotal genes: {len(genes_df)}")
    print(genes_df['Category'].value_counts().to_string())

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = f"gene_summary_{timestamp}.csv"
    genes_df.to_csv(output_path, index=False)
    print(f"\nSaved to {output_path}")


if __name__ == "__main__":
    main()
