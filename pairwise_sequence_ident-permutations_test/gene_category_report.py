#!/usr/bin/env python3
"""
Gene category summary: rows=species, columns=gene families.
Cell format: X Y (w) Z
  X   = genes in palindrome only
  Y   = genes in array only
  (w) = genes in both palindrome and array
  Z   = genes in neither
"""

import os
import pandas as pd

from analyze_array_palindrome import (
    fetch_google_sheet_by_gid,
    parse_gff_for_gene,
    find_array_for_position,
    GID_S1, GID_S13_PALINDROMES, GID_S14_ARRAY,
    GFF_DOWNLOAD_DIR,
)
from gene_summary import find_palindrome_any_overlap

OUTPUT_PATH = "gene_category_report.csv"


def classify_genes(all_genes_df: pd.DataFrame, arrays_df: pd.DataFrame,
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

        in_array = find_array_for_position(gene_start, gene_info['chrom'], gene_family, sp_arrays) is not None
        in_palindrome = find_palindrome_any_overlap(gene_start, gene_end, gene_info['chrom'], sp_palindromes) is not None

        records.append({
            'Species': species,
            'GeneFamily': gene_family,
            'in_palindrome': in_palindrome,
            'in_array': in_array,
        })
    return pd.DataFrame(records)


def format_cell(group: pd.DataFrame) -> str:
    pal   = (group['in_palindrome'] & ~group['in_array']).sum()
    arr   = (~group['in_palindrome'] & group['in_array']).sum()
    both  = (group['in_palindrome'] & group['in_array']).sum()
    neither = (~group['in_palindrome'] & ~group['in_array']).sum()
    return f"{pal} {arr} ({both}) {neither}"


def main():
    all_genes_df = fetch_google_sheet_by_gid(GID_S1, "S1. YAG genes accessions")
    all_genes_df = all_genes_df[all_genes_df['Use'] == 'yes']

    arrays_df = fetch_google_sheet_by_gid(GID_S14_ARRAY, "S14. Array", 1)
    arrays_df = arrays_df[arrays_df["species"].notna()].rename(columns={'species': 'Species'})

    palindromes_df = fetch_google_sheet_by_gid(GID_S13_PALINDROMES, "S13. Palindromes")

    genes_df = classify_genes(all_genes_df, arrays_df, palindromes_df)

    species_list = sorted(genes_df['Species'].unique())
    gene_families = sorted(genes_df['GeneFamily'].unique())

    rows = []
    for species in species_list:
        row = {'Species': species}
        for gf in gene_families:
            group = genes_df[(genes_df['Species'] == species) & (genes_df['GeneFamily'] == gf)]
            row[gf] = format_cell(group) if len(group) > 0 else "-"
        row['Total'] = format_cell(genes_df[genes_df['Species'] == species])
        rows.append(row)

    # Summary row across all species per gene family
    summary = {'Species': 'Total'}
    for gf in gene_families:
        summary[gf] = format_cell(genes_df[genes_df['GeneFamily'] == gf])
    summary['Total'] = format_cell(genes_df)
    rows.append(summary)

    table = pd.DataFrame(rows).set_index('Species')

    print("\nGene Category Summary  |  format: palindrome  array  (both)  neither")
    print("=" * 80)
    print(table.to_string())

    table.to_csv(OUTPUT_PATH)
    print(f"\nSaved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
