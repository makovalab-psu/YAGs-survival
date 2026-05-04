#!/usr/bin/env python3
"""
Summary table: rows=species, columns=gene families.
Cell format: Xg/Np  Yg/Na
  Xg = genes in palindromes, Np = unique palindromes (main Q number)
  Yg = genes in arrays,      Na = unique array regions
"""

import re
from datetime import datetime

import pandas as pd

from analyze_array_palindrome import (
    fetch_google_sheet_by_gid,
    GID_S1, GID_S13_PALINDROMES, GID_S14_ARRAY,
)
from gene_summary import build_gene_table


def main_q(palindrome_name: str) -> str:
    """Extract main Q number from palindrome name, e.g. 'Q1' from 'NC_...Q1.1.A'."""
    match = re.search(r'\.Q(\d+)', palindrome_name)
    return f"Q{match.group(1)}" if match else palindrome_name


def format_cell(group: pd.DataFrame) -> str:
    pal = group[group['InPalindrome']]
    arr = group[group['InArray']]
    n_pal_genes = len(pal)
    n_arr_genes = len(arr)
    n_pal_regions = pal['PalindromeName'].dropna().map(main_q).nunique()
    n_arr_regions = arr[['ArrayStart', 'ArrayEnd']].dropna().drop_duplicates().shape[0]
    return f"{n_pal_genes}g/{n_pal_regions}p  {n_arr_genes}g/{n_arr_regions}a"


def main():
    all_genes_df = fetch_google_sheet_by_gid(GID_S1, "S1. YAG genes accessions")
    all_genes_df = all_genes_df[all_genes_df['Use'] == 'yes']

    arrays_df = fetch_google_sheet_by_gid(GID_S14_ARRAY, "S14. Array", 1)
    arrays_df = arrays_df[arrays_df["species"].notna()].rename(columns={'species': 'Species'})

    palindromes_df = fetch_google_sheet_by_gid(GID_S13_PALINDROMES, "S13. Palindromes")

    genes_df = build_gene_table(all_genes_df, arrays_df, palindromes_df)

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

    summary = {'Species': 'Total'}
    for gf in gene_families:
        summary[gf] = format_cell(genes_df[genes_df['GeneFamily'] == gf])
    summary['Total'] = format_cell(genes_df)
    rows.append(summary)

    table = pd.DataFrame(rows).set_index('Species')

    print("\nPalindrome & Array Summary  |  format: Xg/Np  Yg/Na")
    print("  Xg = genes in palindromes, Np = unique palindromes (Q number)")
    print("  Yg = genes in arrays,      Na = unique array regions")
    print("=" * 80)
    print(table.to_string())

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = f"gene_palindrome_array_summary_{timestamp}.csv"
    table.to_csv(output_path)
    print(f"\nSaved to {output_path}")


if __name__ == "__main__":
    main()
