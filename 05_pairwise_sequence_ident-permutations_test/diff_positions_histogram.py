#!/usr/bin/env python3
"""
For each of the 5 comparison categories, plot a histogram of mutation positions
(Diff_positions column) — where along the ungapped alignment each mismatch occurs,
expressed as % through the gene.
"""

import ast
import os
import glob
from datetime import datetime
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUTPUT_DIR = "./output"

CATEGORIES = [
    ("Tandem\nrepeat",           "ARRAY"),
    ("Palindrome\nopposite\narms", "palindrome_sameQ_trans"),
    ("Palindrome\nsame arm",     "palindrome_sameQ_cis"),
    ("Palindrome\ndifferent",    "palindrome_diffQ"),
    ("Outside",                  "OTHER"),
]

COLORS = ["#e07b39", "#6aab6a", "#5b8ec4", "#c45b8e", "#8e5bc4"]


def load_combined_data() -> pd.DataFrame:
    files = glob.glob(os.path.join(OUTPUT_DIR, "*COMBINED*.csv"))
    print(f"Found {len(files)} COMBINED files")
    if not files:
        return pd.DataFrame()
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    print(f"Total comparisons loaded: {len(df)}")
    return df


def parse_positions(val) -> List[float]:
    if pd.isna(val) or val == "" or val == "[]":
        return []
    if isinstance(val, list):
        return val
    try:
        return ast.literal_eval(val)
    except (ValueError, SyntaxError):
        return []


def get_category_mask(df: pd.DataFrame, tag: str) -> pd.Series:
    if tag == "palindrome_sameQ_trans":
        return (df["PalindromeTag"] == tag) & (df["BothGenesFullyInPalindrome"] == True)
    elif tag in ("palindrome_sameQ_cis", "palindrome_diffQ"):
        return (df["Type"] == "PALINDROME") & (df["PalindromeTag"] == tag)
    elif tag == "ARRAY":
        return (df["Type"] == "ARRAY") & (df["SameArrayRegion"] == True)
    else:
        return df["Type"] == "OTHER"


def get_positions(df: pd.DataFrame, tag: str) -> np.ndarray:
    mask = get_category_mask(df, tag)
    all_pos = []
    for val in df.loc[mask, "Diff_positions"]:
        all_pos.extend(parse_positions(val))
    return np.array(all_pos)


def main():
    df = load_combined_data()
    df = df[df["Relative_lenght_difference"] <= 0.1]

    if df.empty:
        print("ERROR: No COMBINED CSV files found.")
        return

    fig, axes = plt.subplots(1, 5, figsize=(16, 4), sharey=False)

    bins = np.linspace(0, 100, 21)  # 20 bins of 5% each

    for ax, (label, tag), color in zip(axes, CATEGORIES, COLORS):
        positions = get_positions(df, tag)
        n_comparisons = get_category_mask(df, tag).sum()

        ax.hist(positions, bins=bins, color=color, alpha=0.75, edgecolor="white", linewidth=0.5)
        ax.set_title(label, fontsize=9, multialignment="center")
        ax.set_xlabel("Position in gene (%)", fontsize=8)
        ax.set_ylabel("# mutations", fontsize=8)
        ax.set_xlim(0, 100)
        ax.spines[["right", "top"]].set_visible(False)
        ax.tick_params(labelsize=7)
        ax.text(0.97, 0.97, f"n={len(positions)}\n({n_comparisons} pairs)",
                transform=ax.transAxes, fontsize=7, ha="right", va="top", color="#444444")

    plt.suptitle("Distribution of mutation positions along the gene", fontsize=11, y=1.02)
    plt.tight_layout()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = os.path.join(OUTPUT_DIR, f"diff_positions_histogram_{timestamp}.pdf")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
