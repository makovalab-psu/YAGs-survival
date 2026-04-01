#!/usr/bin/env python3
"""
Permutation test: all-by-all comparison of similarity categories.
Matrix output: upper triangle = median p-values, lower triangle = mean p-values.
Violin plot output: 5 categories with mean (blue dashed) and median (orange dashed) lines.
"""

import os
import glob
from datetime import datetime
from typing import Tuple
from itertools import combinations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import FixedLocator

OUTPUT_DIR = "./output"
N_PERMUTATIONS = 10_000

CATEGORIES = [
    ("Tandem\nrepeat\n(same array)", "tandem_same"),
    ("Tandem\nrepeat\n(different array)", "tandem_diff"),
    ("Palindrome\nopposite\narms", "pal_opposite_arms"),
    ("Palindrome\nsame arm", "pal_same_arm"),
    ("Palindrome\ndifferent", "pal_diff_q"),
    ("Outside", "outside"),
]


def load_combined_data() -> pd.DataFrame:
    files = glob.glob(os.path.join(OUTPUT_DIR, "*COMBINED*.csv"))
    print(f"Found {len(files)} COMBINED files")
    if not files:
        return pd.DataFrame()
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    print(f"Total comparisons after filter: {len(df)}")
    return df


def get_category_mask(df: pd.DataFrame, tag: str) -> pd.Series:
    if tag == "tandem_same":
        return (df["Type"] == "ARRAY") & (df["SameArrayRegion"] == True)
    elif tag == "tandem_diff":
        return (df["Type"] == "ARRAY") & (df["SameArrayRegion"] == False)
    elif tag == "pal_opposite_arms":
        return (df["PalindromeTag"] == "palindrome_sameQ_trans") & (df["BothGenesFullyInPalindrome"] == True)
    elif tag == "pal_same_arm":
        return (df["Type"] == "PALINDROME") & (df["PalindromeTag"] == "palindrome_sameQ_cis")
    elif tag == "pal_diff_q":
        return (df["Type"] == "PALINDROME") & (df["PalindromeTag"] == "palindrome_diffQ")
    else:  # outside
        return df["Type"] == "OTHER"


def get_category_data(df: pd.DataFrame, tag: str) -> np.ndarray:
    return df.loc[get_category_mask(df, tag), "Similarity"].values


def get_category_lengths(df: pd.DataFrame, tag: str) -> np.ndarray:
    mask = get_category_mask(df, tag)
    genes_a = df.loc[mask, ["GeneA", "GeneA_length"]].rename(columns={"GeneA": "Gene", "GeneA_length": "Length"})
    genes_b = df.loc[mask, ["GeneB", "GeneB_length"]].rename(columns={"GeneB": "Gene", "GeneB_length": "Length"})
    return pd.concat([genes_a, genes_b]).drop_duplicates(subset="Gene")["Length"].values


def count_category_genes(df: pd.DataFrame, tag: str) -> int:
    mask = get_category_mask(df, tag)
    return pd.concat([df.loc[mask, "GeneA"], df.loc[mask, "GeneB"]]).nunique()


def permutation_test(group1: np.ndarray, group2: np.ndarray,
                     stat_func, n_permutations: int = N_PERMUTATIONS
                     ) -> Tuple[float, float]:
    observed_diff = stat_func(group1) - stat_func(group2)
    combined = np.concatenate([group1, group2])
    n1 = len(group1)
    perm_diffs = np.zeros(n_permutations)
    for i in range(n_permutations):
        np.random.shuffle(combined)
        perm_diffs[i] = stat_func(combined[:n1]) - stat_func(combined[n1:])
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))
    return observed_diff, p_value


def plot_length_distributions(df: pd.DataFrame, timestamp: str) -> None:
    labels = [c[0] for c in CATEGORIES]
    COLORS = ["#e07b39", "#c4a23a", "#6aab6a", "#5b8ec4", "#c45b8e", "#8e5bc4"]
    VIOLIN_WIDTH = 0.7
    LINE_HALF = VIOLIN_WIDTH / 2

    length_data = [get_category_lengths(df, tag) for _, tag in CATEGORIES]

    fig, ax = plt.subplots(figsize=(6, 6))
    n = len(labels)

    parts = ax.violinplot(length_data, positions=range(n), widths=VIOLIN_WIDTH,
                          showmeans=False, showmedians=False, showextrema=True)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(COLORS[i])
        pc.set_edgecolor('none')
        pc.set_alpha(0.6)
    for part_name in ('cbars', 'cmins', 'cmaxes'):
        if part_name in parts:
            parts[part_name].set_edgecolor('#555555')
            parts[part_name].set_linewidth(0.8)

    for i, vals in enumerate(length_data):
        if len(vals) == 0:
            continue
        mean_val = np.mean(vals)
        median_val = np.median(vals)
        ax.plot([i - LINE_HALF, i + LINE_HALF], [mean_val, mean_val], color='black', lw=1, ls='--')
        ax.annotate(f"mean {mean_val:.0f}", xy=(i - LINE_HALF, mean_val),
                    xytext=(0, -3), textcoords='offset points',
                    ha='left', va='top', fontsize=7, color='black')
        ax.plot([i - LINE_HALF, i + LINE_HALF], [median_val, median_val], color='black', lw=1, ls='-')
        ax.annotate(f"median {median_val:.0f}", xy=(i + LINE_HALF, median_val),
                    xytext=(0, 3), textcoords='offset points',
                    ha='right', va='bottom', fontsize=7, color='black', fontweight='bold')

    ax.set_xticks(range(n))
    ax.set_xticklabels(
        [f"{l}\n(n={len(length_data[i])})" for i, l in enumerate(labels)],
        fontsize=6
    )
    ax.set_ylabel("Gene length (bp)", fontsize=7)
    ax.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()

    path = os.path.join(OUTPUT_DIR, f"gene_length_distributions_{timestamp}.pdf")
    plt.savefig(path, dpi=1200, bbox_inches="tight")
    plt.close()
    print(f"Length distribution plot saved: {path}")


def significance_label(p: float) -> str:
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return "ns"


def main():
    print("=" * 60)
    print("Permutation Tests: All-by-All Category Comparison")
    print("=" * 60)

    df = load_combined_data()
    # filter comparisons between genes where length difference is smaller than 10%
    df = df[df['Relative_lenght_difference'] <= 0.1]

    if df.empty:
        print("ERROR: No COMBINED CSV files found.")
        return

    labels = [c[0] for c in CATEGORIES]

    n = len(labels)

    data = {}
    for label, tag in CATEGORIES:
        vals = get_category_data(df, tag)
        data[label] = vals
        n_genes = count_category_genes(df, tag)
        print(f"  {label:20s}  comparisons={len(vals):5d}  genes={n_genes:4d}  median={np.median(vals):.3f}  mean={np.mean(vals):.3f}")

    # Matrices for p-values and diffs
    median_p = np.full((n, n), np.nan)
    mean_p = np.full((n, n), np.nan)
    median_diff = np.full((n, n), np.nan)
    mean_diff = np.full((n, n), np.nan)

    for i, j in combinations(range(n), 2):
        g1, g2 = data[labels[i]], data[labels[j]]
        if len(g1) == 0 or len(g2) == 0:
            continue
        print(f"\nRunning: {labels[i]} vs {labels[j]}...")

        md, mp = permutation_test(g1, g2, np.median)
        median_diff[i, j] = median_diff[j, i] = md
        median_p[i, j] = median_p[j, i] = mp

        mnd, mnp = permutation_test(g1, g2, np.mean)
        mean_diff[i, j] = mean_diff[j, i] = mnd
        mean_p[i, j] = mean_p[j, i] = mnp

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # 1. Save simple CSV table
    rows = []
    for i, j in combinations(range(n), 2):
        g1, g2 = data[labels[i]], data[labels[j]]
        rows.append({
            "Group1": labels[i].replace('\n', ' '),
            "Group2": labels[j].replace('\n', ' '),
            "n1": len(g1),
            "n2": len(g2),
            "median1": np.median(g1),
            "median2": np.median(g2),
            "median_diff": median_diff[i, j],
            "median_p": median_p[i, j],
            "median_sig": significance_label(median_p[i, j]),
            "mean1": np.mean(g1),
            "mean2": np.mean(g2),
            "mean_diff": mean_diff[i, j],
            "mean_p": mean_p[i, j],
            "mean_sig": significance_label(mean_p[i, j]),
        })
    results_df = pd.DataFrame(rows)
    csv_path = os.path.join(OUTPUT_DIR, f"permutation_test_v3_results_{timestamp}.csv")
    results_df.to_csv(csv_path, index=False)
    print(f"\nResults CSV saved: {csv_path}")

    # Print table
    print(f"\n{results_df.to_string(index=False)}")

    # 2. Heatmap: upper triangle = mean, lower triangle = median
    display_matrix = np.full((n, n), np.nan)
    p_matrix = np.full((n, n), np.nan)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            elif j > i:
                display_matrix[i, j] = mean_diff[i, j]
                p_matrix[i, j] = mean_p[i, j]
            else:
                display_matrix[i, j] = median_diff[j, i]
                p_matrix[i, j] = median_p[j, i]

    sig_matrix = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            p = p_matrix[i, j]
            if p < 0.001:
                sig_matrix[i, j] = 3
            elif p < 0.01:
                sig_matrix[i, j] = 2
            elif p < 0.05:
                sig_matrix[i, j] = 1
            else:
                sig_matrix[i, j] = 0

    # colors = ["#e0e0e0", "#fff9c4", "#ffcc80", "#ef5350"]
    # # cmap = mcolors.ListedColormap(colors)
    # # bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    # # # norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # fig, ax = plt.subplots(figsize=(20, 8))
    # im = ax.imshow(sig_matrix, cmap=cmap, norm=norm, aspect='equal')

    # for i in range(n):
    #     ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=True, color='white', zorder=2))

    # for i in range(n):
    #     for j in range(n):
    #         if i == j:
    #             continue
    #         diff = display_matrix[i, j]
    #         p = p_matrix[i, j]
    #         sig = significance_label(p)
    #         text = f"{diff:+.3f}\n{sig} (p={p:.4f})"
    #         ax.text(j, i, text, ha='center', va='center', fontsize=7.5, zorder=3)

    # ax.set_xticks([])
    # ax.set_yticks([])
    # for spine in ax.spines.values():
    #     spine.set_visible(False)

    # for i, label in enumerate(labels):
    #     vals = data[label]
    #     med, mean_val = np.median(vals), np.mean(vals)
    #     full = f"{label}\n{med:.1f} ({mean_val:.1f})"
    #     ax.text(i, -0.55, full, ha='center', va='bottom', fontsize=9, rotation=90,
    #             multialignment='center', clip_on=False)

    # for i, label in enumerate(labels):
    #     vals = data[label]
    #     med, mean_val = np.median(vals), np.mean(vals)
    #     full = f"{label}\n{med:.1f} ({mean_val:.1f})"
    #     ax.text(-0.55, i, full, ha='right', va='center', fontsize=9,
    #             multialignment='left', clip_on=False)

    # ax.plot([-0.8, n - 0.2], [-0.8, n - 0.2], 'k--', lw=1.5, clip_on=False, zorder=5)
    # ax.text(n - 0.2, n, r"Permutation test $\bf{median}$",
    #         ha='right', va='top', fontsize=12, fontstyle='normal', rotation=0, clip_on=False)
    # ax.text(n, n - 0.2, r"Permutation test $\bf{mean}$",
    #         ha='left', va='bottom', fontsize=12, fontstyle='normal', rotation=90, clip_on=False)

    # for i in range(n + 1):
    #     ax.axhline(i - 0.5, color='white', lw=2)
    #     ax.axvline(i - 0.5, color='white', lw=2)

    # cbar = fig.colorbar(im, ax=ax, orientation='vertical', ticks=[0, 1, 2, 3],
    #                     shrink=0.6, pad=0.04, aspect=20)
    # cbar.ax.set_yticklabels(["ns", "p<0.05 *", "p<0.01 **", "p<0.001 ***"])
    # cbar.outline.set_visible(False)

    # plt.tight_layout()
    # fig_path = os.path.join(OUTPUT_DIR, f"permutation_test_v3_heatmap_{timestamp}.pdf")
    # plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    # plt.close()
    # print(f"Heatmap saved: {fig_path}")

    # 3. Violin plot: 5 categories side by side
    COLORS = [ "#7bb0df", "#1964b0",  "#f4a637",  "#db5829", "#894b45" , "#b2b2b2"]
    VIOLIN_WIDTH = 0.7
    LINE_HALF = VIOLIN_WIDTH / 2


    # Conversion factor: 1 inch = 25.4 mm
    MM_PER_INCH = 25.4
    fig_width_mm = 120
    fig_height_mm = 85


    fig, ax = plt.subplots(figsize=(fig_width_mm/MM_PER_INCH, fig_height_mm/MM_PER_INCH))
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"

    violin_data = [data[label] for label in labels]
    parts = ax.violinplot(violin_data, positions=range(n), widths=VIOLIN_WIDTH,
                          showmeans=False, showmedians=False, showextrema=True)

    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(COLORS[i])
        pc.set_edgecolor('none')
        pc.set_alpha(0.75)
    for part_name in ('cbars', 'cmins', 'cmaxes'):
        if part_name in parts:
            parts[part_name].set_edgecolor('#555555')
            parts[part_name].set_linewidth(0.5)

    for i, label in enumerate(labels):
        vals = data[label]
        mean_val = np.mean(vals)
        median_val = np.median(vals)

        # Mean: thin dashed blue line, label left-aligned with line, 3pt below
        ax.plot([i - LINE_HALF, i + LINE_HALF], [mean_val, mean_val],
                color='black', lw=0.5, ls='--')
        ax.annotate(f"mean {mean_val:.2f}%", xy=(i - LINE_HALF, mean_val),
                    xytext=(0, -3), textcoords='offset points',
                    ha='left', va='top', fontsize=5, color='black')

        # Median: thin dashed orange line, label right-aligned with line, 3pt above
        ax.plot([i - LINE_HALF, i + LINE_HALF], [median_val, median_val],
                color='black', lw=0.5, ls='-')
        ax.annotate(f"median {median_val:.2f}%", xy=(i + LINE_HALF, median_val),
                    xytext=(0, 3), textcoords='offset points',
                    ha='right', va='bottom', fontsize=5, color='black', fontweight='bold' )

    # Significance brackets for median p-values
    sig_pairs = sorted(
        [(i, j, median_p[i, j], significance_label(median_p[i, j]))
         for i, j in combinations(range(n), 2)
         if not np.isnan(median_p[i, j]) and median_p[i, j] < 0.05],
        key=lambda x: (x[0], x[1])
    )
    sig_pairs = reversed(sig_pairs)

    if sig_pairs:
        y_max = max(np.max(v) for v in violin_data if len(v) > 0) +0.7
        y_min = min(np.min(v) for v in violin_data if len(v) > 0) 
        step = (y_max - y_min) * 0.04
        # tick_h = step * 0.3

        col_heights = [y_max + step * 1.5] * n 

        for i, j, _, stars in sig_pairs:
            y = max(col_heights[i:j + 1]) 
            mid = (i + j) / 2
            gap = len(stars) * 0.08
            ax.plot([i, mid - gap], [y, y], color='#333333', lw=0.5)
            ax.plot([mid + gap, j], [y, y], color='#333333', lw=0.5)
            # ax.plot([i, i, mid - gap], [y - tick_h, y, y], color='#333333', lw=0.8)
            # ax.plot([mid + gap, j, j], [y, y, y - tick_h], color='#333333', lw=0.8)
            ax.text(mid, y - 0.15, stars, ha='center', va='center', fontsize=5, color='#333333')
            for k in range(i, j + 1):
                col_heights[k] = y + step

        ax.set_ylim(bottom=87, top=max(col_heights) + step +2)
    else:
        ax.set_ylim(bottom=87, top=101)
    # ax.labelsize(5)
    ax.tick_params(axis='y', labelsize =5)
    # ax.set_yticks(range(87, 101, 5))
    ax.set_yticks([87, 90, 95, 100])
    ax.set_yticklabels(['87', '90', '95', '100'])

    # Minor ticks at every integer, no labels
    ax.yaxis.set_minor_locator(FixedLocator(range(87, 100)))
    ax.spines['left'].set_bounds(87, 100)
    ax.set_xticks(range(n))
    ax.set_xticklabels(
        [f"{l}" for l in labels],
        fontsize=6
    )
    ax.set_ylabel("Sequence identity (%)", fontsize=7)
    ax.spines[['right', 'top']].set_visible(False)

    ax.text(0.01, 0.01,
            "* p < 0.05\n** p < 0.01\n*** p < 0.001",
            transform=ax.transAxes, fontsize=5, va='bottom', ha='left',
            linespacing=1.8)

    plt.tight_layout()
    violin_path = os.path.join(OUTPUT_DIR, f"permutation_test_violins_{timestamp}.pdf")
    plt.savefig(violin_path, dpi=1200, bbox_inches="tight")
    plt.close()
    print(f"Violin plot saved: {violin_path}")

    plot_length_distributions(df, timestamp)

    print(f"\n{'=' * 60}")
    print("Analysis complete!")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
