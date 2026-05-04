import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- Load copy number medians ---
df = pd.read_csv("human_copies.csv")
gene_cols = [c for c in df.columns if c not in ("species", "subspecies", "individual")]
medians: dict[str, float] = {gene: df[gene].median() for gene in gene_cols}

# --- Load isoform counts for HomSap ---
with open("results.json") as f:
    results = json.load(f)

isoform_counts: dict[str, int] = {
    d["gene"]: len(d["isoforms"])
    for d in results
    if d["species"] == "HomSap"
}

# --- Intersect genes present in both ---
genes = sorted(set(medians) & set(isoform_counts))
x = [isoform_counts[g] for g in genes]
y = [medians[g] for g in genes]

# --- Plot ---
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(x, y, s=80, color="steelblue", zorder=3)

for gene, xi, yi in zip(genes, x, y):
    ax.annotate(gene, (xi, yi), textcoords="offset points", xytext=(6, 4), fontsize=9)

ax.set_xlabel("Number of structural isoforms (HomSap)")
ax.set_ylabel("Median copy number (HomSap individuals)")
ax.set_title("Structural isoform diversity vs. copy number")
ax.grid(True, linestyle="--", alpha=0.4)

Path("output").mkdir(exist_ok=True)
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
out = f"output/copy_vs_isoforms_{timestamp}.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
plt.show()
