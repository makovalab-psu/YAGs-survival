# Y-chromosome Selection Analysis Pipeline (Round 6)

This repository contains the pipeline and results for quantifying selection on Y-chromosome primate alignments using HyPhy.

## Pipeline Steps

### 1. Sequence De-duplication and Tree Pruning
The first step involves identifying identical sequences within each alignment and removing them to reduce computational overhead and avoid biases in selection analyses. The phylogenetic tree is pruned accordingly to match the unique sequences.

**Tool:** `hyphy-analyses/remove-duplicates/remove-duplicates.bf`

**Execution:**
For each gene directory, the following command was run:
```bash
hyphy /Users/sergei/Development/hyphy-analyses/remove-duplicates/remove-duplicates.bf \
      --msa [ALIGNMENT_FILE] \
      --tree [TREE_FILE] \
      --output [GENE].unique.nexus
```

**Outputs:**
- A `.unique.nexus` file for each gene containing the unique sequences and the pruned Newick tree.

### 2. Alignment Cleaning (Terminal Stop Codon Removal)
HyPhy selection analyses do not permit terminal stop codons. The `cln` tool is used to remove them.

**Tool:** `hyphy cln`

**Execution:**
For each gene directory, the following command was run:
```bash
hyphy cln --alignment [GENE].unique.nexus --output [GENE].cln.fasta
```

**Outputs:**
- A `.cln.fasta` file for each gene.

### Alignment Statistics

| Gene | Original Seqs | Unique Seqs | Codons |
|------|---------------|-------------|--------|
| BPY2 | 8 | 4 | 107 |
| CDY | 56 | 36 | 542 |
| CDY_CDYL | 63 | 43 | 540 |
| DAZ_DAZL_repeats | 256 | 54 | 30 |
| DAZ_DAZL_RRM | 73 | 36 | 165 |
| DAZ_repeats | 249 | 51 | 30 |
| DAZ_RRM | 66 | 31 | 165 |
| HSFY | 21 | 13 | 403 |
| RBMY | 76 | 50 | 400 |
| RBMY_RBMX | 83 | 55 | 400 |
| TSPY | 163 | 59 | 276 |
| VCY | 63 | 41 | 147 |

### Internal Stop Codon Statistics
These are premature stop codons found at positions other than the final codon of the sequence.

| Gene | Internal Stop Codons | Sequences with Internal Stops |
|------|-------------------|----------------------|
| BPY2 | 0 | 0 |
| CDY | 1 | 1 |
| CDY_CDYL | 1 | 1 |
| DAZ_DAZL_repeats | 0 | 0 |
| DAZ_DAZL_RRM | 0 | 0 |
| DAZ_repeats | 0 | 0 |
| DAZ_RRM | 0 | 0 |
| HSFY | 0 | 0 |
| RBMY | 8 | 8 |
| RBMY_RBMX | 8 | 8 |
| TSPY | 0 | 0 |
| VCY | 0 | 0 |

### 3. Selection Analysis (FitMG94)
A global FitMG94 model was used to estimate the mean dN/dS ratio across the entire tree and test for deviations from neutrality (dN/dS != 1).

**Tool:** `hyphy-analyses/FitMG94/FitMG94.bf`

**Execution:**
```bash
hyphy FitMG94.bf --alignment [GENE].unique.nexus --lrt Yes --output MG94fit.json
```

**Results:**

| Gene | dN/dS | 95% CI | p-value (dN/dS != 1) |
|------|-------|--------|----------------------|
| BPY2 | 0.630 | [0.195, 1.470] | 0.5613 |
| CDY | 0.725 | [0.626, 0.833] | 0.01393 |
| CDY_CDYL | 0.397 | [0.358, 0.438] | 0 |
| DAZ_DAZL_repeats | 0.966 | [0.759, 1.204] | 0.8696 |
| DAZ_DAZL_RRM | 0.599 | [0.498, 0.713] | 0.001031 |
| DAZ_repeats | 0.888 | [0.693, 1.114] | 0.6135 |
| DAZ_RRM | 0.637 | [0.519, 0.773] | 0.0102 |
| HSFY | 0.433 | [0.330, 0.557] | 0.0002236 |
| RBMY | 0.769 | [0.694, 0.847] | 0.01054 |
| RBMY_RBMX | 0.494 | [0.451, 0.540] | 1.11e-16 |
| TSPY | 0.768 | [0.640, 0.912] | 0.1181 |
| VCY | 1.048 | [0.871, 1.248] | 0.8022 |

### 4. BUSTED Model Selection
A model selection procedure was performed using BUSTED to identify the best-fitting number of rate classes for both non-synonymous (N) and synonymous (M) substitution rates. After determining the optimal rate classes, the model was further tested for improvement by allowing for multiple nucleotide substitutions (Multiple Hits: Double+Triple). The selection was based on minimizing the small-sample Akaike Information Criterion (c-AIC). Anomalies were checked by inspecting the weight of the "error absorption" rate class.

**Tool:** `hyphy-analyses/busted`

**Procedure:**
1.  Vary `rates` (N) and `syn-rates` (M), starting from N=1, M=1.
2.  If M=1, run with `--srv No`. If M>1, run with `--srv Yes --syn-rates M`.
3.  Iterate by increasing N or M as long as c-AIC improves.
4.  Once the best N and M are found, test if adding `--multiple-hits Double+Triple` further improves c-AIC.
5.  "Anomaly" is defined as a non-zero weight for the error sink component (dN/dS ~ 100).

**Results:**

| Gene | Rates (N) | Syn-Rates (M) | MH | c-AIC | Anomaly Score | Anomaly? | p-value |
|------|-----------|---------------|----|-------|---------------|----------|---------|
| BPY2 | 1 | 1 | None | 994.67 | 0.0000 | No | 0.5 |
| CDY | 2 | 2 | None | 8415.25 | 0.0000 | No | 0.5 |
| CDY_CDYL | 2 | 2 | None | 11631.37 | 0.0000 | No | 0.1598 |
| DAZ_DAZL_repeats | 1 | 1 | None | 1547.95 | 0.0000 | No | 0.5 |
| DAZ_DAZL_RRM | 2 | 2 | None | 3907.83 | 0.0000 | No | 0.5 |
| DAZ_repeats | 1 | 1 | None | 1514.13 | 0.0000 | No | 0.5 |
| DAZ_RRM | 2 | 2 | None | 3460.88 | 0.0000 | No | 0.5 |
| HSFY | 2 | 1 | None | 4451.90 | 0.0000 | No | 0.5 |
| RBMY | 2 | 2 | Double+Triple | 9493.91 | 0.0000 | No | 0.5 |
| RBMY_RBMX | 2 | 2 | None | 11439.60 | 0.0000 | No | 0.0316 |
| TSPY | 1 | 2 | None | 4769.46 | 0.0000 | No | 0.5 |
| VCY | 1 | 2 | None | 3216.28 | 0.0000 | No | 0.4995 |

### 5. Branch Heterogeneity (Global vs Local MG94)
To test for variation in selection pressure across lineages, a local MG94 model (separate dN/dS for each branch) was fitted and compared to the global MG94 model (single dN/dS) using a Likelihood Ratio Test (LRT).

**Tool:** `hyphy-analyses/FitMG94/FitMG94.bf --type local`

**Results:**

| Gene | Global LogL | Global Params | Local LogL | Local Params | LRT | df | p-value |
|------|-------------|---------------|------------|--------------|-----|----|---------|
| BPY2 | -473.91 | 20 | -472.43 | 24 | 2.97 | 4 | 0.5621 |
| CDY | -4135.00 | 76 | -4098.82 | 136 | 72.36 | 60 | 0.1317 |
| CDY_CDYL | -5739.44 | 90 | -5664.81 | 164 | 149.27 | 74 | 5.23e-07 |
| DAZ_DAZL_repeats | -668.45 | 97 | -633.62 | 178 | 69.66 | 81 | 0.8114 |
| DAZ_DAZL_RRM | -1909.83 | 80 | -1865.98 | 144 | 87.70 | 64 | 0.02632 |
| DAZ_repeats | -653.30 | 95 | -620.18 | 174 | 66.25 | 79 | 0.8464 |
| DAZ_RRM | -1690.52 | 70 | -1654.48 | 124 | 72.09 | 54 | 0.05055 |
| HSFY | -2187.66 | 34 | -2178.83 | 52 | 17.68 | 18 | 0.477 |
| RBMY | -4673.60 | 104 | -4635.22 | 192 | 76.77 | 88 | 0.7981 |
| RBMY_RBMX | -5638.80 | 113 | -5548.82 | 210 | 179.96 | 97 | 6.427e-07 |
| TSPY | -2281.71 | 104 | -2244.86 | 192 | 73.71 | 88 | 0.8623 |
| VCY | -1533.55 | 83 | -1506.67 | 150 | 53.76 | 67 | 0.8791 |

### 6. Episodic Branch Selection (aBSREL)
aBSREL was used to test for episodic positive selection on specific branches of the phylogeny.

**Tool:** `hyphy absrel`

**Results:**

| Gene | Selected Branches | Tested Branches |
|------|-------------------|-----------------|
| BPY2 | 0 | 5 |
| CDY | 0 | 61 |
| CDY_CDYL | 0 | 75 |
| DAZ_DAZL_repeats | 0 | 82 |
| DAZ_DAZL_RRM | 0 | 65 |
| DAZ_repeats | 0 | 80 |
| DAZ_RRM | 0 | 55 |
| HSFY | 0 | 19 |
| RBMY | 1 | 89 |
| RBMY_RBMX | 1 | 98 |
| TSPY | 0 | 89 |
| VCY | 0 | 68 |

### 7. Episodic Selection at Sites (MEME)
MEME (Mixed Effects Model of Evolution) was used to detect individual sites subject to episodic positive selection (p < 0.05).

**Tool:** `hyphy meme`

**Results:**

| Gene | Selected Sites | Total Sites |
|------|----------------|-------------|
| BPY2 | 0 | 107 |
| CDY | 1 | 542 |
| CDY_CDYL | 6 | 540 |
| DAZ_DAZL_repeats | 1 | 30 |
| DAZ_DAZL_RRM | 1 | 165 |
| DAZ_repeats | 1 | 30 |
| DAZ_RRM | 1 | 165 |
| HSFY | 0 | 403 |
| RBMY | 8 | 400 |
| RBMY_RBMX | 15 | 400 |
| TSPY | 2 | 276 |
| VCY | 2 | 147 |

### 8. Summary
A comprehensive summary of all analysis results, including alignment statistics, stop codon counts, and selection test p-values, has been generated in HTML format.

**Output:** `summary_table.html`

---
*Last updated: Tuesday, January 27, 2026*
