#!/usr/bin/env bash
# Selection-on-domains pipeline wrapper
# All commands are run from this directory.
#
# STEP OVERVIEW
# -------------
#   align_domains.py  [CDY|RBMY|HSFY|TSPY]   step 00 — MAFFT domain alignments
#   00_import_alignments.py                   step 00b — copy pre-built alignments from test_alignments/
#   01_clean_stops.py                         replace in-frame stop codons with ---
#   02_deduplicate.py                         remove exact duplicate sequences
#   03_build_trees.py                         IQtree2 ML trees with outgroup rooting
#   04_fitmg94.py                             MG94 global + local dN/dS (LRT)
#   05_busted.py                              BUSTED episodic selection (AIC-optimised)
#   06_absrel.py                              aBSREL branch-specific selection
#   07_meme.py                                MEME site-specific selection
#   08_summarize.py                           collect all results → output/summary.tsv
  09_plot_trees.R / 09_plot_trees.py        one tree plot per domain (PDF + PNG)

VENV="../venv/bin/activate"
source "$VENV"

# ─── FULL PIPELINE ──────────────────────────────────────────────────────────

# Build all domain alignments from scratch, then run the full selection pipeline:
# python align_domains.py CDY
# python align_domains.py RBMY
# python align_domains.py HSFY
# python align_domains.py TSPY
# python run_all.py 01 02 03 04 05 06 07 08

# Import pre-built alignments, then run the full pipeline:
# python run_all.py 00b 01 02 03 04 05 06 07 08


# ─── PARTIAL RUNS ───────────────────────────────────────────────────────────

# Re-run a single gene family's domain alignment (overwrites output/CDY_*/alignment.fasta):
# python align_domains.py CDY

# Run only the prep steps (clean + dedup) — fast, safe to re-run:
# python run_all.py 01 02

# Build trees only (after alignments are ready):
# python run_all.py 03

# Run all HyPhy analyses without re-building trees:
# python run_all.py 04 05 06 07

# Run a single analysis step across all domains:
# python 05_busted.py
# python 06_absrel.py

# Run specific steps for one domain only (substring match on directory name):
# python run_all.py 05 06 07 --only DAZ
# python run_all.py 03 04 05 --only CDY_Chrom

# Regenerate summary table from existing JSON results (no re-running):
# python run_all.py 08
# python 08_summarize.py   # same thing, directly


# ─── ACTIVE COMMAND ─────────────────────────────────────────────────────────

python run_all.py 03
