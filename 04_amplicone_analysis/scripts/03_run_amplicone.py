#!/usr/bin/env python3
"""
03_run_amplicone.py
===================
Batch-run AmpliCoNE-count.py on all individuals with mapped reads available.

For each species/individual with a merged.100.bam file:
  1. Creates the output directory if absent.
  2. Indexes the BAM file if no .bai index exists.
  3. Runs AmpliCoNE-count.py from within the individual's output directory
     (AmpliCoNE writes its output files to the current working directory).

All species parameters (accession, length, gene definition filename,
annotation filename, alignments root) are read from config.yaml.

AmpliCoNE-count.py outputs per individual:
  OutputAmpliconic_Summary.txt  — estimated copy number per gene family
  OutputXDG_CopyNumber.txt      — depth at single-copy control regions

Usage:
    python3 03_run_amplicone.py \\
        --config          /path/to/config.yaml \\
        --panels_dir      /path/to/panels \\
        --gene_defs_dir   /path/to/gene_definitions/targets \\
        --results_root    /path/to/results \\
        --amplicone_tool  /path/to/AmpliCoNE-tool
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import yaml


def load_config(config_path: str) -> dict:
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def run_individual(
    species: str,
    individual: str,
    bam: str,
    out_dir: str,
    panels_dir: str,
    gene_defs_dir: str,
    amplicone_tool: str,
    sp_cfg: dict,
) -> None:
    annotation = os.path.join(panels_dir, sp_cfg["annotation_tab"])
    gene_defs  = os.path.join(gene_defs_dir, sp_cfg["gene_definitions"])
    count_py   = os.path.join(amplicone_tool, "AmpliCoNE-count.py")

    os.makedirs(out_dir, exist_ok=True)

    # Index BAM if needed
    if not os.path.exists(bam + ".bai"):
        print(f"  Indexing {bam} ...")
        subprocess.run(["samtools", "index", bam], check=True)

    cmd = (
        f"cd {out_dir} && "
        f"python3 {count_py} "
        f"--LENGTH {sp_cfg['chry_length']} "
        f"--CHR {sp_cfg['chry_accession']} "
        f"--GENE_DEF {gene_defs} "
        f"--ANNOTATION {annotation} "
        f"--BAM {bam}"
    )
    print(f"  Running AmpliCoNE-count for {species}/{individual} ...")
    subprocess.run(cmd, shell=True, check=True)


def main():
    # Default config path: two levels up from this script (pipeline root)
    default_config = Path(__file__).parent.parent / "config.yaml"

    parser = argparse.ArgumentParser(description="Batch AmpliCoNE copy-number estimation.")
    parser.add_argument("--config", default=str(default_config),
                        help="Path to config.yaml (default: ../config.yaml relative to this script)")
    parser.add_argument("--panels_dir", required=True,
                        help="Directory containing annotation .tab files")
    parser.add_argument("--gene_defs_dir", required=True,
                        help="Directory containing gene definition .tsv files")
    parser.add_argument("--results_root", required=True,
                        help="Root directory for AmpliCoNE output: {species}/{individual}/")
    parser.add_argument("--amplicone_tool", required=True,
                        help="Path to AmpliCoNE-tool directory (contains AmpliCoNE-count.py)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    alignments_root = cfg["alignments_root"]
    species_cfg = cfg["species"]

    for species, sp_cfg in species_cfg.items():
        species_align_dir = os.path.join(alignments_root, species)
        if not os.path.isdir(species_align_dir):
            print(f"WARNING: alignment directory not found: {species_align_dir}", file=sys.stderr)
            continue

        for individual in os.listdir(species_align_dir):
            bam = os.path.join(species_align_dir, individual, "merged.100.bam")
            if not os.path.exists(bam):
                continue

            out_dir = os.path.join(args.results_root, species, individual)
            print(f"Processing {species}/{individual}")
            run_individual(
                species=species,
                individual=individual,
                bam=bam,
                out_dir=out_dir,
                panels_dir=args.panels_dir,
                gene_defs_dir=args.gene_defs_dir,
                amplicone_tool=args.amplicone_tool,
                sp_cfg=sp_cfg,
            )

    print("Batch run complete.")


if __name__ == "__main__":
    main()
