#!/usr/bin/env python3
"""
02_generate_gene_definitions.py
================================
Parse NCBI RefSeq GFF files to extract coordinates of single-copy Y-chromosome
control genes and append them as CONTROL entries to the species-specific gene
definition TSV/BED files in gene_definitions/targets/.

The multicopy gene family entries (TSPY, RBMY, CDY, DAZ, HSFY, BPY2, VCY) were
curated manually and must already be present in the target files before running
this script.  Only the CONTROL entries are generated automatically.

GFF paths and chrY accessions are read from config.yaml (no paths are
hardcoded in this script).

Environment: Python 3  (conda env: 02_Amp_py3)
Requires   : pyyaml, bcbio-gff  (pip install pyyaml bcbio-gff)

Usage:
    python3 02_generate_gene_definitions.py [--config /path/to/config.yaml]

Outputs appended to:
    gene_definitions/targets/{species}.tsv
    gene_definitions/targets/{species}.bed

Single-copy control genes used:
    Canonical: DDX3Y, KDM5D, NLGN4Y, TBL1Y, USP9Y, UTY, ZFY, SRY, AMELY
    Species-specific LOC aliases:
        LOC129530426  = TBL1Y  (Gorilla / GorGor)
"""

import argparse
import os
from pathlib import Path

import yaml
from BCBio import GFF

SINGLE_COPY_GENES = [
    # Canonical Y-linked single-copy genes
    "DDX3Y", "KDM5D", "NLGN4Y", "TBL1Y", "USP9Y", "UTY", "ZFY", "SRY", "AMELY",
    # Species-specific LOC aliases (manually identified from NCBI annotations)
    "LOC129530426",   # TBL1Y  — Gorilla (GorGor)
]


def load_config(config_path: str) -> dict:
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def main():
    default_config = Path(__file__).parent.parent / "config.yaml"
    default_target_dir = Path(__file__).parent.parent / "gene_definitions" / "targets"

    parser = argparse.ArgumentParser(description="Append CONTROL entries to gene definition files.")
    parser.add_argument("--config", default=str(default_config),
                        help="Path to config.yaml (default: ../config.yaml relative to this script)")
    parser.add_argument("--target_dir", default=str(default_target_dir),
                        help="Directory containing gene definition TSV/BED files to append to")
    args = parser.parse_args()

    cfg = load_config(args.config)
    species_cfg = cfg["species"]
    target_dir = args.target_dir

    os.makedirs(target_dir, exist_ok=True)

    # Track occurrences for a summary printout
    found_counts = {gene: {sp: 0 for sp in species_cfg} for gene in SINGLE_COPY_GENES}

    for species, sp_cfg in species_cfg.items():
        gff_path = sp_cfg["gff"]
        accession = sp_cfg["chry_accession"]
        print(f"Processing {species} ({gff_path}) ...")

        tsv_out = os.path.join(target_dir, f"{species}.tsv")
        bed_out = os.path.join(target_dir, f"{species}.bed")

        with open(gff_path) as gff_fh, \
             open(tsv_out, "a") as tsv_fh, \
             open(bed_out, "a") as bed_fh:

            for record in GFF.parse(gff_fh):
                for feature in record.features:
                    if feature.type != "gene":
                        continue
                    gene_name = feature.qualifiers.get("Name", [None])[0]
                    if gene_name not in SINGLE_COPY_GENES:
                        continue

                    strand = "+" if feature.location.strand == 1 else "-"
                    start  = int(feature.location.start)
                    end    = int(feature.location.end)

                    print(f"  {gene_name}  {start}-{end} ({strand})")
                    found_counts[gene_name][species] += 1

                    # TSV for AmpliCoNE-build: start, end, label
                    tsv_fh.write(f"{start}\t{end}\tCONTROL\n")
                    # BED for visualisation: chrom, start, end, name, score, strand
                    bed_fh.write(f"{accession}\t{start}\t{end}\tCONTROL_{gene_name}\t.\t{strand}\n")

    print("\n=== Control gene copy counts per species ===")
    for gene in SINGLE_COPY_GENES:
        counts = "  ".join(f"{sp}:{found_counts[gene][sp]}" for sp in species_cfg)
        print(f"  {gene:20s}  {counts}")


if __name__ == "__main__":
    main()
