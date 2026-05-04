#!/usr/bin/env python3
"""
04_collect_results.py
=====================
Aggregate per-individual AmpliCoNE results into a single summary TSV.

Reads OutputAmpliconic_Summary.txt from each individual's results directory:
  {results_root}/{species}/{individual}/OutputAmpliconic_Summary.txt

Each file has tab-separated rows:  GENE_FAMILY  <copy_number>

Writes one TSV with columns:
  species  individual  TSPY  VCY  CDY  HSFY  RBMY  BPY2  DAZ

Missing gene values are left as empty cells.
The species list is read from config.yaml (all entries under 'species:').

Usage:
    python3 04_collect_results.py \\
        --config       /path/to/config.yaml \\
        --results_root /path/to/results \\
        --output_dir   /path/to/output
"""

import argparse
import datetime
import os
from pathlib import Path

import yaml

GENES = ["TSPY", "VCY", "CDY", "HSFY", "RBMY", "BPY2", "DAZ"]


def load_config(config_path: str) -> dict:
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def collect(results_root: str, species_list: list) -> dict:
    data = {sp: {} for sp in species_list}

    for species in species_list:
        sp_dir = os.path.join(results_root, species)
        if not os.path.isdir(sp_dir):
            print(f"WARNING: results directory not found: {sp_dir}")
            continue

        for individual in os.listdir(sp_dir):
            summary_file = os.path.join(sp_dir, individual, "OutputAmpliconic_Summary.txt")
            data[species][individual] = {}

            if not os.path.exists(summary_file):
                print(f"  Missing: {summary_file}")
                continue

            with open(summary_file) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if parts[0] in GENES:
                        data[species][individual][parts[0]] = parts[1]

    return data


def write_tsv(data: dict, species_list: list, output_dir: str) -> None:
    date = datetime.date.today().strftime("%Y%m%d")
    out_path = os.path.join(output_dir, f"ampl_results.{date}.tsv")

    with open(out_path, "w") as fh:
        fh.write("species\tindividual\t" + "\t".join(GENES) + "\n")
        for species in species_list:
            for individual, gene_vals in data[species].items():
                row = f"{species}\t{individual}"
                for gene in GENES:
                    row += "\t" + gene_vals.get(gene, "")
                fh.write(row + "\n")

    print(f"Results written to: {out_path}")


def main():
    default_config = Path(__file__).parent.parent / "config.yaml"

    parser = argparse.ArgumentParser(description="Collect AmpliCoNE results into summary TSV.")
    parser.add_argument("--config", default=str(default_config),
                        help="Path to config.yaml (default: ../config.yaml relative to this script)")
    parser.add_argument("--results_root", required=True,
                        help="Root directory with {species}/{individual}/OutputAmpliconic_Summary.txt")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to write ampl_results.YYYYMMDD.tsv")
    args = parser.parse_args()

    cfg = load_config(args.config)
    species_list = list(cfg["species"].keys())

    data = collect(args.results_root, species_list)
    write_tsv(data, species_list, args.output_dir)


if __name__ == "__main__":
    main()
