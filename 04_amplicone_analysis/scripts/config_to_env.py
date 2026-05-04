#!/usr/bin/env python3
"""
config_to_env.py
================
Read config.yaml and print bash variable assignments that the pipeline
shell script can eval-source.

Usage (inside amplicone_pipeline.sh):
    eval "$(python3 "${SCRIPTS_DIR}/config_to_env.py" "${PIPELINE_DIR}/config.yaml")"

Exported variables:
    TRF_BIN, GEM_BIN_DIR, BEDOPS_BIN_DIR
    ALIGNMENTS_ROOT
    SPECIES                          — space-separated list of active species
    {SPECIES}_REFERENCE              — e.g. GORGOR_REFERENCE
    {SPECIES}_GFF
    {SPECIES}_REPEATMASKER_OUT
    {SPECIES}_CHRY_RM_TAG
    {SPECIES}_CHRY_ACCESSION
    {SPECIES}_CHRY_LENGTH
    {SPECIES}_GENE_DEFINITIONS
    {SPECIES}_ANNOTATION_TAB
"""

import sys
import yaml   # pip install pyyaml


def bash_var(species: str, field: str) -> str:
    """Return the bash variable name for a species field, e.g. GORGOR_REFERENCE."""
    return f"{species.upper()}_{field.upper()}"


def quote(value) -> str:
    return f'"{value}"'


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <config.yaml>", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as fh:
        cfg = yaml.safe_load(fh)

    lines = []

    # Tools
    tools = cfg.get("tools", {})
    lines.append(f'TRF_BIN={quote(tools["trf_bin"])}')
    lines.append(f'GEM_BIN_DIR={quote(tools["gem_bin_dir"])}')
    lines.append(f'BEDOPS_BIN_DIR={quote(tools["bedops_bin_dir"])}')

    # Alignments root
    lines.append(f'ALIGNMENTS_ROOT={quote(cfg["alignments_root"])}')

    # Per-species variables
    species_cfg = cfg.get("species", {})
    lines.append(f'SPECIES={quote(" ".join(species_cfg.keys()))}')

    for sp, vals in species_cfg.items():
        for field in ("reference", "gff", "repeatmasker_out", "chry_rm_tag",
                      "chry_accession", "chry_length", "gene_definitions", "annotation_tab"):
            var = bash_var(sp, field)
            lines.append(f'{var}={quote(vals[field])}')

    print("\n".join(lines))


if __name__ == "__main__":
    main()
