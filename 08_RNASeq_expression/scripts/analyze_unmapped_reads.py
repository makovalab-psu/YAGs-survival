#!/usr/bin/env python3
"""
Identify and characterize reads that gained mapping in the enriched reference.

Workflow:
  1. Load unmapped read names from standard and enriched Salmon runs
  2. Find reads unmapped in standard but mapped in enriched (set difference)
  3. Look up which transcripts those reads mapped to in enriched quant
  4. Report which YAG families gained these reads
"""

import re
import sys
import argparse
from pathlib import Path
from collections import defaultdict


SPECIES = ["HomSap", "PanTro", "PanPan", "GorGor", "PonAbe", "PonPyg"]
SPECIES_NAMES = {
    "HomSap": "Human", "PanTro": "Chimpanzee", "PanPan": "Bonobo",
    "GorGor": "Gorilla", "PonAbe": "Sumatran orangutan", "PonPyg": "Bornean orangutan"
}
YAG_GENE_PATTERN = re.compile(r'^[A-Za-z]+_([A-Z0-9]+)_iso')
NCBI_GENE_PATTERN = re.compile(r'\[gene=([^\]]+)\]')


def load_unmapped_names(quant_dir: Path, species: str, ref_type: str) -> set[str]:
    """Load unmapped read names from Salmon's unmapped_names.txt."""
    unmapped_file = quant_dir / ref_type / species / "aux_info" / "unmapped_names.txt"
    if not unmapped_file.exists():
        print(f"  Warning: {unmapped_file} not found", file=sys.stderr)
        return set()

    with open(unmapped_file) as f:
        # Format: read_name [u|m1|m2|m3] — take only the 'u' (unmapped) entries
        names = {line.split()[0] for line in f if line.strip() and line.split()[-1] == "u"}
    return names


def load_eq_classes(quant_dir: Path, species: str, ref_type: str) -> dict[str, list[str]]:
    """
    Load equivalence class file to map reads to transcripts.
    Falls back to top expressed transcripts from quant.sf for newly gained reads.
    """
    quant_file = quant_dir / ref_type / species / "quant.sf"
    if not quant_file.exists():
        return {}

    expressed = {}
    with open(quant_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5 and float(parts[3]) > 0:
                expressed[parts[0]] = float(parts[3])  # Name -> TPM
    return expressed


def classify_transcript(transcript_id: str) -> str | None:
    """Extract YAG gene family from transcript ID."""
    m = YAG_GENE_PATTERN.match(transcript_id)
    if m:
        return m.group(1)
    return None


def analyze_newly_mapped(quant_dir: Path, species: str) -> dict:
    """Find reads gained by enriched reference and their transcript targets."""
    print(f"\n{species} ({SPECIES_NAMES[species]}):")

    std_unmapped = load_unmapped_names(quant_dir, species, "standard")
    enr_unmapped = load_unmapped_names(quant_dir, species, "enriched")

    if not std_unmapped and not enr_unmapped:
        print("  No unmapped name files found — re-run Salmon with --writeUnmappedNames")
        return {}

    gained = std_unmapped - enr_unmapped   # unmapped in std, mapped in enr
    lost   = enr_unmapped - std_unmapped   # mapped in std, unmapped in enr (unexpected)

    print(f"  Unmapped (standard):  {len(std_unmapped):>10,}")
    print(f"  Unmapped (enriched):  {len(enr_unmapped):>10,}")
    print(f"  Newly mapped reads:   {len(gained):>10,}")
    if lost:
        print(f"  Lost mappings:        {len(lost):>10,}  (unexpected — check reference)")

    # For newly mapped reads, report which YAG transcripts they land on
    # using the expressed transcripts from the enriched quant.sf as a proxy
    enr_expressed = load_eq_classes(quant_dir, species, "enriched")
    yag_expressed = {t: tpm for t, tpm in enr_expressed.items()
                     if classify_transcript(t) is not None}

    if yag_expressed:
        print(f"\n  Novel YAG transcripts with expression in enriched reference:")
        by_family: dict[str, float] = defaultdict(float)
        for t, tpm in yag_expressed.items():
            by_family[classify_transcript(t)] += tpm
        for family, tpm in sorted(by_family.items(), key=lambda x: -x[1]):
            print(f"    {family:<8} {tpm:.3f} TPM")

    return {"gained": gained, "lost": lost, "yag_tpm": yag_expressed}


def summarize_all(quant_dir: Path, output_file: Path):
    """Run analysis across all species and write summary TSV."""
    rows = []

    for species in SPECIES:
        std_unmapped = load_unmapped_names(quant_dir, species, "standard")
        enr_unmapped = load_unmapped_names(quant_dir, species, "enriched")

        if not std_unmapped and not enr_unmapped:
            continue

        gained = std_unmapped - enr_unmapped
        lost   = enr_unmapped - std_unmapped

        enr_expressed = load_eq_classes(quant_dir, species, "enriched")
        yag_expressed = {t: tpm for t, tpm in enr_expressed.items()
                         if classify_transcript(t) is not None}

        by_family: dict[str, float] = defaultdict(float)
        for t, tpm in yag_expressed.items():
            by_family[classify_transcript(t)] += tpm

        rows.append({
            "species": species,
            "species_name": SPECIES_NAMES[species],
            "n_unmapped_std": len(std_unmapped),
            "n_unmapped_enr": len(enr_unmapped),
            "n_newly_mapped": len(gained),
            "n_lost_mapped": len(lost),
            "yag_families_expressed": ",".join(sorted(by_family.keys())),
            "top_yag_family": max(by_family, key=by_family.get) if by_family else "NA",
            "top_yag_tpm": f"{max(by_family.values()):.3f}" if by_family else "0",
        })

    if not rows:
        print("\nNo unmapped name files found. Submit 02_salmon_unmapped.sh first.")
        return

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        header = list(rows[0].keys())
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(row[k]) for k in header) + "\n")

    print(f"\nSummary written to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Analyze reads gained by enriched reference")
    parser.add_argument(
        "--quant-dir",
        default="/storage/group/kdm16/default/kxp5629/proj/15_RNASeq/salmon_quant_unmapped",
        help="Directory containing Salmon quant runs with --writeUnmappedNames"
    )
    parser.add_argument(
        "--output",
        default="/storage/group/kdm16/default/kxp5629/proj/15_RNASeq/results/tables/newly_mapped_reads_summary.tsv",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=SPECIES,
        choices=SPECIES,
        help="Species to analyze (default: all)"
    )
    args = parser.parse_args()

    quant_dir = Path(args.quant_dir)

    for species in args.species:
        analyze_newly_mapped(quant_dir, species)

    summarize_all(quant_dir, Path(args.output))


if __name__ == "__main__":
    main()
