"""Shared configuration for the selection-on-domains pipeline."""

from pathlib import Path

# Tool paths
IQTREE   = "/Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2"
HYPHY    = "/usr/local/bin/hyphy"
FITMG94  = "/Users/kxp5629/proj/Y/src/29_recreate_selection/FitMG94.bf"

OUTPUT_DIR = Path("output")


def get_family(subdir: Path) -> str:
    """Derive gene family name from output subdirectory name.

    Uses the FAMILY_PREFIXES map for imported alignments with compound names,
    falls back to the first underscore-delimited token for domain-level outputs.
    """
    name = subdir.name
    for prefix in FAMILY_PREFIXES:
        if name.startswith(prefix):
            return prefix
    return name.split("_")[0]


# Compound family prefixes for imported alignments (longest-match first)
FAMILY_PREFIXES = [
    "DAZ-DAZL",
    "RBMY-RBMX",
    "CDY_CDYL",
    "CDY",
    "RBMY",
    "HSFY",
    "TSPY",
    "DAZ",
]


def is_outgroup(seq_id: str, family: str) -> bool:
    """Return True if seq_id should be used as an outgroup for the given family."""
    if family in ("CDY", "CDY_CDYL"):
        return "CDYL" in seq_id and "CDYL2" not in seq_id
    if family in ("RBMY", "RBMY-RBMX", "RBMY_RBMX"):
        return "RBMX" in seq_id
    if family in ("HSFY", "TSPY"):
        return seq_id.startswith("SymSyn")
    if family in ("DAZ", "DAZ-DAZL"):
        return "DAZL" in seq_id
    return False


ANALYSIS_ALIGNMENT = "alignment.dedup.fasta"


def iter_domains() -> list[Path]:
    """Return all output subdirectories (one per family+domain combination)."""
    return sorted(p for p in OUTPUT_DIR.iterdir() if p.is_dir())
