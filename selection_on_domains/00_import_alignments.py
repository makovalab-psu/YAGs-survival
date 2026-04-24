"""
Step 00b — Import pre-built alignments into the pipeline output structure.
Copies each fasta from test_alignments/ into output/{name}/alignment.fasta.

The subdirectory name is derived from the filename stem (cleaned up),
or can be overridden with a name→subdir mapping below.

Usage:
    python 00_import_alignments.py
"""

import os
import shutil
from pathlib import Path

os.chdir(Path(__file__).parent)
from pipeline_config import OUTPUT_DIR

SOURCE_DIR = Path("test_alignments")

# filename stem → output subdirectory name
# Add entries here to override the auto-derived name.
NAME_MAP = {
    "DAZ_DAZL_RBD.aligmnent.20240924":           "DAZ-DAZL_RBD",
    "RBMY_with_merged_chry_chrx_v6.output.renamed": "RBMY-RBMX_RRM",
}


def derive_name(stem: str) -> str:
    return NAME_MAP.get(stem, stem)


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    fastas = [f for f in SOURCE_DIR.iterdir() if f.suffix in (".fasta", ".fa", ".fna")]
    if not fastas:
        print("No fasta files found in test_alignments/")
        return

    for src in sorted(fastas):
        name = derive_name(src.stem)
        dest_dir = OUTPUT_DIR / name
        dest_dir.mkdir(exist_ok=True)
        dest = dest_dir / "alignment.fasta"
        shutil.copy2(src, dest)
        print(f"  {src.name} → output/{name}/alignment.fasta")


if __name__ == "__main__":
    main()
