"""
Step 02 — Remove exact duplicate sequences using HyPhy remove-duplicates.bf.
Input:  output/{family}_{domain}/alignment.clean.fasta
Output: output/{family}_{domain}/alignment.dedup.fasta
"""

import os
import re
import subprocess
import tempfile
from pathlib import Path

os.chdir(Path(__file__).parent)
from pipeline_config import iter_domains, HYPHY

REMOVE_DUPLICATES_BF = "/Users/kxp5629/proj/Y/hyphy-analyses/remove-duplicates/remove-duplicates.bf"


def parse_hyphy_nexus(nexus_path: Path) -> list[tuple[str, str]]:
    """Parse HyPhy NEXUS (NOLABELS) → list of (name, sequence)."""
    text = nexus_path.read_text()
    # Extract taxon names
    taxa_block = re.search(r"TAXLABELS\s+(.*?);", text, re.DOTALL)
    if not taxa_block:
        raise ValueError("No TAXLABELS block found")
    names = re.findall(r"'([^']+)'", taxa_block.group(1))

    # Extract matrix sequences
    matrix_block = re.search(r"MATRIX\s*(.*?)\s*;", text, re.DOTALL)
    if not matrix_block:
        raise ValueError("No MATRIX block found")
    seqs = [line.strip() for line in matrix_block.group(1).splitlines() if line.strip()]

    if len(names) != len(seqs):
        raise ValueError(f"TAXLABELS ({len(names)}) vs MATRIX rows ({len(seqs)}) mismatch")
    return list(zip(names, seqs))


def nexus_to_fasta(nexus_path: Path, fasta_path: Path) -> None:
    pairs = parse_hyphy_nexus(nexus_path)
    with open(fasta_path, "w") as fh:
        for name, seq in pairs:
            fh.write(f">{name}\n{seq}\n")


def run_dedup(input_path: Path, output_path: Path) -> None:
    with tempfile.NamedTemporaryFile(suffix=".nex", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    cmd = [
        HYPHY, REMOVE_DUPLICATES_BF,
        "--msa", str(input_path.resolve()),
        "--output", str(tmp_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        tmp_path.unlink(missing_ok=True)
        print(f"  ERROR: {result.stderr}\n{result.stdout}")
        raise RuntimeError(f"HyPhy dedup failed for {input_path}")

    nexus_to_fasta(tmp_path, output_path)
    tmp_path.unlink(missing_ok=True)


def main() -> None:
    for d in iter_domains():
        inp = d / "alignment.clean.fasta"
        out = d / "alignment.dedup.fasta"
        if not inp.exists():
            print(f"  SKIP {d.name}: alignment.clean.fasta not found")
            continue
        run_dedup(inp, out)
        print(f"  {d.name}: done → {out.name}")


if __name__ == "__main__":
    main()
