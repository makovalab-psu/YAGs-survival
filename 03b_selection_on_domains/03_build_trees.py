"""
Step 03 — Build rooted ML trees with IQtree2.
Outgroups are detected automatically from sequence IDs based on gene family:
  CDY  → CDYL sequences (not CDYL2)
  RBMY → RBMX sequences
  HSFY → SymSyn sequences
  TSPY → SymSyn sequences

Input:  output/{family}_{domain}/{ANALYSIS_ALIGNMENT}
Output: output/{family}_{domain}/tree.treefile  (+ IQtree ancillary files)
"""

import os
import subprocess
from pathlib import Path
from Bio import SeqIO

os.chdir(Path(__file__).parent)
from pipeline_config import IQTREE, ANALYSIS_ALIGNMENT, iter_domains, get_family, is_outgroup


def get_outgroups(alignment: Path, family: str) -> list[str]:
    # Names are already sanitized in alignment.clean.fasta (step 01)
    return [
        rec.id for rec in SeqIO.parse(alignment, "fasta")
        if is_outgroup(rec.id, family)
    ]


def run_iqtree(alignment: Path, outgroups: list[str], prefix: Path) -> bool:
    cmd = [
        IQTREE, "-s", str(alignment),
        "-bb", "1000",
        "--prefix", str(prefix),
        "-keep-ident",
        "--redo",
    ]
    if outgroups:
        cmd += ["-o", ",".join(outgroups)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    IQtree FAILED:\n{result.stderr[-500:]}")
        return False
    return True


def main() -> None:
    for d in iter_domains():
        inp = d / ANALYSIS_ALIGNMENT
        treefile = d / "tree.treefile"
        if not inp.exists():
            print(f"  SKIP {d.name}: {ANALYSIS_ALIGNMENT} not found")
            continue
        if treefile.exists():
            print(f"  SKIP {d.name}: tree already exists")
            continue

        family = get_family(d)
        outgroups = get_outgroups(inp, family)
        if not outgroups:
            print(f"  WARNING {d.name}: no outgroup sequences found — tree will be unrooted")
        else:
            print(f"  {d.name}: {len(outgroups)} outgroup(s)")

        prefix = d / "tree"
        ok = run_iqtree(inp, outgroups, prefix)
        if ok:
            print(f"  {d.name}: tree saved → {treefile.name}")
        else:
            print(f"  {d.name}: FAILED")


if __name__ == "__main__":
    main()
