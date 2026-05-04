"""
Step 01 — Replace in-frame stop codons with gap triplets (---) and sanitize
sequence names to match IQtree's internal renaming (replaces ()[],:;' with _).
Both operations must use the same names so the alignment and treefile agree.

Input:  output/{family}_{domain}/alignment.fasta
Output: output/{family}_{domain}/alignment.clean.fasta
"""

import os
from pathlib import Path

os.chdir(Path(__file__).parent)
from pipeline_config import iter_domains

STOPS = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}


def sanitize_name(seq_id: str) -> str:
    """Mimic IQtree's sequence name sanitization: replace ()[],:;' with _."""
    for ch in "()[],:;' ":
        seq_id = seq_id.replace(ch, "_")
    return seq_id


def clean_stops(input_path: Path, output_path: Path) -> int:
    """Replace stop codons with --- and sanitize sequence names.
    Returns number of stops replaced."""
    replaced = 0
    with open(input_path) as fin, open(output_path, "w") as fout:
        header, seq_id, seq_parts = None, None, []

        def flush():
            nonlocal replaced
            seq = "".join(seq_parts)
            codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
            cleaned = []
            for c in codons:
                if len(c) == 3 and c in STOPS:
                    cleaned.append("---")
                    replaced += 1
                else:
                    cleaned.append(c)
            fout.write(f">{sanitize_name(seq_id)}\n{''.join(cleaned)}\n")

        for line in fin:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    flush()
                header = line
                seq_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            flush()
    return replaced


def main() -> None:
    domains = iter_domains()
    for d in domains:
        inp = d / "alignment.fasta"
        out = d / "alignment.clean.fasta"
        if not inp.exists():
            print(f"  SKIP {d.name}: alignment.fasta not found")
            continue
        n = clean_stops(inp, out)
        print(f"  {d.name}: {n} stop codons replaced → {out.name}")


if __name__ == "__main__":
    main()
