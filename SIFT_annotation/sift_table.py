#!/usr/bin/env python3
"""Parse sift_input and sift_result from <gene> folder, output TSV table
with mutations in gapped (MSA) coordinates.

Usage: python sift_table.py <gene>
"""

import re
import sys
from pathlib import Path

MUT_RE = re.compile(r'^([A-Z])(\d+)([A-Z])$')
SUBST_RE = re.compile(
    r"Substitution at pos (\d+) from (\w) to (\w) is predicted to "
    r"(AFFECT PROTEIN FUNCTION|be TOLERATED) with a score of ([\d.]+)\."
)
LOWCONF_RE = re.compile(r"LOW CONFIDENCE")
NOTALLOWED_RE = re.compile(r"WARNING: Original amino acid \w at position (\d+) is not allowed")


def parse_mutations(text: str) -> list[tuple[int, str, str]]:
    return [
        (int(m.group(2)), m.group(1), m.group(3))
        for line in text.splitlines()
        if (m := MUT_RE.match(line.strip()))
    ]


def parse_sift_input(path: Path) -> tuple[list, list]:
    """Return (msa_mutations, gapless_mutations) as lists of (pos, orig, alt)."""
    text = path.read_text()
    blocks = re.split(r'>[^\n]+\n[^\n]+\n', text)
    msa_muts = parse_mutations(blocks[1] if len(blocks) > 1 else '')
    gapless_muts = parse_mutations(blocks[2] if len(blocks) > 2 else '')
    return msa_muts, gapless_muts


def parse_sift_result(path: Path) -> dict[int, dict]:
    """Return {gapless_pos: {orig, alt, prediction, score, low_confidence, not_allowed}}."""
    text = path.read_text()
    not_allowed = {int(m.group(1)) for m in NOTALLOWED_RE.finditer(text)}

    results = {}
    current_pos = None
    for line in text.splitlines():
        m = SUBST_RE.search(line)
        if m:
            pos = int(m.group(1))
            prediction = m.group(4).replace('be TOLERATED', 'TOLERATED')
            current_pos = pos
            results[pos] = {
                'orig': m.group(2),
                'alt': m.group(3),
                'prediction': prediction,
                'score': m.group(5),
                'low_confidence': False,
                'not_allowed': pos in not_allowed,
            }
        elif LOWCONF_RE.search(line) and current_pos is not None:
            results[current_pos]['low_confidence'] = True

    return results


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <gene>")

    gene = sys.argv[1]
    gene_dir = Path(gene)

    if not gene_dir.is_dir():
        sys.exit(f"Error: directory '{gene}' not found")

    sift_input_file = gene_dir / 'sift_input'
    sift_result_file = gene_dir / 'sift_result'
    for f in (sift_input_file, sift_result_file):
        if not f.exists():
            sys.exit(f"Error: {f} not found")

    msa_muts, gapless_muts = parse_sift_input(sift_input_file)
    sift = parse_sift_result(sift_result_file)

    gapless_to_msa = {g[0]: m[0] for m, g in zip(msa_muts, gapless_muts)}

    print("MSA_pos\tOrig_AA\tChanged_AA\tSIFT_prediction\tScore\tConfidence")
    for gpos, entry in sorted(sift.items()):
        msa_pos = gapless_to_msa.get(gpos, '?')
        confidence = 'Low' if entry['low_confidence'] else 'High'
        if entry['not_allowed']:
            confidence += '*'
        print(f"{msa_pos}\t{entry['orig']}\t{entry['alt']}\t{entry['prediction']}\t{entry['score']}\t{confidence}")

    sys.stderr.write("# Confidence: Low = insufficient sequence diversity; * = original AA not allowed by SIFT prediction\n")


if __name__ == "__main__":
    main()
