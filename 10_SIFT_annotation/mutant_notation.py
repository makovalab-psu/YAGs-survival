#!/usr/bin/env python3
"""For each codon position in mutated_positions.txt, report the two most common amino acids.
Prints consensus protein FASTA (MSA-length, gaps as '-') followed by mutation lines.

Usage: python mutant_notation.py <gene>

Output: FASTA block then one AA1<position>AA2 per position (e.g. A45T)
"""

import sys
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def codon_aa_counts(records: list, codon_site: int) -> Counter:
    """Return amino acid counts at 1-indexed codon_site across all sequences."""
    counts: Counter = Counter()
    for rec in records:
        codon = str(rec.seq)[(codon_site - 1) * 3: codon_site * 3]
        if '-' not in codon and len(codon) == 3:
            aa = str(Seq(codon).translate())
            if aa != '*':
                counts[aa] += 1
    return counts


def build_consensus(records: list, n_codons: int) -> str:
    """Return MSA-length consensus protein (one char per codon column, '-' at gap columns)."""
    n_seqs = len(records)
    chars = []
    for i in range(1, n_codons + 1):
        counts = codon_aa_counts(records, i)
        total_non_gap = sum(counts.values())
        if total_non_gap > n_seqs / 2 and counts:
            chars.append(counts.most_common(1)[0][0])
        else:
            chars.append('-')
    return ''.join(chars)


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <gene>")

    gene = sys.argv[1]
    gene_dir = Path(gene)

    if not gene_dir.is_dir():
        sys.exit(f"Error: directory '{gene}' not found")

    fasta_files = list(gene_dir.glob("*.clean.fasta"))
    if len(fasta_files) != 1:
        sys.exit(f"Error: expected 1 *.clean.fasta, found {len(fasta_files)}")

    pos_file = gene_dir / "mutated_positions.txt"
    if not pos_file.exists():
        sys.exit(f"Error: {pos_file} not found")

    records = list(SeqIO.parse(fasta_files[0], "fasta"))
    if not records:
        sys.exit("Error: no sequences in FASTA file")

    n_codons = len(records[0].seq) // 3

    consensus = build_consensus(records, n_codons)

    # msa_to_gapless[i] = 1-indexed position in gapless consensus for MSA codon i (0-indexed)
    # -1 if that column is a gap
    msa_to_gapless = []
    gapless_pos = 0
    for ch in consensus:
        if ch != '-':
            gapless_pos += 1
            msa_to_gapless.append(gapless_pos)
        else:
            msa_to_gapless.append(-1)

    positions = [int(line.strip()) for line in pos_file.read_text().splitlines() if line.strip()]

    def mutation_lines(pos_to_label: callable) -> list[str]:
        lines = []
        for pos in positions:
            if pos < 1 or pos > n_codons:
                sys.stderr.write(f"Warning: position {pos} out of range (max {n_codons}), skipping\n")
                continue
            counts = codon_aa_counts(records, pos)
            top2 = counts.most_common(2)
            if len(top2) < 2:
                sys.stderr.write(f"Warning: position {pos} has fewer than 2 distinct AAs, skipping\n")
                continue
            aa1, aa2 = top2[0][0], top2[1][0]
            label = pos_to_label(pos)
            if label is not None:
                lines.append(f"{aa1}{label}{aa2}")
        return lines

    # MSA-indexed output
    print(f">{gene}_consensus\n{consensus}\n")
    for line in mutation_lines(lambda pos: pos):
        print(line)

    # Gapless output
    gapless_consensus = consensus.replace('-', '')
    print(f"\n>{gene}_consensus_gapless\n{gapless_consensus}\n")
    def gapless_label(pos: int) -> int | None:
        gl = msa_to_gapless[pos - 1]
        if gl == -1:
            sys.stderr.write(f"Warning: position {pos} is a gap column, skipping in gapless output\n")
            return None
        return gl
    for line in mutation_lines(gapless_label):
        print(line)


if __name__ == "__main__":
    main()
