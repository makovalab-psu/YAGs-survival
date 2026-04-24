"""
Remove domain columns from full codon alignments.

For each gene family:
  1. Reads output_full/{GENE}/protein_alignment.fasta to map alignment
     columns → reference protein positions (via __REFERENCE__ row).
  2. Local-aligns each domain protein (from domains/{GENE}) against the
     unaligned reference protein to find domain boundaries.
  3. Removes those columns from the codon alignment.
  4. Writes output/{GENE}_nodomain/alignment.fasta so the existing
     selection pipeline (steps 01-08) can process it directly.

Usage:
    python mask_domains.py CDY     # one family
    python mask_domains.py         # all families
"""

import sys
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

os.chdir(Path(__file__).parent)

FULL_OUT_DIR = Path("output_full")
DOMAINS_DIR  = Path("domains")
OUTPUT_DIR   = Path("output")
REFERENCE_ID = "__REFERENCE__"
FULL_SEQ_DIR = Path("full_sequence")

FAMILIES = ["CDY", "HSFY", "RBMY", "TSPY"]


def load_reference_protein(family: str) -> str:
    text = (FULL_SEQ_DIR / family).read_text()
    return "".join(c for c in text.upper() if c.isalpha())


def build_col_to_refpos(ref_aligned: str) -> dict[int, int]:
    """Map alignment column index → reference protein position (0-based)."""
    mapping: dict[int, int] = {}
    ref_pos = 0
    for col, aa in enumerate(ref_aligned):
        if aa != "-":
            mapping[col] = ref_pos
            ref_pos += 1
    return mapping


def find_domain_columns(
    ref_protein: str, domain_seq: str,
    col_to_refpos: dict[int, int],
    min_score_frac: float = 0.5,
) -> set[int]:
    """Return alignment columns that fall within the domain region."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    alns = aligner.align(ref_protein.lower(), domain_seq.lower())
    try:
        best = next(iter(alns))
    except StopIteration:
        return set()
    if best.score < min_score_frac * len(domain_seq):
        print(f"    WARNING: low alignment score ({best.score:.1f} < {min_score_frac * len(domain_seq):.1f})")
        return set()
    aa_start = best.aligned[0][0][0]
    aa_end   = best.aligned[0][-1][1]
    print(f"    domain maps to ref positions {aa_start}–{aa_end} (score={best.score:.1f})")
    return {col for col, pos in col_to_refpos.items() if aa_start <= pos < aa_end}


def remove_columns(records: list[SeqRecord], cols_to_remove: set[int]) -> list[SeqRecord]:
    """Remove specified columns (by codon triplet) and all-gap columns."""
    if not records:
        return records
    aln_len = len(records[0].seq)
    n_codons = aln_len // 3

    # Convert protein-space columns to nucleotide-space codon triplets to keep
    keep_codons = [
        c for c in range(n_codons)
        if c not in cols_to_remove
    ]

    result = []
    for rec in records:
        seq = str(rec.seq)
        new_seq = "".join(seq[c * 3: c * 3 + 3] for c in keep_codons)
        result.append(SeqRecord(Seq(new_seq), id=rec.id, description=""))

    # Remove all-gap codon columns
    n_codons_new = len(keep_codons)
    all_gap = {
        c for c in range(n_codons_new)
        if all(str(r.seq)[c * 3: c * 3 + 3] == "---" for r in result)
    }
    if all_gap:
        print(f"    removing {len(all_gap)} all-gap codon columns after masking")
        keep2 = [c for c in range(n_codons_new) if c not in all_gap]
        result = [
            SeqRecord(
                Seq("".join(str(r.seq)[c * 3: c * 3 + 3] for c in keep2)),
                id=r.id, description=""
            )
            for r in result
        ]

    return result


def main(families: list[str]) -> None:
    for family in families:
        full_dir  = FULL_OUT_DIR / family
        prot_path = full_dir / "protein_alignment.fasta"
        aln_path  = full_dir / "alignment.fasta"
        domain_file = DOMAINS_DIR / family

        if not prot_path.exists() or not aln_path.exists():
            print(f"[{family}] SKIP: run align_full.py first (missing {full_dir}/)")
            continue
        if not domain_file.exists():
            print(f"[{family}] SKIP: no domain file at {domain_file}")
            continue

        ref_protein = load_reference_protein(family)

        # Build column → reference protein position mapping
        prot_aln = {r.id: str(r.seq) for r in SeqIO.parse(prot_path, "fasta")}
        if REFERENCE_ID not in prot_aln:
            print(f"[{family}] SKIP: {REFERENCE_ID} not found in protein_alignment.fasta")
            continue
        col_to_refpos = build_col_to_refpos(prot_aln[REFERENCE_ID])
        print(f"[{family}] Reference spans {len(col_to_refpos)} alignment columns")

        # Find all domain columns
        domains = list(SeqIO.parse(domain_file, "fasta"))
        domain_cols: set[int] = set()
        for dom in domains:
            dom_seq = str(dom.seq)
            print(f"  Domain: {dom.id} ({len(dom_seq)} aa)")
            cols = find_domain_columns(ref_protein, dom_seq, col_to_refpos)
            print(f"    → {len(cols)} alignment columns masked")
            domain_cols |= cols

        print(f"  Total masked: {len(domain_cols)} codon columns")

        codon_recs = list(SeqIO.parse(aln_path, "fasta"))
        n_before = len(str(codon_recs[0].seq)) // 3 if codon_recs else 0

        masked = remove_columns(codon_recs, domain_cols)
        n_after = len(str(masked[0].seq)) // 3 if masked else 0
        print(f"  Alignment: {n_before} → {n_after} codon positions ({n_before - n_after} removed)")

        out_dir = OUTPUT_DIR / f"{family}_nodomain"
        out_dir.mkdir(exist_ok=True)
        out_path = out_dir / "alignment.fasta"
        SeqIO.write(masked, out_path, "fasta")
        print(f"  Saved → {out_path}")


if __name__ == "__main__":
    all_families = FAMILIES
    families = [f.upper() for f in sys.argv[1:]] if len(sys.argv) > 1 else all_families
    main(families)
