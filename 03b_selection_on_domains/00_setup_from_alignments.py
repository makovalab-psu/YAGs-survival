"""
Set up output dirs for CDY_CDYL, HSFY, RBMY_RBMX from pre-built clean alignments.

For each gene:
  - Loads .clean.fasta from alignments_and_trees/
  - Maps domain protein positions to alignment columns via a reference sequence
  - Writes domain-only and nodomain alignment.fasta to output/{GENE}_{domain}/ etc.
  - Steps 01-08 then run as normal (clean_stops, dedup, trees, HyPhy)

Usage:
    python 00_setup_from_alignments.py
"""

import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

os.chdir(Path(__file__).parent)

ALN_DIR     = Path("alignments_and_trees")
DOMAINS_DIR = Path("domains")
OUTPUT_DIR  = Path("output")

# (clean_fasta_stem, domain_family, ref_seq_keyword)
GENES = [
    ("CDY_CDYL", "CDY",  "CDY1"),
    ("HSFY",     "HSFY", "HSFY1"),
    ("RBMY_RBMX","RBMY", "RBMY1"),
    ("TSPY",     "TSPY", "TSPY1"),
]


def pick_reference(records: list[SeqRecord], keyword: str) -> SeqRecord:
    kw = keyword.upper()
    for rec in records:
        if kw in rec.id.upper():
            return rec
    return records[0]


def ref_col_to_aa(ref_row: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    aa_pos = 0
    for col in range(len(ref_row) // 3):
        codon = ref_row[col * 3: col * 3 + 3]
        if "-" not in codon:
            mapping[col] = aa_pos
            aa_pos += 1
    return mapping


def find_domain_columns(
    ref_protein: str, domain_seq: str, col_to_aa: dict[int, int],
    min_score_frac: float = 0.5,
) -> set[int]:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    alns = aligner.align(ref_protein.lower(), domain_seq.lower())
    try:
        best = next(iter(alns))
    except StopIteration:
        return set()
    if best.score < min_score_frac * len(domain_seq):
        print(f"    WARNING: low score ({best.score:.1f})")
        return set()
    aa_start = best.aligned[0][0][0]
    aa_end   = best.aligned[0][-1][1]
    print(f"    ref aa {aa_start}–{aa_end} (score={best.score:.1f})")
    return {col for col, aa in col_to_aa.items() if aa_start <= aa < aa_end}


def filter_columns(records: list[SeqRecord], keep_cols: list[int]) -> list[SeqRecord]:
    result = [
        SeqRecord(Seq("".join(str(r.seq)[c * 3: c * 3 + 3] for c in keep_cols)), id=r.id, description="")
        for r in records
    ]
    # drop all-gap codon columns
    n = len(keep_cols)
    all_gap = {c for c in range(n) if all(str(r.seq)[c * 3: c * 3 + 3] == "---" for r in result)}
    if all_gap:
        keep2 = [c for c in range(n) if c not in all_gap]
        result = [
            SeqRecord(Seq("".join(str(r.seq)[c * 3: c * 3 + 3] for c in keep2)), id=r.id, description="")
            for r in result
        ]
    return result


def write_dir(gene: str, tag: str, records: list[SeqRecord]) -> None:
    d = OUTPUT_DIR / f"{gene}_{tag}"
    d.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, d / "alignment.fasta", "fasta")
    print(f"  → {d.name}: {len(records)} seqs, {len(records[0].seq) // 3} codons")


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    for gene, domain_family, ref_kw in GENES:
        fasta_path  = ALN_DIR / f"{gene}.clean.fasta"
        domain_file = DOMAINS_DIR / domain_family

        print(f"\n[{gene}]")
        if not fasta_path.exists():
            print(f"  SKIP: {fasta_path} not found")
            continue
        if not domain_file.exists():
            print(f"  SKIP: {domain_file} not found")
            continue

        records = list(SeqIO.parse(fasta_path, "fasta"))
        print(f"  Loaded {len(records)} sequences, {len(records[0].seq) // 3} codons")

        ref = pick_reference(records, ref_kw)
        print(f"  Reference: {ref.id}")
        ref_nt = str(ref.seq).upper().replace("-", "")
        ref_protein = str(Seq(ref_nt[: len(ref_nt) - len(ref_nt) % 3]).translate(to_stop=False))
        col_to_aa = ref_col_to_aa(str(ref.seq).upper())

        domains = list(SeqIO.parse(domain_file, "fasta"))
        all_domain_cols: set[int] = set()
        for dom in domains:
            dom_seq = str(dom.seq).replace("-", "")
            print(f"  Domain: {dom.id} ({len(dom_seq)} aa)")
            cols = find_domain_columns(ref_protein, dom_seq, col_to_aa)
            print(f"    {len(cols)} codon columns")
            if cols:
                n_codons = len(records[0].seq) // 3
                write_dir(gene, dom.id, filter_columns(records, sorted(cols)))
            all_domain_cols |= cols

        print(f"  Total domain columns: {len(all_domain_cols)}")
        n_codons = len(records[0].seq) // 3
        nodomain_cols = [c for c in range(n_codons) if c not in all_domain_cols]
        write_dir(gene, "nodomain", filter_columns(records, nodomain_cols))


if __name__ == "__main__":
    main()
