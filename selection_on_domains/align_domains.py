"""
Align gene family sequences to domain sequences via local alignment,
extract matching regions, and produce domain-level MSA with MAFFT.

Usage:
    python align_domains.py CDY
    python align_domains.py RBMY
    python align_domains.py HSFY
    python align_domains.py TSPY
"""

import sys
import os
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

GENE_DATA = Path("sequences")
DOMAINS_DIR = Path("domains")
OUTPUT_DIR = Path("output")

FAMILY_CONFIG = {
    "CDY":  ("CDY",  []),
    "RBMY": ("RBMY", []),
    "HSFY": ("HSFY", []),
    "TSPY": ("TSPY", []),
}


def get_sequence_files(family: str) -> list[Path]:
    folder_name, excludes = FAMILY_CONFIG[family]
    folder = GENE_DATA / folder_name
    return [
        f for f in folder.iterdir()
        if f.suffix in (".fasta", ".fa", ".faa")
        and f.name not in excludes
        and not f.name.startswith(".")
    ]


def translate_cds(nt_seq: str) -> tuple[str, str]:
    """Return (cleaned_nt, protein). Cleans gaps, trims to codon boundary."""
    nt_clean = nt_seq.upper().replace("-", "")
    trim = len(nt_clean) - (len(nt_clean) % 3)
    nt_clean = nt_clean[:trim]
    prot = str(Seq(nt_clean).translate(to_stop=True))
    return nt_clean, prot


def load_sequences(seq_files: list[Path]) -> list[tuple[SeqRecord, SeqRecord]]:
    """Return list of (nt_record, protein_record) pairs, deduplicated by sequence ID.
    First occurrence wins — handles folders with both per-species files and combined
    alignment files that would otherwise produce duplicate IDs.
    """
    records = []
    seen_ids: set[str] = set()
    for f in seq_files:
        for rec in SeqIO.parse(f, "fasta"):
            if rec.id in seen_ids:
                continue
            seq_str = str(rec.seq)
            non_nt = sum(1 for c in seq_str.upper() if c not in "ACGTN-")
            if non_nt / max(len(seq_str), 1) > 0.05:
                continue
            nt_clean, prot = translate_cds(seq_str)
            if len(prot) < 10:
                continue
            seen_ids.add(rec.id)
            nt_rec = SeqRecord(Seq(nt_clean), id=rec.id, description="")
            prot_rec = SeqRecord(Seq(prot.upper()), id=rec.id, description="")
            records.append((nt_rec, prot_rec))
    return records


def extract_domain_nt(
    nt_seq: str, protein: str, domain: str,
    aligner: PairwiseAligner, min_score_frac: float = 0.5
) -> str | None:
    """Local-align protein vs domain; return nucleotide substring for the aligned region.
    min_score_frac: minimum score as fraction of domain length (default 0.5).
    """
    alns = aligner.align(protein.lower(), domain.lower())
    try:
        best = next(iter(alns))
    except StopIteration:
        return None
    if best.score < min_score_frac * len(domain):
        return None
    aa_start = best.aligned[0][0][0]
    aa_end = best.aligned[0][-1][1]
    return nt_seq[aa_start * 3 : aa_end * 3]


def run_mafft(sequences: list[SeqRecord]) -> list[SeqRecord]:
    with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w", delete=False) as fh:
        SeqIO.write(sequences, fh, "fasta")
        tmp_in = fh.name
    tmp_out = tmp_in + ".aligned"
    result = subprocess.run(
        ["mafft", "--auto", "--nuc", "--quiet", tmp_in],
        capture_output=True, text=True
    )
    os.unlink(tmp_in)
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed: {result.stderr}")
    with open(tmp_out, "w") as fh:
        fh.write(result.stdout)
    aligned = list(SeqIO.parse(tmp_out, "fasta"))
    os.unlink(tmp_out)
    return aligned


def main(family: str) -> None:
    if family not in FAMILY_CONFIG:
        raise ValueError(f"Unknown family '{family}'. Choose from: {list(FAMILY_CONFIG)}")

    domain_file = DOMAINS_DIR / family
    if not domain_file.exists():
        raise FileNotFoundError(f"Domain file not found: {domain_file}")

    domains = list(SeqIO.parse(domain_file, "fasta"))
    seq_files = get_sequence_files(family)
    print(f"[{family}] Loading sequences from {len(seq_files)} file(s)...")
    seq_pairs = load_sequences(seq_files)
    print(f"[{family}] Loaded {len(seq_pairs)} nucleotide sequences")

    aligner = PairwiseAligner()
    aligner.mode = "local"

    OUTPUT_DIR.mkdir(exist_ok=True)

    for domain_rec in domains:
        domain_name = domain_rec.id.replace(" ", "_")
        domain_seq = str(domain_rec.seq)
        print(f"  Domain: {domain_name} ({len(domain_seq)} aa)")

        extracted: list[SeqRecord] = []
        for nt_rec, prot_rec in seq_pairs:
            region_nt = extract_domain_nt(
                str(nt_rec.seq), str(prot_rec.seq), domain_seq, aligner
            )
            if region_nt:
                extracted.append(SeqRecord(Seq(region_nt), id=nt_rec.id, description=""))
            else:
                print(f"    WARNING: no alignment for {nt_rec.id}")

        if not extracted:
            print(f"    SKIP: no sequences extracted for {domain_name}")
            continue

        print(f"  Aligning {len(extracted)} sequences with MAFFT...")
        aligned = run_mafft(extracted)

        subdir = OUTPUT_DIR / f"{family}_{domain_name}"
        subdir.mkdir(exist_ok=True)
        out_path = subdir / "alignment.fasta"
        SeqIO.write(aligned, out_path, "fasta")
        print(f"  Saved: {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    os.chdir(Path(__file__).parent)
    main(sys.argv[1].upper())
