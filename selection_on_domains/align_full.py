"""
Protein-guided codon alignment of full gene family CDS sequences.

For each family, translates all CDS to protein, aligns proteins with MAFFT
(reference protein from full_sequence/ included as anchor), and back-translates
to produce a codon-aware nucleotide MSA.

Outputs per family in output_full/{GENE}/:
  alignment.fasta           — codon MSA (no reference row)
  protein_alignment.fasta   — protein MSA including __REFERENCE__ row
                              (used by mask_domains.py for column mapping)

Usage:
    python align_full.py CDY     # one family
    python align_full.py         # all families
"""

import sys
import os
import io
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

os.chdir(Path(__file__).parent)

GENE_DATA    = Path("sequences")
FULL_SEQ_DIR = Path("full_sequence")
OUTPUT_DIR   = Path("output_full")
REFERENCE_ID = "__REFERENCE__"

FAMILY_CONFIG = {
    "CDY":  ("CDY",  []),
    "RBMY": ("RBMY", []),
    "HSFY": ("HSFY", []),
    "TSPY": ("TSPY", []),
}

# Sequence IDs to exclude per family (poorly aligning / misannotated)
SEQ_EXCLUDES: dict[str, set[str]] = {
    "RBMY": {"PanTro_chrY_LOC100615847_-_100615847"},
}

# Exon sequences to remove per family before alignment (e.g. tandem repeat exons).
# All occurrences are removed iteratively via local alignment.
# Reading frame is preserved when the exon length is a multiple of 3 and the
# splice sites are at consistent phase (as with RBMY C-terminal repeat exons).
EXON_MASKS: dict[str, list[str]] = {
}


# Anchor-based region cuts per family: remove everything between the end of
# left_anchor and the start of right_anchor in each translated protein (and the
# corresponding nucleotides). More robust than exon-by-exon masking for
# highly variable tandem repeat regions.
# RBMY: the C-terminal SRGY repeat array sits between GYATNDG and YHDGYGE.
ANCHOR_CUTS: dict[str, tuple[str, str]] = {
    "RBMY": ("GYATNDG", "YHDGYGE"),
}

# Trim alignment start to begin at first occurrence of motif in consensus.
START_ANCHORS: dict[str, str] = {
    "CDY": "VESIVDKR",
}

# Per-family 3' tail trims: list of (taxon_prefix, nt_suffix_to_remove).
# Applied to raw CDS before translation; first matching suffix wins per sequence.
TAIL_TRIMS: dict[str, list[tuple[str, str]]] = {
    "TSPY": [
        ("GorGor", "TCCCAGATGT"),  # 6 gorilla sequences with spurious 3' repeat codons
        ("GorGor", "TTCC"),        # shorter gorilla sequence (LOC129530254_trim)
    ],
}


def consensus_protein(aligned_prots: list[SeqRecord], exclude_id: str = REFERENCE_ID) -> str:
    """Majority-vote consensus (gaps ignored) from protein alignment, excluding reference."""
    from collections import Counter
    seqs = [str(r.seq).upper() for r in aligned_prots if r.id != exclude_id]
    result = []
    for col in zip(*seqs):
        non_gap = [c for c in col if c != "-"]
        result.append(Counter(non_gap).most_common(1)[0][0] if non_gap else "-")
    return "".join(result)


def cut_alignment_by_anchors(
    codon_recs: list[SeqRecord],
    consensus: str,
    left_anchor: str,
    right_anchor: str,
    aligner: PairwiseAligner,
    min_score_frac: float = 0.85,
) -> tuple[list[SeqRecord], bool]:
    """Remove alignment columns between end-of-left_anchor and start-of-right_anchor
    as located in the alignment consensus. Returns (trimmed_records, success).
    """
    # Map alignment columns to consensus positions (strip gaps)
    col_to_pos: list[int] = []
    pos = 0
    for aa in consensus:
        if aa != "-":
            col_to_pos.append(pos)
            pos += 1
        else:
            col_to_pos.append(-1)
    consensus_nogap = consensus.replace("-", "")

    left_aln = next(iter(aligner.align(consensus_nogap.lower(), left_anchor.lower())), None)
    if left_aln is None or left_aln.score < min_score_frac * len(left_anchor):
        print(f"    WARNING: left anchor '{left_anchor}' not found in consensus")
        return codon_recs, False
    left_end_pos = left_aln.aligned[0][-1][1]  # position in gap-free consensus

    right_search = consensus_nogap[left_end_pos:]
    if not right_search:
        print(f"    WARNING: nothing after left anchor")
        return codon_recs, False
    right_aln = next(iter(aligner.align(right_search.lower(), right_anchor.lower())), None)
    if right_aln is None or right_aln.score < min_score_frac * len(right_anchor):
        print(f"    WARNING: right anchor '{right_anchor}' not found in consensus")
        return codon_recs, False
    right_start_pos = left_end_pos + right_aln.aligned[0][0][0]

    # Find alignment columns corresponding to these positions
    left_end_col  = next(c for c, p in enumerate(col_to_pos) if p == left_end_pos - 1) + 1
    right_start_col = next(c for c, p in enumerate(col_to_pos) if p == right_start_pos)

    print(f"    Anchor cut: columns {left_end_col}–{right_start_col} "
          f"({right_start_col - left_end_col} protein cols = "
          f"{(right_start_col - left_end_col) * 3} nt removed)")

    keep_codons = (
        list(range(left_end_col)) +
        list(range(right_start_col, len(consensus)))
    )
    trimmed = [
        SeqRecord(
            Seq("".join(str(r.seq)[c * 3: c * 3 + 3] for c in keep_codons)),
            id=r.id, description=""
        )
        for r in codon_recs
    ]
    return trimmed, True


def load_reference_protein(family: str) -> str:
    """Parse plain-text protein from full_sequence/{family}. Strips spaces/newlines/=."""
    text = (FULL_SEQ_DIR / family).read_text()
    return "".join(c for c in text.upper() if c.isalpha())


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
    """Return (cleaned_nt, protein). Strips gaps, trims to codon boundary."""
    nt_clean = nt_seq.upper().replace("-", "")
    trim = len(nt_clean) - (len(nt_clean) % 3)
    nt_clean = nt_clean[:trim]
    prot = str(Seq(nt_clean).translate(to_stop=True))
    return nt_clean, prot


def remove_exon_occurrences(nt_seq: str, exon_queries: list[str], min_score_frac: float = 0.85) -> tuple[str, int]:
    """Iteratively remove all local-alignment matches of each exon query from nt_seq.
    Returns (cleaned_sequence, total_n_removed).
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.gap_score = -2

    n_removed = 0
    for query in exon_queries:
        threshold = min_score_frac * len(query)
        while True:
            try:
                best = next(iter(aligner.align(nt_seq, query)))
            except StopIteration:
                break
            if best.score < threshold:
                break
            start = best.aligned[0][0][0]
            end   = best.aligned[0][-1][1]
            nt_seq = nt_seq[:start] + nt_seq[end:]
            n_removed += 1
    return nt_seq, n_removed


def load_cds_pairs(
    seq_files: list[Path],
    exclude_ids: set[str] | None = None,
    exon_masks: list[str] | None = None,
    tail_trims: list[tuple[str, str]] | None = None,
) -> list[tuple[SeqRecord, SeqRecord]]:
    """Return (nt_record, protein_record) pairs, deduplicated by ID. First occurrence wins."""
    records = []
    seen_ids: set[str] = set()
    total_exons_removed = 0

    for f in seq_files:
        for rec in SeqIO.parse(f, "fasta"):
            if rec.id in seen_ids:
                continue
            if exclude_ids and rec.id in exclude_ids:
                continue
            seq_str = str(rec.seq).upper().replace("-", "")
            non_nt = sum(1 for c in seq_str if c not in "ACGTN")
            if non_nt / max(len(seq_str), 1) > 0.05:
                continue
            if tail_trims:
                for taxon_prefix, suffix in tail_trims:
                    if rec.id.startswith(taxon_prefix):
                        # Allow a terminal stop codon (3 nt) after the suffix
                        for stop_len in (3, 0):
                            tail = seq_str[len(seq_str) - len(suffix) - stop_len:]
                            if tail.startswith(suffix):
                                seq_str = seq_str[: len(seq_str) - len(suffix) - stop_len] + seq_str[len(seq_str) - stop_len:]
                                break
                        else:
                            continue
                        break
            if exon_masks:
                seq_str, n = remove_exon_occurrences(seq_str, exon_masks)
                total_exons_removed += n
            nt_clean, prot = translate_cds(seq_str)
            if len(prot) < 10:
                continue
            seen_ids.add(rec.id)
            nt_rec   = SeqRecord(Seq(nt_clean), id=rec.id, description="")
            prot_rec = SeqRecord(Seq(prot.upper()), id=rec.id, description="")
            records.append((nt_rec, prot_rec))

    if exon_masks and total_exons_removed:
        print(f"  Exon masking: {total_exons_removed} repeat exon(s) removed across all sequences")
    return records


def mafft_protein(sequences: list[SeqRecord]) -> list[SeqRecord]:
    """Run MAFFT amino acid alignment; return aligned records."""
    with tempfile.NamedTemporaryFile(suffix=".faa", mode="w", delete=False) as fh:
        SeqIO.write(sequences, fh, "fasta")
        tmp_in = fh.name
    result = subprocess.run(
        ["mafft", "--amino", "--globalpair", "--maxiterate", "1000", "--quiet", tmp_in],
        capture_output=True, text=True
    )
    os.unlink(tmp_in)
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed: {result.stderr}")
    return list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))


def trim_to_coverage(records: list[SeqRecord], min_frac: float = 0.70) -> tuple[list[SeqRecord], int, int]:
    """Trim alignment to codon columns where >= min_frac of sequences are non-gap.
    Returns (trimmed_records, start_codon, end_codon).
    """
    n = len(records)
    seq_strs = [str(r.seq) for r in records]
    n_codons = len(seq_strs[0]) // 3

    def coverage(c: int) -> float:
        return sum(1 for s in seq_strs if s[c * 3: c * 3 + 3] != "---") / n

    # First codon with >= min_frac coverage
    start = next((c for c in range(n_codons) if coverage(c) >= min_frac), 0)
    # Last codon with >= min_frac coverage
    end = next((c for c in range(n_codons - 1, -1, -1) if coverage(c) >= min_frac), n_codons - 1) + 1

    trimmed = [
        SeqRecord(Seq(str(r.seq)[start * 3: end * 3]), id=r.id, description="")
        for r in records
    ]
    return trimmed, start, end


def back_translate(aligned_protein: str, nt_seq: str) -> str:
    """Thread nucleotides through protein alignment gaps."""
    result = []
    nt_pos = 0
    for aa in aligned_protein:
        if aa == "-":
            result.append("---")
        else:
            result.append(nt_seq[nt_pos:nt_pos + 3])
            nt_pos += 3
    return "".join(result)


def main(families: list[str]) -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    for family in families:
        ref_protein = load_reference_protein(family)
        print(f"[{family}] Reference protein: {len(ref_protein)} aa")

        seq_files = get_sequence_files(family)
        print(f"[{family}] Loading from {len(seq_files)} file(s)...")
        pairs = load_cds_pairs(seq_files, exclude_ids=SEQ_EXCLUDES.get(family),
                               exon_masks=EXON_MASKS.get(family),
                               tail_trims=TAIL_TRIMS.get(family))
        print(f"[{family}] {len(pairs)} sequences loaded")

        ref_rec = SeqRecord(Seq(ref_protein), id=REFERENCE_ID, description="")
        prot_recs = [ref_rec] + [p for _, p in pairs]

        print(f"[{family}] Aligning {len(prot_recs)} proteins with MAFFT...")
        aligned_prots = mafft_protein(prot_recs)

        prot_by_id = {r.id: str(r.seq) for r in aligned_prots}
        nt_by_id   = {nt.id: str(nt.seq) for nt, _ in pairs}

        # Codon back-translation (exclude reference)
        codon_recs = []
        for nt_rec, _ in pairs:
            seq_id = nt_rec.id
            codon_seq = back_translate(prot_by_id[seq_id], nt_by_id[seq_id])
            codon_recs.append(SeqRecord(Seq(codon_seq), id=seq_id, description=""))

        # Anchor-based repeat region cut (applied on consensus after alignment)
        anchor_cut = ANCHOR_CUTS.get(family)
        if anchor_cut:
            cons = consensus_protein(aligned_prots)
            anchor_aligner = PairwiseAligner()
            anchor_aligner.mode = "local"
            codon_recs, _ = cut_alignment_by_anchors(
                codon_recs, cons, anchor_cut[0], anchor_cut[1], anchor_aligner
            )

        # Remove all-gap codon columns (reference-only regions)
        n_codons = len(codon_recs[0].seq) // 3
        keep = [
            c for c in range(n_codons)
            if not all(str(r.seq)[c * 3: c * 3 + 3] == "---" for r in codon_recs)
        ]
        n_removed = n_codons - len(keep)
        if n_removed:
            codon_recs = [
                SeqRecord(
                    Seq("".join(str(r.seq)[c * 3: c * 3 + 3] for c in keep)),
                    id=r.id, description=""
                )
                for r in codon_recs
            ]
            print(f"[{family}] Removed {n_removed} all-gap codon columns (reference-only regions)")

        # Trim to 70% coverage start/end
        n_before_trim = len(codon_recs[0].seq) // 3
        codon_recs, trim_start, trim_end = trim_to_coverage(codon_recs, min_frac=0.70)
        n_after_trim = len(codon_recs[0].seq) // 3
        print(f"[{family}] Coverage trim: codons {trim_start}–{trim_end} "
              f"({n_before_trim - n_after_trim} removed, {n_after_trim} remain)")

        # Start-anchor trim: find motif in consensus of trimmed alignment, drop before it
        start_motif = START_ANCHORS.get(family)
        if start_motif:
            from collections import Counter as _Counter
            seqs_str = [str(r.seq) for r in codon_recs]
            n_codons_cur = len(seqs_str[0]) // 3
            trimmed_cons = ""
            for c in range(n_codons_cur):
                codons = [s[c * 3: c * 3 + 3] for s in seqs_str if s[c * 3: c * 3 + 3] != "---"]
                if codons:
                    aas = [str(Seq(cod).translate()) for cod in codons if len(cod) == 3]
                    trimmed_cons += _Counter(aas).most_common(1)[0][0] if aas else "-"
                else:
                    trimmed_cons += "-"
            motif_pos = trimmed_cons.replace("-", "").find(start_motif)
            if motif_pos < 0:
                print(f"[{family}] WARNING: start anchor '{start_motif}' not found in consensus")
            else:
                pos, start_col = 0, None
                for c, aa in enumerate(trimmed_cons):
                    if aa != "-":
                        if pos == motif_pos:
                            start_col = c
                            break
                        pos += 1
                if start_col is not None:
                    codon_recs = [
                        SeqRecord(Seq(str(r.seq)[start_col * 3:]), id=r.id, description="")
                        for r in codon_recs
                    ]
                    n_after_trim = len(codon_recs[0].seq) // 3
                    print(f"[{family}] Start-anchor trim '{start_motif}': dropped {start_col} codon columns, {n_after_trim} remain")

        out_dir = OUTPUT_DIR / family
        out_dir.mkdir(exist_ok=True)

        aln_path  = out_dir / "alignment.fasta"
        prot_path = out_dir / "protein_alignment.fasta"

        SeqIO.write(codon_recs, aln_path, "fasta")
        SeqIO.write(aligned_prots, prot_path, "fasta")

        print(f"[{family}] {len(codon_recs)} sequences, {n_after_trim} codon positions → {out_dir}/")


if __name__ == "__main__":
    all_families = list(FAMILY_CONFIG.keys())
    families = [f.upper() for f in sys.argv[1:]] if len(sys.argv) > 1 else all_families
    main(families)
