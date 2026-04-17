#!/usr/bin/env python3
"""Convert MSA codon site (in-frame alignment) to PDB residue number.

Usage: python codon_to_residue.py <gene> <codon_site>

Looks in ./<gene>/ for one *.clean.fasta (MSA) and one *.pdb file.
Codon site is 1-indexed, counting all columns including gaps.
"""

import sys
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1 as three_to_one
from Bio.Align import PairwiseAligner


def build_consensus_protein(records: list, n_codons: int) -> tuple[str, list[bool]]:
    """Build a consensus protein from an in-frame MSA.

    Returns:
        protein_seq: gapless consensus amino acid sequence
        is_gap_col:  list of length n_codons; True if column is a consensus gap
                     (i.e. >50% of sequences have a gap at that codon position)
    """
    n_seqs = len(records)
    is_gap_col = []
    protein_chars = []

    for i in range(n_codons):
        codons = [str(r.seq)[i * 3:i * 3 + 3] for r in records]
        non_gap = [c for c in codons if '-' not in c]
        if len(non_gap) <= n_seqs / 2:
            is_gap_col.append(True)
        else:
            aa_counts = Counter(str(Seq(c).translate()) for c in non_gap)
            is_gap_col.append(False)
            protein_chars.append(aa_counts.most_common(1)[0][0])

    return ''.join(protein_chars), is_gap_col


def get_pdb_chains(pdb_file: Path) -> dict[str, tuple[str, list[tuple[int, str]]]]:
    """Return {chain_id: (aa_seq, [(resnum, icode), ...])} from ATOM records."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_file)
    chains = {}
    for model in structure:
        for chain in model:
            seq, resnums = [], []
            for res in chain:
                if res.id[0] == ' ':  # skip HETATM
                    try:
                        seq.append(three_to_one(res.resname))
                    except Exception:
                        seq.append('X')
                    resnums.append((res.id[1], res.id[2].strip()))
            if seq:
                chains[chain.id] = (''.join(seq), resnums)
        break  # first model only
    return chains


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def map_query_pos_to_target(aln, query_pos_0: int) -> int | None:
    """Map 0-indexed query position to 0-indexed target position via alignment blocks.

    Returns None if query position falls in a gap (unaligned) region.
    """
    for (qs, qe), (ts, te) in zip(aln.aligned[0], aln.aligned[1]):
        if qs <= query_pos_0 < qe:
            return ts + (query_pos_0 - qs)
    return None


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <gene> <codon_site>")

    gene = sys.argv[1]
    site = int(sys.argv[2])
    gene_dir = Path(gene)

    if not gene_dir.is_dir():
        sys.exit(f"Error: directory '{gene}' not found")

    fasta_files = list(gene_dir.glob("*.clean.fasta"))
    pdb_files = list(gene_dir.glob("*.pdb"))

    if len(fasta_files) != 1:
        sys.exit(f"Error: expected 1 *.clean.fasta, found {len(fasta_files)}")
    if len(pdb_files) != 1:
        sys.exit(f"Error: expected 1 *.pdb, found {len(pdb_files)}")

    records = list(SeqIO.parse(fasta_files[0], "fasta"))
    if not records:
        sys.exit("Error: no sequences in FASTA file")

    msa_len = len(records[0].seq)
    n_codons = msa_len // 3
    if site < 1 or site > n_codons:
        sys.exit(f"Error: site {site} out of range (MSA has {n_codons} codon columns)")

    pdb_chains = get_pdb_chains(pdb_files[0])
    if not pdb_chains:
        sys.exit("Error: no standard residues found in PDB")

    consensus_protein, is_gap_col = build_consensus_protein(records, n_codons)

    if is_gap_col[site - 1]:
        sys.exit(f"Error: codon site {site} is a consensus gap (>50% of sequences have a gap there)")

    # 1-indexed position of this site in the gapless consensus protein
    aa_pos = sum(1 for i in range(site) if not is_gap_col[i])
    consensus_aa = consensus_protein[aa_pos - 1]

    # Align consensus protein to PDB once; pick best chain
    aligner = make_aligner()
    best_chain, best_aln, best_resnums = None, None, None
    best_score = -float('inf')

    for chain_id, (pdb_seq, pdb_resnums) in pdb_chains.items():
        alignments = aligner.align(consensus_protein, pdb_seq)
        try:
            aln = next(iter(alignments))
        except StopIteration:
            continue
        if aln.score > best_score:
            best_score, best_chain, best_aln, best_resnums = aln.score, chain_id, aln, pdb_resnums

    if best_chain is None:
        sys.exit("Error: could not align consensus protein to any PDB chain")

    target_pos = map_query_pos_to_target(best_aln, aa_pos - 1)
    if target_pos is None or target_pos >= len(best_resnums):
        sys.exit(f"Error: consensus aa position {aa_pos} does not align to any PDB residue")

    resnum, icode = best_resnums[target_pos]
    resnum_str = f"{resnum}{icode}"
    print(f"Gene: {gene}, MSA codon site: {site}")
    print(f"Consensus aa: {consensus_aa}")
    print(f"PDB chain {best_chain}, residue {resnum_str} ({consensus_aa})")
    print("--")
    print("PyMol command to view residue:")
    print(f"select res{resnum_str}, resi {resnum_str}")
    print(f"show sticks, res{resnum_str}")
    print(f"color yellow, res{resnum_str}")


if __name__ == "__main__":
    main()
