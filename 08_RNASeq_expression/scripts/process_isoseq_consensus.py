#!/usr/bin/env python3
"""
Process IsoSeq reads to generate consensus sequences and extract CDS regions.

Requirements:
    - BioPython: pip install biopython
    - SPOA: conda install -c bioconda spoa
    - exonerate: conda install -c bioconda exonerate

Usage:
    python process_isoseq_consensus.py --input-dir <path> --protein-fasta <path> --output-dir <path>
"""

import argparse
import os
import re
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def parse_filename(filepath):
    """Extract species, gene, fingerprint from filename."""
    basename = os.path.basename(filepath)

    # Non-DAZ: Species_Gene_ORFs.fa_Gene_isoform_N_FINGERPRINT.fasta
    match = re.match(r'([A-Za-z]+)_([A-Z0-9]+)_ORFs\.fa_\2_isoform_(\d+)_(.+)\.fasta', basename)
    if match:
        return {'species': match[1], 'gene': match[2], 'isoform': match[3],
                'fingerprint': match[4], 'type': 'fasta', 'path': filepath}

    # DAZ: Species_DAZ-tail.REPEAT-PATTERN.fq
    match = re.match(r'([A-Za-z]+)_DAZ-tail\.(.+)\.fq', basename)
    if match:
        return {'species': match[1], 'gene': 'DAZ', 'isoform': None,
                'fingerprint': match[2], 'type': 'fastq', 'path': filepath}
    return None


def parse_protein_header(header, gene_from_file=None):
    """
    Parse protein header formats:
    - Standard: Species::Gene_isoform_N::Fingerprint::N_reads[_modORF]
      Example: GorGor::CDY_isoform_2::CCGAAGACCGTA::53_reads
    - DAZ: Species::Fingerprint::N_reads
      Example: GorGor::C-J-G-B-I-H-B-D-A-A-E-A-F::54_reads
    """
    parts = header.split('::')
    if len(parts) < 2:
        return None

    species = parts[0]

    # Check if second part is gene_isoform or fingerprint
    if '_isoform_' in parts[1]:
        # Standard format: Species::Gene_isoform_N::Fingerprint::N_reads
        gene_iso = parts[1]
        fingerprint = parts[2] if len(parts) > 2 else ''
        gene_match = re.match(r'([A-Z0-9]+)_isoform_(\d+)', gene_iso)
        if gene_match:
            gene = gene_match[1]
            isoform = gene_match[2]
        else:
            gene = gene_iso
            isoform = None
    else:
        # DAZ format: Species::Fingerprint::N_reads (no gene in header)
        fingerprint = parts[1]
        gene = gene_from_file if gene_from_file else 'DAZ'
        isoform = None

    return {
        'species': species,
        'gene': gene,
        'isoform': isoform,
        'fingerprint': fingerprint,
        'full_header': header
    }


def find_input_files(input_dir):
    """Find all FASTA/FASTQ files to process."""
    files = []
    for pattern in ['*_ORFs.fa_*.fasta', '*_DAZ-tail.*.fq']:
        for f in Path(input_dir).rglob(pattern):
            info = parse_filename(str(f))
            if info:
                files.append(info)
    return files


def run_spoa(input_file, file_type='fasta'):
    """Run SPOA to generate consensus. Returns consensus sequence string."""
    # Convert FASTQ to FASTA if needed
    if file_type == 'fastq':
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            for rec in SeqIO.parse(input_file, 'fastq'):
                tmp.write(f">{rec.id}\n{str(rec.seq)}\n")
            input_for_spoa = tmp.name
    else:
        input_for_spoa = input_file

    # print(input_for_spoa)

    result = subprocess.run(['spoa', input_for_spoa, '-r', '0'],
                            capture_output=True, text=True)

    # Cleanup temp file
    if file_type == 'fastq':
        os.remove(input_for_spoa)

    if result.returncode != 0:
        return None

    # Parse FASTA output from SPOA
    lines = result.stdout.strip().split('\n')
    return ''.join(lines[1:]) if len(lines) > 1 else None


def run_exonerate(protein_seq, dna_seq, seq_id):
    """
    Use exonerate protein2dna to find CDS region.
    Returns (cds_sequence, start, end, strand) or None.
    """

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as prot_tmp:
        prot_tmp.write(f">protein\n{protein_seq}\n")
        prot_file = prot_tmp.name

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as dna_tmp:
        dna_tmp.write(f">{seq_id}\n{dna_seq}\n")
        dna_file = dna_tmp.name

    try:
        result = subprocess.run([
            'exonerate', '--model', 'protein2dna',
            '--query', prot_file, '--target', dna_file,
            '--showalignment', 'no', '--showvulgar', 'no',
            '--ryo', '%tab\t%tae\t%qal\t%tal\n'  # target begin, end, query alignment, target alignment
        ], capture_output=True, text=True)

        if result.returncode != 0 or not result.stdout.strip():
            return None

        # Parse best hit (first line)
        line = result.stdout.strip().split('\n')[3]

        parts = line.split('\t')
        if len(parts) >= 2:
            start, end = int(parts[0]), int(parts[1])
            if start > end:  # Reverse strand
                cds = str(Seq(dna_seq[end:start]).reverse_complement())
                return cds, end, start, '-'
            else:
                return dna_seq[start:end], start, end, '+'

    finally:
        os.remove(prot_file)
        os.remove(dna_file)


    return None


def generate_id(info):
    """Generate clean sequence ID."""
    fp = info['fingerprint'].replace('.', '').replace('_', '-')
    if info.get('isoform'):
        return f"{info['species']}_{info['gene']}_iso{info['isoform']}_{fp}"
    return f"{info['species']}_{info['gene']}_{fp}"


def load_proteins(protein_dir):
    """
    Load all protein sequences from directory.
    Returns dict: (species, gene, fingerprint) -> (protein_seq, full_header)
    """
    proteins = {}

    # Find all .fa files in directory
    for fa_file in Path(protein_dir).glob('*.fa'):
        # Extract gene from filename: collect_all_iso_GENE.fa
        fname = fa_file.stem  # e.g., "collect_all_iso_CDY"
        gene_from_file = None
        if fname.startswith('collect_all_iso_'):
            gene_from_file = fname.replace('collect_all_iso_', '').replace('_fixed', '')

        for rec in SeqIO.parse(fa_file, 'fasta'):
            parsed = parse_protein_header(rec.id, gene_from_file)
            if parsed:
                # Create lookup key
                key = (parsed['species'], parsed['gene'], parsed['fingerprint'], parsed['isoform'])
                proteins[key] = (str(rec.seq), rec.id)

                # Also store with dots replaced (file names have dots replaced)
                fp_nodots = parsed['fingerprint'].replace('.', '')
                if fp_nodots != parsed['fingerprint']:
                    key2 = (parsed['species'], parsed['gene'], fp_nodots)
                    proteins[key2] = (str(rec.seq), rec.id)

    return proteins


def main():
    parser = argparse.ArgumentParser(description='Process IsoSeq to consensus + CDS')
    parser.add_argument('--input-dir', required=True, help='Input directory with FASTA/FASTQ files')
    parser.add_argument('--protein-dir', required=True, help='Directory with protein FASTA files')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load protein sequences
    print("Loading protein sequences...")
    proteins = load_proteins(args.protein_dir)
    print(f"  Loaded {len(proteins)} protein entries")

    # Find input files
    print("Finding input files...")
    input_files = find_input_files(args.input_dir)
    print(f"  Found {len(input_files)} files")

    # Process each file
    results_by_species = defaultdict(list)
    stats = {'matched': 0, 'no_protein': 0, 'no_cds': 0, 'spoa_failed': 0}

    for i, info in enumerate(input_files):
        # print(info)
        # assert False
        seq_id = generate_id(info)
        print(f"[{i+1}/{len(input_files)}] {seq_id}")

        # Step 1: SPOA consensus
        consensus = run_spoa(info['path'], info['type'])
        if not consensus:
            print(f"  SKIP: SPOA failed")
            stats['spoa_failed'] += 1
            continue

        # Step 2: Find matching protein
        # Build lookup key from file info
        fp = info['fingerprint']
        fp_clean = fp.replace('.', '')

        protein_seq = None
        protein_header = None
        # print(proteins)
        # assert False
        for key in [
            (info['species'], info['gene'], fp, info['isoform']),
            (info['species'], info['gene'], fp_clean, info['isoform']),
        ]:
            if key in proteins:
                protein_seq, protein_header = proteins[key]
                break

        if protein_seq:
            cds_result = run_exonerate(protein_seq, consensus, seq_id)
            if cds_result:
                cds_seq, start, end, strand = cds_result
                print(f"  CDS: {len(cds_seq)}bp ({strand} strand)")
                results_by_species[info['species']].append((seq_id, cds_seq, consensus))
                stats['matched'] += 1
            else:
                print(f"  WARN: exonerate found no match, using full consensus")
                results_by_species[info['species']].append((seq_id, consensus, consensus))
                stats['no_cds'] += 1
        else:
            print(f"  WARN: No protein match for ({info['species']}, {info['gene']}, {fp})")
            results_by_species[info['species']].append((seq_id, consensus, consensus))
            stats['no_protein'] += 1

    # Write output files per species
    print("\nWriting output files...")
    for species, seqs in results_by_species.items():
        # CDS file
        cds_file = os.path.join(args.output_dir, f"{species}_YAG_cds.fasta")
        with open(cds_file, 'w') as f:
            for seq_id, cds, _ in seqs:
                f.write(f">{seq_id}\n{cds}\n")
        print(f"  {cds_file} ({len(seqs)} sequences)")

        # Full consensus file
        cons_file = os.path.join(args.output_dir, f"{species}_YAG_consensus.fasta")
        with open(cons_file, 'w') as f:
            for seq_id, _, cons in seqs:
                f.write(f">{seq_id}\n{cons}\n")
        print(f"  {cons_file}")

    # Print summary statistics
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    print(f"  Total files processed: {len(input_files)}")
    print(f"  Successfully matched:  {stats['matched']}")
    print(f"  No protein match:      {stats['no_protein']}")
    print(f"  No CDS found:          {stats['no_cds']}")
    print(f"  SPOA failed:           {stats['spoa_failed']}")
    print("="*50)
    print("Done!")


if __name__ == '__main__':
    main()
