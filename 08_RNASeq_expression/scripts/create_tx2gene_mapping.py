#!/usr/bin/env python3
"""
Extract transcript ID to gene symbol mapping from NCBI CDS FASTA files.
Creates a tx2gene mapping file for use with Salmon/tximport.
"""

import re
import argparse
from pathlib import Path


SPECIES_CDS = {
    "HomSap": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/HomSap/data/GCF_009914755.1/cds_from_genomic.fna",
    "PanTro": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/PanTro/data/GCF_028858775.2/cds_from_genomic.fna",
    "PanPan": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/PanPan/data/GCF_029289425.2/cds_from_genomic.fna",
    "GorGor": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/GorGor/data/GCF_029281585.2/cds_from_genomic.fna",
    "PonAbe": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/PonAbe/data/GCF_028885655.2/cds_from_genomic.fna",
    "PonPyg": "/storage/home/kxp5629/group_storage/shared/T2Tv2.assemblies/NCBI_RefSeq/PonPyg/data/GCF_028885625.2/cds_from_genomic.fna",
}

GENE_PATTERN = re.compile(r'\[gene=([^\]]+)\]')


def parse_fasta_headers(fasta_path: Path, species: str) -> list[tuple[str, str, str]]:
    """Extract transcript_id and gene_symbol from FASTA headers."""
    mappings = []

    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                continue

            header = line[1:].strip()
            transcript_id = header.split()[0]

            gene_match = GENE_PATTERN.search(header)
            gene_symbol = gene_match.group(1) if gene_match else "NA"

            mappings.append((transcript_id, gene_symbol, species))

    return mappings


def create_tx2gene_mapping(output_path: Path, include_novel: bool = True, novel_dir: Path = None):
    """Create tx2gene mapping file for all species."""
    all_mappings = []

    for species, cds_path in SPECIES_CDS.items():
        cds_file = Path(cds_path)
        if not cds_file.exists():
            print(f"Warning: CDS file not found for {species}: {cds_path}")
            continue

        print(f"Processing {species}...")
        mappings = parse_fasta_headers(cds_file, species)
        all_mappings.extend(mappings)
        print(f"  Found {len(mappings)} transcripts")

    if include_novel and novel_dir:
        print("\nProcessing novel YAG transcripts...")
        for species in SPECIES_CDS.keys():
            novel_fasta = novel_dir / f"{species}_YAG_cds.fasta"
            if novel_fasta.exists():
                mappings = parse_novel_yag_headers(novel_fasta, species)
                all_mappings.extend(mappings)
                print(f"  {species}: {len(mappings)} novel transcripts")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write("transcript_id\tgene_symbol\tspecies\n")
        for tx_id, gene, species in all_mappings:
            f.write(f"{tx_id}\t{gene}\t{species}\n")

    print(f"\nWrote {len(all_mappings)} mappings to {output_path}")


def parse_novel_yag_headers(fasta_path: Path, species: str) -> list[tuple[str, str, str]]:
    """Parse novel YAG transcript headers (format: Species_GENE_isoX_variant)."""
    mappings = []
    yag_pattern = re.compile(r'^[A-Za-z]+_([A-Z0-9]+)_iso')

    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                continue

            header = line[1:].strip()
            transcript_id = header.split()[0]

            match = yag_pattern.match(transcript_id)
            if match:
                gene_symbol = match.group(1)
            else:
                gene_symbol = transcript_id.split('_')[1] if '_' in transcript_id else transcript_id

            mappings.append((transcript_id, gene_symbol, species))

    return mappings


def main():
    parser = argparse.ArgumentParser(description="Create tx2gene mapping from NCBI CDS FASTA files")
    parser.add_argument(
        "-o", "--output",
        default="/storage/group/kdm16/default/kxp5629/proj/15_RNASeq/data/tx2gene_mapping.tsv",
        help="Output TSV file path"
    )
    parser.add_argument(
        "--novel-dir",
        default="/storage/home/kxp5629/proj/15_RNASeq/data/novel_transcripts",
        help="Directory containing novel YAG transcript FASTAs"
    )
    parser.add_argument(
        "--no-novel",
        action="store_true",
        help="Skip novel YAG transcripts"
    )
    args = parser.parse_args()

    novel_dir = Path(args.novel_dir) if not args.no_novel else None
    create_tx2gene_mapping(
        Path(args.output),
        include_novel=not args.no_novel,
        novel_dir=novel_dir
    )


if __name__ == "__main__":
    main()
