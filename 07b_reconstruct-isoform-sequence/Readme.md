# Estimate sequence diversity of YAGs from IsoSeq reads.

## Pipeline and data

The analysis consists of the following steps:

  1.  Identify reads belonging to a YAG family by mapping IsoSeq data to the
whole genome and collecting reads overlapping annotated genes from one YAG family.
  2.  Collect all reads identified as belonging to one gene family and map them
  to a single representative copy.
  3.  Map all sequences of annotated genes (*reference sequences* - extracted from the assembly) from one gene family to the representative copy.
  4.  Extract all unique positions present in the annotated gene (reference) sequences.
  5.  Identify additional unique positions based on variant calling on the IsoSeq reads.
  6.  Use unique positions to reconstruct a *signature* sequence.
  7.  Use the signature sequence to characterize and group each read.


The data folder contains:

  1.  Bed file with the locations of all annotated copies for each gene family.
  2.  Single representative gene sequence.
  3.  Bed file with annotation for the representative sequence (exon locations).
  4.  Fasta file with coding sequences of all annotated gene copies.

Additional data needed:

  1. Genomic references:
  ```
    HomSap -> GCF_009914755.1
    GorGor -> GCF_029281585.2
    PanPan -> GCF_029289425.2
    PanTro -> GCF_028858775.2
    PonAbe -> GCF_028885655.2
    PonPyg -> GCF_028885625.2
    SymSyn -> GCF_028878055.2
  ```
  2. IsoSeq sequencing reads:
  Most reads can be downloaded through the `download_reads.sh` script (requires path to sratoolkit).
  Untargeted human sequencing reads are not on SRA yet.


## How to run.

Update the config.yaml file so that the paths are valid.

All steps of the analysis are coordinated by the `pipeline.py` script, and one gene family
per species can be run with a single command. (See wrapper.sh script for examples).

## Output

1. The main output of the analysis is a table with the occurrences of each signature per input sample (each replicate):

|signature|sample1|sample2|sum|
|---------|-------|-------|---|
|AGCT|6|7|15|
|ACCT|5|4|9|

2. A table is produced where each annotated (reference) gene sequence is annotated with a signature:

|ID|coding sequnece| signature|
|--|---------------|----------|
|PanTro_chrY_YAG1_CDS|	ATGGTAGAAGCAGATCATCCTGGCAAGCTTTTCATTGGT|	AGCT|
|PanTro_chrY_YAG1_CDS|	ATGCTAGAAGCAGATCATCCTGGCAAGCTTTTCATTGGT|	ACCT|

3... The analysis produces one json file per sample. This is then used to produce the main output table (1).


## Requirements:

### Tools:
 - minimap2
 - samtools
 - bedtools
 - bcftools
 - freebayes

### Python packages:
(I'm using python 3.11)
  - argparse
  - pysam
  - pyvcf
