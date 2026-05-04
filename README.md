# How and why multicopy genes survive on the human Y chromosome

Accompanying code and scripts for Pal et al. — integrating T2T genome assemblies, long-read transcriptomics, protein structure predictions, and selection analyses to characterize the seven Y ampliconic gene families (BPY2, CDY, DAZ, HSFY, RBMY, TSPY, VCY) across great apes.

## Repository structure

| Directory | Description |
|-----------|-------------|
| `00_palindrome_identification/` | Palindrome detection with PALINDROVER; tandem array detection with AMPLICOVER |
| `01_multicopy_wholegenome/` | Whole-genome copy number estimation across ape species |
| `02_MSA_trees/` | Multiple sequence alignments (MAFFT) and phylogenetic trees (IQ-TREE) |
| `03_selection_full/` | Genome-wide selection analysis (HYPHY: MG94, BUSTED-E, aBSREL, MEME) |
| `03_selection_on_domains/` | Domain-specific selection analysis (crotonase, chromodomain, RRM, HSF, NAP-S) |
| `04_amplicone_analysis/` | AMPLICONE-based copy number estimation for chimpanzee and gorilla short-read data |
| `04b_copy_number_human/` | Copy number variation across 45 T2T human Y chromosomes |
| `05_pairwise_sequence_ident-permutations_test/` | Pairwise sequence identity among gene copies; permutation tests comparing palindrome vs. array homogenization |
| `06_gene_family_size_n_struct_var/` | Gene family size and structural variant summaries |
| `07_stuctural-isoforms/` | Structural isoform reconstruction from PacBio Iso-Seq + Illumina (STRINGTIE, GFFCOMPARE) |
| `08_RNASeq_expression/` | Expression quantification with Salmon; differential expression (R) |
| `09_seq_struct_prot_clust_ribbon/` | Protein structure clustering (COLABFOLD + FOLDSEEK); Sankey/ribbon visualization |
| `10_SIFT_annotation/` | SIFT functional impact annotation of non-synonymous substitutions |

