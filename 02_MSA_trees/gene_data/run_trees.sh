#!/usr/bin/env bash

IQTREE_CMD="/Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2" 


# -s renamed/VCY_VCX.alignment.trimmed.20240924.renamed.fasta -b 1000
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/RBMY_with_merged_chry_chrx_v6.output.renamed.fasta -b 1000
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/HSFY_alignment_v3.renamed.fasta -b 1000
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/CDY.aligned.renamed.fasta -b 1000
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/BPY2_alignments.renamed.fasta -b 1000
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/CDY.aligned.renamed.fasta -b 1000 
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/HSFY_alignment_v3.renamed.fasta -b 1000 
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/CDY_CDYL_CDYL2.aligned_trimmed.fasta -bb 1000     
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/CDY_CDYL.aligned_trimmed.fasta -bb 1000 

## DAZ repeats

# mkdir -p output

# alignment_file="gene_data/renamed/CDY.aligned_trimmed.renamed.fasta"
# outgroups="SymSyn_chrY_LOC129476678_19362995_Q6.1A_nonC_-,SymSyn_chrY_LOC129476738_21408483_Q6.5B_nonC_t,SymSyn_chrY_LOC129476749_20615009_Q6.2B_nonC_-,SymSyn_chrY_LOC129476751_21915157_Q6.9B_nonC_-,SymSyn_chrY_LOC129476752_20108307_Q6.2A_nonC_t"
# $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/CDY -keep-ident


# alignment_file="gene_data/renamed/CDY_CDYL.aligned_trimmed.renamed.fasta"
# outgroups="GorGor_chr5_CDYL_22394017_NA_nonC_t,HomSap_chr6_CDYL_4575080_NA_nonC_t,PanPan_chr5_CDYL_19400677_NA_nonC_t,PanTro_chr5_CDYL_10660151_NA_nonC_t,PonAbe_chr5_CDYL_4587803_NA_nonC_t,PonPyg_chr5_CDYL_4592538_NA_nonC_t,SymSyn_chr23_CDYL_65930269_NA_nonC_-"
# # -keep-ident flag to keep identical sequences, otherwise he removes some of the DAZL sequences I use as outgrups 
# $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/CDY_CDYL -keep-ident


# alignment_file="gene_data/renamed/BPY2_alignments.renamed.fasta"
# outgroups="GorGor_chrY_LOC129530363_42645394_Q17.1A_C2_t,GorGor_chrY_LOC129530366_43077068_Q17.1B_C2_-"
# $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/BPY2 -keep-ident

# alignment_file="gene_data/renamed/HSFY_alignment_v3.renamed.fasta"
# outgroups="SymSyn_chrY_LOC129476741_21345383_Q6.5B_nonC_t,SymSyn_chrY_LOC129476748_19426501_Q6.1A_nonC_-,SymSyn_chrY_LOC129476753_20045205_Q6.2A_nonC_t,SymSyn_chrY_LOC129476760_20678509_Q6.2B_nonC_-,SymSyn_chrY_LOC129476768_21978653_Q6.9B_nonC_-"
# $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/HSFY -keep-ident

# alignment_file="gene_data/renamed/RBMY_with_merged_chry_chrx_v6.output.renamed.fasta"
# outgroups="GorGor_chrX_RBMX_148580455_nonP-closest-Q18.2B-at-983454_nonC_-,HomSap_chrX_RBMX_135177914_nonP-closest-Q17.6B-at-979846_nonC_-,PanPan_chrX_RBMX_136287749_nonP-closest-Q18B-at-924545_nonC_-,PanTro_chrX_RBMX_134308342_nonP-closest-Q17B-at-918369_nonC_-,PonAbe_chrX_RBMX_142960192_nonP-closest-Q15.1B-at-1624736_nonC_-,PonPyg_chrX_RBMX_141378543_nonP-closest-Q14B-at-1624274_nonC_-,SymSyn_chrX_RBMX_139740568_nonP-closest-Q21.1B-at-957459_nonC_-"
# $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/RBMY_RBMX -keep-ident

# alignment_file="gene_data/renamed/RBMY_with_merged_chry_chrx_v6.output.renamed.fasta.Y-only"
# outgroups="SymSyn_chrY_LOC129476690_19589548_Q6.1A_nonC_t,SymSyn_chrY_LOC129476701_20841562_Q6.2B_nonC_t,SymSyn_chrY_LOC129476740_21171504_Q6.5B_nonC_-,SymSyn_chrY_LOC129476703_19871324_Q6.2A_nonC_-,SymSyn_chrY_LOC129476739_21123360_Q6.2B_nonC_-"
# $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/RBMY -keep-ident


# alignment_file="gene_data/renamed/VCY_VCX.alignment.trimmed.20240924.renamed.fasta"
# outgroups="SymSyn_chrX_LOC129475165_21983856_Q6.9B_nonC_t,SymSyn_chrX_LOC129475166_22219286_nonP-closest-Q7A-at-9303_nonC_-,SymSyn_chrX_LOC129475731_21900132_Q6.9B_nonC_t,SymSyn_chrX_LOC129475809_20612358_Q6.2B_nonC_-,SymSyn_chrX_LOC129475853_20631349_Q6.2B_nonC_-,SymSyn_chrX_LOC129475865_20593405_Q6.2B_nonC_-"
# $IQTREE_CMD -s $alignment_file -bb 1000 -redo  --prefix output/VCY_VCX -keep-ident

# alignment_file="/Users/kxp5629/proj/Y/src/20_align_tree/gene_data/DAZ/renamed/DAZ_DAZL_RBD.aligmnent.20240924.renamed.fasta"
# outgroups="GorGor_chr2_DAZL-RBD1_26423032_NA_nonC_-,HomSap_chr3_DAZL-RBD1_16587918_NA_nonC_-,PanPan_chr2_DAZL-RBD1_16570582_NA_nonC_-,PanTro_chr2_DAZL-RBD1_20590218_NA_nonC_-,PonAbe_chr2_DAZL-RBD1_138140947_NA_nonC_t,PonPyg_chr2_DAZL-RBD1_137333547_NA_nonC_t,SymSyn_chr10_DAZL-RBD1_109881390_NA_nonC_-"
# $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/DAZ_DAZL_RBD -keep-ident

alignment_file="/Users/kxp5629/proj/Y/src/20_align_tree/gene_data/DAZ/renamed/DAZ_RBD.aligmnent.20240924.renamed.fasta"
outgroups="SymSyn_chrY_LOC129476757-RBD1_18127859_Q5B_nonC_t,SymSyn_chrY_LOC129476757-RBD2_18127859_Q5B_nonC_t,SymSyn_chrY_LOC129476758-RBD1_18081065_Q5A_nonC_-,SymSyn_chrY_LOC129476758-RBD2_18081065_Q5A_nonC_-"
$IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/DAZ_RBD -keep-ident


# alignment_file="gene_data/DAZ/renamed/DAZ_DAZL_repeats.alignment.20240924.renamed.fasta"
# outgroups="GorGor_chr2_DAZL-E72_26423032_NA_nonC_-,HomSap_chr3_DAZL-E72_16587918_NA_nonC_-,PanPan_chr2_DAZL-E72_16570582_NA_nonC_-,PanTro_chr2_DAZL-E72_20590218_NA_nonC_-,PonAbe_chr2_DAZL-E72_138140947_NA_nonC_t,PonPyg_chr2_DAZL-E72_137333547_NA_nonC_t,SymSyn_chr10_DAZL-E72_109881390_NA_nonC_-"
# # model="TPM3u+F+G4"
# $IQTREE_CMD -s $alignment_file  -bb 1000 -o $outgroups --prefix output/DAZ_DAZL_repeats -keep-ident

alignment_file="gene_data/DAZ/renamed/DAZ_repeats.alignment.20240924.renamed.fasta"
outgroups="SymSyn_chrY_LOC129476757-E72.1_18127859_Q5B_nonC_t,SymSyn_chrY_LOC129476758-E72.1_18081065_Q5A_nonC_-,SymSyn_chrY_LOC129476757-E72.3_18127859_Q5B_nonC_t,SymSyn_chrY_LOC129476758-E72.3_18081065_Q5A_nonC_-,SymSyn_chrY_LOC129476757-E72.2_18127859_Q5B_nonC_t,SymSyn_chrY_LOC129476758-E72.2_18081065_Q5A_nonC_-"
# model="TPM3u+F+G4"
$IQTREE_CMD -s $alignment_file  -bb 1000 -o $outgroups --prefix output/DAZ_repeats -keep-ident



# # # # alignment_file="/Users/kxp5629/proj/Y/src/20_align_tree/gene_data/TSPY/renamed/TSPY.alignment.20240924.renamed.fasta"
# # # # outgroups="SymSyn_chrY_LOC129476679_19477741_Q6.1A_nonC_-,SymSyn_chrY_LOC129476680_22029903_Q6.9B_nonC_-,SymSyn_chrY_LOC129476681_20729773_Q6.2B_nonC_-,SymSyn_chrY_LOC129476682_22010366_Q6.9B_nonC_-,SymSyn_chrY_LOC129476683_20710209_Q6.2B_nonC_-,SymSyn_chrY_LOC129476684_19458206_Q6.1A_nonC_-,SymSyn_chrY_LOC129476742_21293330_Q6.5B_nonC_t,SymSyn_chrY_LOC129476743_21312997_Q6.5B_nonC_t,SymSyn_chrY_LOC129476754_19993157_Q6.2A_nonC_t,SymSyn_chrY_LOC129476756_20012818_Q6.2A_nonC_t"
# # # # $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/TSPY -keep-ident

# #added missing sequences to TSPY
# alignment_file="/Users/kxp5629/proj/Y/src/20_align_tree/gene_data/TSPY_sequences_20260103/TSPY_20260103_trimmed.renamed.fasta"
# outgroups="SymSyn_chrY_LOC129476679_19477741_Q6.1A_nonC_-,SymSyn_chrY_LOC129476680_22029903_Q6.9B_nonC_-,SymSyn_chrY_LOC129476681_20729773_Q6.2B_nonC_-,SymSyn_chrY_LOC129476682_22010366_Q6.9B_nonC_-,SymSyn_chrY_LOC129476683_20710209_Q6.2B_nonC_-,SymSyn_chrY_LOC129476684_19458206_Q6.1A_nonC_-,SymSyn_chrY_LOC129476742_21293330_Q6.5B_nonC_t,SymSyn_chrY_LOC129476743_21312997_Q6.5B_nonC_t,SymSyn_chrY_LOC129476754_19993157_Q6.2A_nonC_t,SymSyn_chrY_LOC129476756_20012818_Q6.2A_nonC_t"
# $IQTREE_CMD -s $alignment_file -bb 1000 -o $outgroups  --prefix output/TSPY -keep-ident
