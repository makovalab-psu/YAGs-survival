MERGED_#!/bin/bash

# List of species
species=("bonobo" "gor" "bor_orang" "sum_orang" "chimp")

# Path to gffcompare
gffcompare="/galaxy/home/ajg7274/tools/gffcompare/gffcompare/gffcompare"

# Loop through each species
for sp in "${species[@]}"
do
    # Merge gff files
    $gffcompare stringtie_targeted/${sp}_techrep1.stringtie.Ychr.gtf \
                 stringtie_targeted/${sp}_techrep2.stringtie.Ychr.gtf \
                 stringtie_untargeted/${sp}_techrep1.stringtie.Ychr.gtf -o MERGED_${sp}
done

#human

#filter by Y chromosome
awk '$1 == "NC_060948.1"' /galaxy/home/ajg7274/isoseq_apes/targeted_isoseq/snakemake-isoseq-full-genome/stringtie/human1_techrep1.stringtie.gtf > /galaxy/home/ajg7274/isoseq_apes/targeted_isoseq/snakemake-isoseq-full-genome/stringtie/human1_techrep1.stringtie.Ychr.gtf
awk '$1 == "NC_060948.1"' /galaxy/home/ajg7274/isoseq_apes/targeted_isoseq/snakemake-isoseq-full-genome/stringtie/human1_techrep2.stringtie.gtf > /galaxy/home/ajg7274/isoseq_apes/targeted_isoseq/snakemake-isoseq-full-genome/stringtie/human1_techrep2.stringtie.Ychr.gtf

$gffcompare /galaxy/home/ajg7274/isoseq_apes/targeted_isoseq/snakemake-isoseq-full-genome/stringtie/human1_techrep1.stringtie.Ychr.gtf /galaxy/home/ajg7274/isoseq_apes/targeted_isoseq/snakemake-isoseq-full-genome/stringtie/human1_techrep2.stringtie.Ychr.gtf -o MERGED_human

# Create folder and move files
# mkdir stringtie_merge/
mv MERGED* stringtie_merge/
mv /galaxy/home/ajg7274/isoseq_apes/mix_union_targeted_untargeted/stringtie_targeted/MERGED* stringtie_merge/
