#!/bin/bash

base_folder="/galaxy/home/ajg7274/isoseq_apes/mix_union_targeted_untargeted/compare_isoforms"

gene_fam=$1

folder="$base_folder/$gene_fam"

for i in $(seq 1 $(ls "$folder"/transcript_*.gff | wc -l))
do
    # Run gffcompare for each transcript file against the reference file
    /galaxy/home/ajg7274/tools/gffcompare/gffcompare/gffcompare "$folder"/"${gene_fam}_transcripts.gff" -r "$folder"/"transcript_${i}.gff" -o "transcript_${i}"
done
