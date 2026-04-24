#!/bin/bash

# Define a list of species
species_list=('bonobo' 'gor' 'bor_orang' 'sum_orang' 'chimp' 'human')

# Full path to gffcompare
gffcompare_path="/galaxy/home/ajg7274/tools/gffcompare/gffcompare/gffcompare"

# Iterate over each species
for species in "${species_list[@]}"; do
    # Define file paths for the current species
    reference_gff="reference/${species}_techrep1.stringtie.annotation.hc_cORF.Ychr.gff"
    targeted1_gff="targeted/${species}_techrep1.stringtie.Ychr.gtf"
    targeted2_gff="targeted/${species}_techrep2.stringtie.Ychr.gtf"
    untargeted_gff="untargeted/${species}_techrep1.stringtie.Ychr.gtf"
    
    # Output directory names
    output_dir_targeted1="${species}_targeted1"
    output_dir_targeted2="${species}_targeted2"
    output_dir_untargeted="${species}_untargeted"
    
    # Run gffcompare for targeted1 GFF
    "$gffcompare_path" -r "$reference_gff" -o "$output_dir_targeted1" "$targeted1_gff"
    
    # Run gffcompare for targeted2 GFF
    "$gffcompare_path" -r "$reference_gff" -o "$output_dir_targeted2" "$targeted2_gff"
    
    # Run gffcompare for untargeted GFF
    "$gffcompare_path" -r "$reference_gff" -o "$output_dir_untargeted" "$untargeted_gff"
done

#human separately

#"$gffcompare_path" -r "reference/human_techrep1.stringtie.annotation.hc_cORF.Ychr.gff" -o "human_targeted1" "targeted/human_techrep1.stringtie.Ychr.gtf"
#"$gffcompare_path" -r "reference/human_techrep1.stringtie.annotation.hc_cORF.Ychr.gff" -o "human_targeted2" "targeted/human_techrep2.stringtie.Ychr.gtf"
