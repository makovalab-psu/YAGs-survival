#!/usr/bin/env bash

source /storage/home/kxp5629/miniconda3/etc/profile.d/conda.sh
conda activate 09_SD

mkdir -p logs

SPECIES_LIST=("HomSap" "PanPan" "PanTro" "PonAbe" "PonPyg" "GorGor")
# SPECIES_LIST=("GorGor")
# SPECIES_LIST=("PonAbe")

# for species in ${SPECIES_LIST[@]};do
#     echo $species
#     # ./pipeline.py --speces "$species" --gene TSPY --config configy.yaml
#     ./pipeline.py --species "$species" --gene TSPY --config config.yaml > "logs/${species}_TSPY.log" 2>&1
#     echo $species done
# done

# SPECIES_LIST=("HomSap")
SPECIES_LIST=("HomSap" "PanPan" "PanTro" "PonAbe" "PonPyg" "GorGor")
# # SPECIES_LIST=("PanTro" "GorGor")
# 

# SPECIES_LIST=("HomSap")

GENE="DAZ"
for species in ${SPECIES_LIST[@]};do
    echo $species
    ./pipeline.py --species "$species" --gene "$GENE" --config config.yaml  > "logs/${species}_${GENE}.log" 2>&1
    echo $species done
done

# GENE="RBMY"
# #GENE="TSPY"
# for species in ${SPECIES_LIST[@]};do
#     echo $species
#     ./pipeline.py --species "$species" --gene "$GENE" --config config.yaml --isoform > "logs/${species}_${GENE}.log" 2>&1
#     echo $species done
# done

# SPECIES_LIST=("PonAbe" "PonPyg") 
# GENE="RBMYB"
# 
# for species in ${SPECIES_LIST[@]};do
#     echo $species
#     ./pipeline.py --species "$species" --gene "$GENE" --config config.yaml --isoform > "logs/${species}_${GENE}.log" 2>&1
#     echo $species done
# done

# SPECIES_LIST=("HomSap" "PanPan" "PanTro" "PonAbe" "PonPyg" "GorGor")
# # SPECIES_LIST=("HomSap" "PonAbe" "PonPyg") 
# GENE="CDY"
# 
# for species in ${SPECIES_LIST[@]};do
#     echo $species
#     ./pipeline.py --species "$species" --gene "$GENE" --config config.yaml --isoform > "logs/${species}_${GENE}.log" 2>&1
#     echo $species done
# done

# SPECIES_LIST=("HomSap" "PonAbe" "PonPyg" "GorGor")
# # SPECIES_LIST=("HomSap") 
# GENE="HSFY"
# 
# for species in ${SPECIES_LIST[@]};do
#     echo $species
#     ./pipeline.py --species "$species" --gene "$GENE" --config config.yaml --isoform > "logs/${species}_${GENE}.log" 2>&1
#     echo $species done
# done

# SPECIES_LIST=("PanPan" "PanTro" "GorGor" "HomSap")
# # SPECIES_LIST=("PanPan")
# 
# GENE="BPY2"
# 
# for species in ${SPECIES_LIST[@]};do
#     echo $species
#     ./pipeline.py --species "$species" --gene "$GENE" --config config.yaml --isoform > "logs/${species}_${GENE}.log" 2>&1
#     echo $species done
# done

# SPECIES_LIST=("PanPan" "PanTro" "HomSap" )
# SPECIES_LIST=("PanTro")
# GENE="VCY"
# 
# for species in ${SPECIES_LIST[@]};do
#     echo $species
#     ./pipeline.py --species "$species" --gene "$GENE" --config config.yaml --isoform > "logs/${species}_${GENE}.log" 2>&1
#     echo $species done
# done
