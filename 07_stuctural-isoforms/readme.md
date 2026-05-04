This repository contains Snakemake file with scripts to analyse IsoSeq data from great ape species (bonobo, gorilla, human, chimpanzee, Bornean orangutan, Sumatran orangutan). 

We analyzed data from three different sources:
- targeted IsoSeq (two technical replicates)
- untargeted IsoSeq (one technical replicate)
- Illumina data (one technical replicate)

We performed reference-based assembly with StringTie in Mix mode and then created a union of transcripts between three sets of transcripts (taregeted IsoSeq technical replicate 1, targeted IsoSeq technical replicate 2, untargeted IsoSeq).
