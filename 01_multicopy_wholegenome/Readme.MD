### Intro

The script `run_whole_genome.py` was used to identify and collect clusters of multicopy genes across the entire set of genes (from all chromosomes) within the 7 ape lineages (bonobo, chimpanzee, gorilla, human, bornean and sumatran orangutans, and the siamang) - version 2 of assemblies. 


### How to run

#### Dependencies 

The script uses biopython and bcbio GFF packages. Install with `pip` or `conga`. E.g. use:
```
pip install bcbio-gff
```
to install bcbio GFF, which also installs biopython as it's own dependency.

Local installation of **BLAST** is expected. On roar you can do:

```
module load blast
```


To run the script first edit paths in config.json. Make sure they point to the assemblies and annotations of the following accessions: 

```
  HomSap -> GCF_009914755.1
  GorGor -> GCF_029281585.2
  PanPan -> GCF_029289425.2
  PanTro -> GCF_028858775.2
  PonAbe -> GCF_028885655.2
  PonPyg -> GCF_028885625.2
  SymSyn -> GCF_028878055.2
```

Then use the command:

```
run_whole_genome.py -w ./path/to/working/directory
```

### Workflow overview:

1. Subset the gene annotations to only contain the longest isoform of each gene. (time: ~2 hours on MacBook M2, ~20 minutes on 16 a core node on ROAR)
2. Create a blast database of all proteins.
3. Blast the proteins of each species against the blast database. (time: ~10 hours on MacBook M2, ~14 hours on roar 16 cores)
4. Construct graphs of homologous genes. Two genes are considered homologs if they reciprocally share at least 50% sequence identity of at least 35% of protein lengths.
5. Homologous genes were clustered into multi-copy gene families using a simplified single linkage approach (if genes A and B were identified as homologues and so were genes B and C, a group of genes A, B and C was created).


### Output files: 

`clusters.tsv` &rarr; overview table contains the following columns

id|gene_symbol|Description(s)|ampliconic|PanPan|PanTro|HomSap|GorGor|PonPyg|PonAbe|SymSyn
|---|---|---|---|---|---|---|---|---|---|---|
cluster ID| gene symbol from one of the gene copies| Descriptions collected from all gene copies| yes/no* | Bonobo gene counts| Chimpanzee gene counts| Human gene counts| Gorilla gene counts| Bornean orangutan gene counts| Sumatran orangutan gene counts| Siamang gene counts

<small> *yes if at least one pairwise identity of gene copies within one species > 97%. </small>


`cluster_XX.faa` &rarr; One fasta file per gene cluster.

### Post processing 

To get a table of chromosome numbers where the multicopy genes occur run:

```
python get_chromosome_counts.py -c {clusters_merged directory} -g {config json file} > {output.table}
e.g.:
python get_chromosome_counts.py -c ~/scratch/workdir/protein_extracted_longest/clusters_merged/ -g roarconfig.json > chromosomes.tsv
```
