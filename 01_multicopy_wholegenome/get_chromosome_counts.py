import argparse
import os
import sys
import json
from collections import defaultdict

def load_config(config_file):
    with open(config_file, 'r') as f:
        return json.load(f)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cluster_dir', help="path to working directory where intermediate and final results will be stored", default="./merged_clusters/")
    parser.add_argument('-g', '--config', help="path to config file with paths to assemblies and annotations", default="./config.json")
    args = parser.parse_args()

    cluster_dir = args.cluster_dir

    config = load_config(args.config)

    data = config["data"]
    species_list = [d['species'] for d in data ]

    #read helper files to get dictionary chrID -> chr nubmer (NC_073224.2 -> 1)
    chr_id_dict = {}
    print("=== load helper files ===", file=sys.stderr)
    for species in species_list:
        print(f"   {species}", file=sys.stderr)
        chr_map_file = f"./helpers/sequence_report_{species}.tsv"
        with open(chr_map_file,'r') as read_file:
            for line in read_file:
                if line.startswith("Assembly"):
                    continue

                line = line.strip()
                line = line.split("\t")
                chr_id_dict[line[8]] = line[2]
            
    
    prot_id_to_chr_dict = {}

    # from the gff files I need to collect a dictionary pointing me from gene "XP" ID to chromosome.
    print("=== process annotation ===", file=sys.stderr)
    for species in species_list:
        print(f"   {species}", file=sys.stderr)
        path_to_gff = [d['data'] for d in data if d['species'] == species][0]['path_to_annotation_NCBI']

        item_separator = ";"
        key_value_separator = "="

        with open(path_to_gff, 'r') as gff_file:
            for line in gff_file:
                line = line.strip()

                if line.startswith("#"):
                    continue

                line = line.split("\t")

                # Split the string into individual key-value pair strings
                pairs = line[8].split(item_separator)

                # Use dictionary comprehension to split each pair and create the dictionary
                result_dict = {
                    key_value_pair.split(key_value_separator)[0]: key_value_pair.split(key_value_separator)[1]
                    for key_value_pair in pairs
                }
                if "protein_id" in result_dict.keys():
                    prot = result_dict['protein_id']
                    chrom = line[0]
                    prot_id_to_chr_dict[prot] = chr_id_dict[chrom]


    #construct and print header
    header_list = []
    header_list.append("cluster_id")
    header_list.append("gene_name")
    header_list.append("gene_desc")

    for species in species_list:
        header_list.append(species)

    header = "\t".join(header_list)
    print(header)

    lines = []

    # for each cluster in clusters.tsv get gene name, description, cluster ID, use cluster ID to open
    # corresponding cluster fasta fila where for each gene translate geneID to chromosome number
    # collect and print
    print("=== processing cluster files ===", file=sys.stderr)
    with open(f"{cluster_dir}/clusters.tsv", 'r') as clustsers_file:
        cluster_dict = {}
        for line in clustsers_file:
            line = line.strip()
            if line.startswith("id"):
                continue
        
            line = line.split("\t")

            cluster_id = line[0]
            gene_name = line[1]
            gene_desc = line[2]

            cluster_dict["cluster_id"] = cluster_id
            cluster_dict['gene_name'] = gene_name
            cluster_dict['gene_desc'] = gene_desc
        
            species_dict = defaultdict(list)

            with open(f"{cluster_dir}/cluster_{cluster_id}.faa", 'r') as single_cluster:
                for line in single_cluster:
                    line = line.strip()
                    if line.startswith(">"):
                        line = line.split("::")
                        gene_id = line[0].replace(">","")
                        chromosome = prot_id_to_chr_dict[gene_id]
                        species = line[3].split(" ")[0]
                        # print(f"{species} -> {gene_id}")
                        species_dict[species].append(chromosome)

            lines.append({**cluster_dict, **species_dict})

        for line in lines:
            out_list = []
            for species in species_list:
                if species in line.keys():
                    out_list.append( ",".join(sorted(line[species])))
                else:
                    out_list.append("")

            out_string = "\t".join(out_list)
            print(f"{line['cluster_id']}\t{line['gene_name']}\t{line['gene_desc']}\t{out_string}")


if __name__ == "__main__":
    main()
