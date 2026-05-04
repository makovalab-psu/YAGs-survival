from BCBio import GFF
from BCBio.GFF import GFFExaminer
examiner = GFFExaminer()

import argparse
from Bio import SeqIO
import copy
import json 
from datetime import datetime as dt
import os
import subprocess


class TimeTracker:
    def __init__(self):
        self.start_time = dt.now()
        self.prev_time = dt.now()
    
    def print_time_diff(self, message=None):
        current_time = dt.now()
        overall_time = current_time - self.start_time
        prev_time_diff = current_time - self.prev_time

        print(f"   ----------")
        print(f"   Overall time elapsed: {overall_time}")
        print(f"   Time since last checkpoint: {prev_time_diff}")
        print(f"   Current time: {current_time}")
        if message:
            print(f"   ----------")
            print(f"   next task: {message}")
        print(f"   ----------")

        self.prev_time = current_time
    
    def reset(self):
        self.start_time = dt.now()
        self.prev_time = dt.now()
        
def load_config(config_file):
    with open(config_file, 'r') as f:
        return json.load(f)
    
def get_annotation_file(species, data):
    return  [d for d in data if d['species'] == species][0]['data']['path_to_annotation_NCBI'] 

def get_single_trans_annot_file(species, data, work_dir):    
    os.makedirs(f"{work_dir}/subset_annotation", exist_ok=True)

    annotation_file_orig =  get_annotation_file(species, data)
    annotation_file_single_trans = work_dir + "/subset_annotation/" + "_".join(annotation_file_orig.split("/")[-2:]) 
    return annotation_file_single_trans.replace(".gff", ".protcoding_singletrans.gff") 

def process_species_annotation(species, data, work_dir):

    annotation_file_orig = get_annotation_file(species, data)
    annotation_file_single_trans = get_single_trans_annot_file(species, data, work_dir)
    print(f"  {species}")
    subset_annotation(annotation_file_orig, annotation_file_single_trans)

def subset_annotation(annotation_file, subset_file):
    handle = open(annotation_file)

    records = GFF.parse(handle)
    recs = []
    for rec in records:
        # protss = set()
        rec_out = copy.copy(rec)
        rec_out.features = []
        

        for feature in rec.features:

            if feature.type != "gene":
                continue

            if not "gene_biotype" in feature.qualifiers:
                continue

            if not "protein_coding" in feature.qualifiers["gene_biotype"]:
                continue

            feature_out = copy.copy(feature)
            feature_out.sub_features = []
            has_mrna = False
            subf_out = None
            for subf in feature.sub_features:


                max_len = 0
                if subf.type == "mRNA":
                    
                    length = 0
                    this_out = copy.copy(subf)

                    for subsubf in subf.sub_features:
                        if subsubf.type != "CDS":
                            continue
                        length += len(subsubf)

                    if length > max_len:
                        max_len = length
                        subf_out = copy.copy(this_out)
                        has_mrna = True

            if (has_mrna ):
                feature_out.sub_features.append(subf_out)
                rec_out.features.append(feature_out)
        recs.append(rec_out)    
    with open(subset_file, "w") as out_handle:
        GFF.write(recs, out_handle)
    handle.close()


def report(clusters,prot_record_dict, counter, total_len,prot_to_species, min_identity_dictionaries, species_list, work_dir):

    results = {}
    cluster = clusters[counter]

    if counter %100 == 0:
        print(f"    {counter}/{total_len}")

    cluster_copy = copy.deepcopy(cluster)

    gene_names = set()
    gene_symbol = ""
    gene_identity_counter_per_species = {}
    max_identity = 0.0
    report = False

    species_count = {}
    for species in species_list:
        species_count[species] = 0
    

    cluster_fasta_file = f"{work_dir}/protein_extracted_longest/clusters_merged/cluster_{counter}.faa"
    with open(cluster_fasta_file, "w") as cluster_outfile:
        for prot in cluster:
            species = ""
            desc = ""

            element = prot_record_dict[prot]
            gene_names.add(element["name"])
            species = element["species"]
            species_count[species] += 1
            # print(element)
            desc = element["name"]
            gene = element["gene_id"]
            chromosome = element["chromosome"]
            report_gene = gene

            
            sequence = element["seq"]

            el_start =  element["start"]
            el_end =  element["end"]

            if species == "HomSap":
                if not gene.startswith("LOC"):
                    gene_symbol = gene

            cluster_copy.remove(prot)

            print(f">{prot}::{report_gene}::chr{chromosome}::{species} description:{desc}\n{sequence}", file=cluster_outfile)   

            for prot_B in cluster_copy:

                species = prot_to_species[prot]

                if species != prot_to_species[prot_B]:
                    continue

                key_pair = tuple(sorted([prot, prot_B]))
                if not key_pair in min_identity_dictionaries[species].keys():
                    continue

                identity = min_identity_dictionaries[species][key_pair]
                if species not in gene_identity_counter_per_species:
                    gene_identity_counter_per_species[species] = [identity]
                else:
                    gene_identity_counter_per_species[species].append(identity)
  
        
    if(sorted(species_count.values()))[-1] > 1:
        report = True  

    if report:
        
        for species in gene_identity_counter_per_species:
            max_identity_in_species = sorted(gene_identity_counter_per_species[species])[-1]
            if max_identity_in_species > max_identity:
                max_identity = max_identity_in_species

        if gene_symbol == "":
            gene_symbol = list(gene_names)[0]
        ampliconic = "no"
        if max_identity > 97.0:
            ampliconic = "yes"
        line = f"{counter}\t{gene_symbol}\t{';'.join(gene_names)}\t{ampliconic}"
        for species in species_list:
            line += f"\t{species_count[species]}"

        results['file_line'] = line

        identity_line = f"{counter}\t{gene_symbol}\t{';'.join(gene_names)}"
        for species in species_list:
            if species not in gene_identity_counter_per_species:
                gene_identity_counter_per_species[species] = [0]
            
            identities_str = [str(x) for x in gene_identity_counter_per_species[species]]
            identity_line += f"\t{';'.join(identities_str)}"
        results['identity_line'] = identity_line

    results['report'] =report
    return results

def main():
    
    #Code transcribed from jupyter notebook (thats why it's a huge monolith)

    time_tracker = TimeTracker()

    parser = argparse.ArgumentParser(
        prog="run_whole_genome.py",
        description="Search for multicopy gene from annotated protein coding sequences from the T2T assembly."
    )
    parser.add_argument('-s', '--skip_blast', action='store_true')
    parser.add_argument('-c', '--chromosome_subset', help='If present run program only on a subst of the GFF file. E.g. "chr_1"')
    parser.add_argument('-w', '--work_dir', help="path to working directory where intermediate and final results will be stored", default="./work_dir/cluster_genes")
    parser.add_argument('-g', '--config', help="path to config file with paths to assemblies and annotations", default="./config.json")
    parser.add_argument('-j', '--jobs', help="number of cores for blast", default="8")
    args = parser.parse_args()

    cores = int(args.jobs)
    work_dir = args.work_dir

    record_dict = {}

    config = load_config(args.config)

    data = config["data"]
    chromosomes_dict = config["chromosomes_dict"]
    species_list = [d['species'] for d in data]

    time_tracker.print_time_diff("subset annotationfiles:")
       
    # subsettting runs for approx 2 hours (on my M2) - results can be reused
    for species in species_list:
        annotation_file = get_annotation_file(species, data)
        annotation_file_single_trans = get_single_trans_annot_file(species, data, work_dir)
        print(f"    {species} start: {dt.now()}")
        subset_annotation(annotation_file, annotation_file_single_trans)

    time_tracker.print_time_diff("create record dictionary:")

    for species in species_list:
        print(f"     {species}")

        annotation_file_orig = get_single_trans_annot_file(species, data, work_dir)
        annotation_file = annotation_file_orig

        if args.chromosome_subset != None:

            chromosome = args.chromosome_subset
            annotation_file = annotation_file_orig.replace(".gff", f"_{chromosome}.gff")
            chr_id = [d for d in data if d['species'] == species][0]['data'][chromosome]

            if not os.path.isfile(annotation_file):
                with open(annotation_file_orig,'r') as in_file:
                    with open(annotation_file,"w") as outfile:
                        lines = in_file.readlines()
                        for line in lines:
                            line = line.strip()
                            if line.startswith("#"):
                                print(line,file=outfile)
                                continue

                            id = line.strip().split("\t")[0]
                            if id != chr_id:
                                continue
                            
                            print(line,file=outfile)
            
        handle = open(annotation_file)

        records = GFF.parse(handle)
        for rec in records:
            protss = set()
            for feature in rec.features:
                if feature.type != "gene":
                    continue
                if "description" in feature.qualifiers:
                    description = feature.qualifiers["description"][0]

                    gene_id = feature.id.replace("gene-", "")
                    gene_spec_id = gene_id + "_" + species
                        
                    # print(feature.qualifiers)
                    for subf in feature.sub_features:
                        if subf.type == "mRNA":
                            mrna_id = subf.id.replace("rna-", "")
                            for subsubf in subf.sub_features:
                                if subsubf.type != "CDS":
                                    continue
                                # cds_id = subsubf.id.replace("cds-", "")
                                prot_id = subsubf.qualifiers["protein_id"][0]
                                if prot_id in protss:
                                    continue

                                protss.add(prot_id)
                                element = {}
                                element["species"] = species
                                element["name"] = description
                                element["id"] = prot_id
                                element["mrna_id"] = mrna_id
                                element["gene_id"] = gene_id
                                element["chromosome"] = chromosomes_dict[rec.id]
                                element["start"] = int(subsubf.location.start)
                                element["end"] = int(subsubf.location.end)
                                element["gene_spec_id"] = gene_spec_id

                                # element["seq"] = str(subf.extract(rec.seq))
                                if gene_id in record_dict:
                                    record_dict[gene_id].append(element)
                                else:
                                    record_dict[gene_id] = [element]

        handle.close()

    prot_dicts = {}
    time_tracker.print_time_diff("extract sequences")
    for species in species_list:
        print(f"     {species}")
        sequence_file = [d for d in data if d['species'] == species][0]['data']['prot']
        sequence_handle = open(sequence_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(sequence_handle, "fasta"))

        prot_dicts[species] = seq_dict

        sequence_handle.close()

    #load protein sequences
    time_tracker.print_time_diff("load protein sequences")
    for record in record_dict:
        for element in record_dict[record]:
            species = element["species"]
            id = element["id"]
            name = element["name"]
            prot_sequences = [prot_dicts[species][x] for x in prot_dicts[species] if (id in x )]
            if len(prot_sequences) == 0:
                print("No sequence found")
                print(id)
                assert False
            seq = str(prot_sequences[0].seq)
            if len(seq) > 1:
                for prot in prot_sequences:
                    if len(prot.seq) > len(seq):
                        seq = str(prot.seq)

            element["seq"] = seq
            
    # reduce list of prot sequences to longest sequence per gene
    time_tracker.print_time_diff("build gene records")
    gene_record_dict = {}
    for record in record_dict:
        gene_spec_ids = set()
        gene_spec_elements = []
        
        for element in record_dict[record]:
            gene_spec_id = element["gene_spec_id"]
            
            if gene_spec_id not in gene_spec_ids:
                gene_spec_elements.append(element)
            else:
                if len(element["seq"]) > len([x for x in gene_spec_elements if x["gene_spec_id"] == gene_spec_id][0]["seq"]):
                    for el in gene_spec_elements:
                        if el["gene_spec_id"] == gene_spec_id:
                            el["seq"] = element["seq"]
                            break

            gene_spec_ids.add(element["gene_spec_id"])
        gene_record_dict[record] = gene_spec_elements

    if (not args.skip_blast):
        time_tracker.print_time_diff("starting blast stuff")
        cmd = f"mkdir -p {work_dir}/protein_extracted_longest"
        subprocess.run(cmd, shell=True)
        for species in species_list:
            print(f"     {species}")
            file = f"{work_dir}/protein_extracted_longest/{species}.faa"


            species_records = {}
            for gene in gene_record_dict:
                species_records[gene] = []
                for element in gene_record_dict[gene]:
                    if element["species"] == species:
                        species_records[gene].append(element)


            with open(file, "w") as outfile:
                for entry in species_records:
                    element = species_records[entry]
                    if len(element) > 1:
                        print(species)
                        assert False
                    if len(element) == 0:
                        continue
                    element = element[0]
                    outfile.write(f">{element['id']} gene:{element['gene_id']} mrna:{element['mrna_id']} {element['species']} {element['name']}\n")
                    outfile.write(f"{element['seq']}\n")

        #make blast db
        cmd = f"mkdir -p {work_dir}/protein_extracted_longest/blastdb"
        subprocess.run(cmd, shell=True, check=True)
        #concat all fasta files
        cmd = f"cat {work_dir}/protein_extracted_longest/*.faa \
            > {work_dir}/protein_extracted_longest/blastdb/all_proteins.faa"
        # run command
        subprocess.run(cmd, shell=True, check=True)

        #make blast db
        cmd = f"makeblastdb -in {work_dir}/protein_extracted_longest/blastdb/all_proteins.faa \
            -dbtype prot\
            -out {work_dir}/protein_extracted_longest/blastdb/all_proteins"
        subprocess.run(cmd, shell=True, check=True)

        # blast all against DB
        time_tracker.print_time_diff("run blast")
        for species in species_list:
            print(f"     {species}")
            cmd = f"blastp -num_threads {cores} -query {work_dir}/protein_extracted_longest/{species}.faa \
            -db {work_dir}/protein_extracted_longest/blastdb/all_proteins \
            -out {work_dir}/protein_extracted_longest/{species}.blastp.tsv \
            -outfmt \"6 qseqid sseqid pident mismatch gapopen gaps qcovs qcovshsp evalue\" "
            subprocess.run(cmd, shell=True, check=True)

    #collect edges for each protein
    time_tracker.print_time_diff("extract edges and save identities")
    edges = {}
    identity_thres = 50
    coverage_thres = 35
    score_thres = 0.001 #this is too high at the moment
    identity_dictionaries = {}
    min_identity_dictionaries = {}

    for species in species_list:
        pair_wise_identities = {} # A->B and B->A identities separately (max for each)
        min_pair_wise_identities = {} # if A->B and B->A then lower of the two identities
                                      # this is useful during reporting, anytime I have genes A and B
                                      # I can get the correct identity to report without having to check
                                      # both directions
        print(f"     {species}")
        file = f"{work_dir}/protein_extracted_longest/{species}.blastp.tsv"
        with open(file, "r") as infile:
            for line in infile:
                line = line.strip().split("\t")

                gene_A = line[0]
                gene_B = line[1]
                identity = float(line[2])
                coverage = int(line[6])
                score = float(line[7])

                if identity >= identity_thres and coverage >= coverage_thres and score < score_thres:
                    if gene_A in edges:
                        edges[gene_A].append(gene_B)
                    else:
                        edges[gene_A] = [gene_B]
                    
                    key_pair_AB = tuple([gene_A,gene_B])

                    if key_pair_AB in pair_wise_identities.keys():
                        if identity > pair_wise_identities[key_pair_AB]:
                            pair_wise_identities[key_pair_AB] = identity
                    else:
                        pair_wise_identities[key_pair_AB] = identity

                    key_pair_BA = tuple([gene_B,gene_A])
                    if key_pair_BA in pair_wise_identities:
                        key_pair = tuple(sorted([gene_A,gene_B]))
                        min_pair_wise_identities[key_pair] = min(pair_wise_identities[key_pair_AB],pair_wise_identities[key_pair_BA])

        identity_dictionaries[species] = pair_wise_identities
        min_identity_dictionaries[species] = min_pair_wise_identities

    #only keep two way edges
    edges_2way = set()
    vertices = set()
    
    time_tracker.print_time_diff("extract two way edges")
    for node_A in edges:
        for node_B in edges[node_A]:
            if node_B in edges and node_A in edges[node_B]:
                tpl = (node_A, node_B)
                tpl = tuple(sorted(tpl))

                edges_2way.add(tpl)
                vertices.add(node_A)
                vertices.add(node_B)   

    #create clusters with transitive clusering
    time_tracker.print_time_diff("create clusters")
    clusters = []
    vertex_to_cluster = {}
    print(f"we have {len(edges_2way)} 2way edges")
    for edge in edges_2way:
        v1, v2 = edge
        # If both vertices are not in any cluster yet
        if v1 not in vertex_to_cluster and v2 not in vertex_to_cluster:
            new_cluster = set(edge)
            clusters.append(new_cluster)
            vertex_to_cluster[v1] = new_cluster
            vertex_to_cluster[v2] = new_cluster
        # If only one vertex is in a cluster
        elif not v1 in vertex_to_cluster:
            cluster = vertex_to_cluster[v2]
            cluster.add(v1)
            vertex_to_cluster[v1] = cluster

        elif not v2 in vertex_to_cluster:
            cluster = vertex_to_cluster[v1]
            cluster.add(v2)
            vertex_to_cluster[v2] = cluster
        
        # If both vertices are in different clusters
        else:
            cluster1 = vertex_to_cluster[v1]
            cluster2 = vertex_to_cluster[v2]
            if cluster1 != cluster2:
                merged_cluster = cluster1.union(cluster2)
                clusters.remove(cluster1)
                clusters.remove(cluster2)
                clusters.append(merged_cluster)
                for vertex in merged_cluster:
                    vertex_to_cluster[vertex] = merged_cluster

    #merge clusters with common elements
    # time_tracker.print_time_diff("merge clusters")
    merged_clusters = clusters

    time_tracker.print_time_diff("assign species")
    prot_to_species = {}
    prot_record_dict = {}
    for gene in gene_record_dict:
        for element in gene_record_dict[gene]:
            my_id = element["id"]
            prot_record_dict[my_id] = element
    
    for cluster in merged_clusters:
        for prot in cluster:
            species = prot_record_dict[prot]['species']
            prot_to_species[prot] = species

    cmd = f"mkdir -p {work_dir}/protein_extracted_longest/clusters_merged"
    subprocess.run(cmd, shell=True, check=True)

    file_clusters = open(f"{work_dir}/protein_extracted_longest/clusters_merged/clusters.tsv", "w")
    file_identities = open(f"{work_dir}/protein_extracted_longest/clusters_merged/clusters_identities.tsv", "w")

    header = "id\tgene_symbol\tDescription(s)"
    species_header = ""
    for species in species_list:
        species_header += f"\t{species}"

    print(f"{header}\tampliconic{species_header}", file=file_clusters)
    print(f"{header}{species_header}", file=file_identities)

    total_len = len(merged_clusters)

    time_tracker.print_time_diff("report clusters")
    
    for counter in range(total_len):

        cluster = clusters[counter]
        results = report(merged_clusters,prot_record_dict, counter, total_len, prot_to_species, min_identity_dictionaries, species_list, work_dir)

        if (results['report']):
            print(results['file_line'], file = file_clusters)
            print(results["identity_line"], file = file_identities)

    file_clusters.close()
    file_identities.close()

    time_tracker.print_time_diff("FINISHED")


if __name__ == "__main__":
    main()
