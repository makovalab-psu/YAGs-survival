from BCBio import GFF
from Bio import SeqIO
import gspread


def main():
    gc = gspread.oauth()
    sheet = "Y chromosome - various tables"
    sh = gc.open(sheet)
    worksheet = sh.worksheet("Multicopy genes accessions")
    data = worksheet.get_all_values()

    gene_list = []

    # print("loading gene list")

    for i in range(len(data)):
        
        row = data[i]
        if row[0] != "HomSap":
            continue
        if row[4].startswith("x"):
            continue

        gene_list.append(f"{row[1]}_{row[3]}")

        

    # print("gene list loades")

    # print(gene_list)

    # print("loading gff file")

    rna_id_dict = {}
    rna_id_list = []
    rna_id_to_gene = {}

    GFF_FILE = "/Users/kxp5629/proj/Y/data/references/HomSap/ncbi_dataset/data/GCF_009914755.1/genomic_chrY.gff"
    with open(GFF_FILE) as in_handle:
        for rec in GFF.parse(in_handle):
            # print(rec)
            # print("")
            for feature in rec.features:
                if(feature.type == 'gene'):
                    for sub_feature in feature.sub_features:
                        if(sub_feature.type == 'mRNA'):
                            # print(sub_feature)
                            # print(sub_feature.qualifiers['Parent'])
                            if any(gene.split("_")[1] == sub_feature.qualifiers['gene'][0] for gene in gene_list):
                                # print("found") 
                                gene_name = sub_feature.qualifiers['Parent'][0]

        
                                rna_id_list.append(sub_feature.qualifiers['ID'][0].split("-")[1])
                                rna_id_to_gene[sub_feature.qualifiers['ID'][0].split("-")[1]] = gene_name
                                if gene_name in rna_id_dict:
                                    rna_id_dict[gene_name].append(sub_feature.qualifiers['ID'][0].split("-")[1])
                                else:
                                    rna_id_dict[gene_name] = [sub_feature.qualifiers['ID'][0].split("-")[1]]
                                
                            
                            
    print(rna_id_dict)

    sequences = {}

    fasta_file = "/Users/kxp5629/proj/Y/data/references/HomSap/ncbi_dataset/data/GCF_009914755.1/rna.fna"
    
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in rna_id_list:
                sequence = record.seq

                gene_name = rna_id_to_gene[record.id]
                if gene_name in sequences:
                    if len(sequences[gene_name]) < len(sequence):
                        sequences[gene_name] = str(sequence)
                else:
                    sequences[gene_name] = str(sequence)

    with open ("human_YAG_sequences.fasta", "w") as handle:

        
        for gene in sequences:
            gene_family = ""
            # which gene is this rna associated with
            for g in gene_list:
                if g.split("_")[1] == gene.split("-")[1]:
                    gene_family = g.split("_")[0]
                    break

            handle.write(f">{gene_family}_{gene}\n")
            handle.write(f"{sequences[gene]}\n")

if __name__ == "__main__":
    main()

