from Bio import SeqIO
from Bio.Seq import Seq
import os
import subprocess

def main():
    header ="sample\tgene_family\tsequence\tprotein\tstrand\tdna_start\tdna_end\tbed12_lengths\tbed12_starts"

    fasta_folder = "/Users/kxp5629/proj/Y/data/43_human_Y/assemblies/verkko_1.0_release/chrY/"

    all_out_file = "alignments/all_genes.txt"

    if not os.path.exists(f"alignments"):
                os.makedirs(f"alignments")

    with open(all_out_file, "w") as all_handle:
        print(header, file=all_handle)

    for fasta_file in os.listdir(fasta_folder):
        if fasta_file.endswith(".fasta"):
            sample = fasta_file.split(".")[0]

            

            file_name = ".".join(fasta_file.split(".")[:-1])

            if not os.path.exists(f"alignments/{file_name}_alignment.sam"):

                output = f"alignments/{file_name}.bed"

                subprocess.run(["minimap2", "-x", "splice", "-t", "8", f"{fasta_folder}{fasta_file}", "human_YAG_sequences.fasta","-a", "-o", f"alignments/{file_name}_alignment.sam"])
                subprocess.run(["bedtools", "bamtobed","-bed12", "-i", f"alignments/{file_name}_alignment.sam"], stdout=open(output, "w"))


                filtered_file = f"alignments/{file_name}_filtered.bed"
                sorted_file = f"alignments/{file_name}_filtered.sorted.bed"
                merged_file = f"alignments/{file_name}_filtered.sorted.merged.bed"
                with open(output, "r") as handle, \
                    open(filtered_file, "w") as filtered_handle, \
                    open(sorted_file, "w") as sorted_handle, \
                    open(merged_file, "w") as merged_handle:
                    for line in handle:
                        line = line.strip()
                        fields = line.split("\t")

                        gene = fields[3].split("_")[0]
                        if (int(fields[2]) - int(fields[1]) > 25000 ) and gene == "BPY2":
                            continue

                        if (int(fields[2]) - int(fields[1]) > 3000 ) and gene == "CDY":
                            continue

                        if (int(fields[2]) - int(fields[1]) > 90000 ) and gene == "DAZ":
                            continue   

                        if (int(fields[2]) - int(fields[1]) > 3500 ) and gene == "HSFY":
                            continue  

                        if (int(fields[2]) - int(fields[1]) > 24000 ) and gene == "RBMY":
                            continue

                        if (int(fields[2]) - int(fields[1]) > 5000 ) and gene == "TSPY":
                            continue

                        if (int(fields[2]) - int(fields[1]) > 2000 ) and gene == "VCY":
                            continue

                        print(line, file=filtered_handle)

                    subprocess.run(["sort", "-k1,1", "-k2,2n", filtered_file], stdout=sorted_handle)
                    subprocess.run(["bedtools", "merge", "-i", sorted_file],   stdout=merged_handle)

                #for each region in the merged bed file, get regions form bed12 file extract the sequences from the reference 
                # translate and make sure there is no stop codon in at least the first half of the sequence
                # collect "good" sequences, infer gene name

                subdir = f"alignments/{file_name}_sub_dir"
                if not os.path.exists(subdir):
                    os.makedirs(subdir)

                with open(f"alignments/{file_name}_genes.txt", "w") as gene_handle:
                    print(header, file=gene_handle)
                with open(merged_file, "r") as handle:
                    i = 0
                    
                    for line in handle:
                        i += 1
                        sub_fasta = f"{subdir}/region_{i}.fasta"
                        intersect_bed = f"{subdir}/region_{i}_intersect.bed"
                        with open (f"{subdir}/region_{i}.bed", "w") as region_handle:
                            print(line, file=region_handle)
                        
                        #intersect bed12

                        try:
                            subprocess.run(["bedtools", "intersect", "-a", f"alignments/{file_name}.bed", "-b", f"{subdir}/region_{i}.bed", "-split"], stdout=open(intersect_bed, "w"))
                        except subprocess.CalledProcessError as e:
                            print(e.stderr)
                            assert False

                        #extract sequences
                        try:
                            subprocess.run(["bedtools", "getfasta", "-split", "-s", "-fi", f"{fasta_folder}{fasta_file}", "-bed", intersect_bed, "-fo", sub_fasta])
                        except subprocess.CalledProcessError as e:
                            print(e.stderr)
                            assert False

                
                        with open(sub_fasta, 'r') as handle, open(intersect_bed, 'r') as bed_handle:
                            #TODO: check if bed and fasta coordinates match
                            all_orfs = []
                            bed_lines = bed_handle.readlines()
                            gene_family = ""
                            for record in SeqIO.parse(handle, "fasta"):
                                bed = bed_lines.pop(0).strip().split("\t")
                                gene_family = bed[3].split("_")[0]
                                for frame in range(3):
                                    prot = record.seq[frame:].translate(table=1)

                                    start_positions = [i for i, aa in enumerate(prot) if aa == 'M']
                                    

                                    for start in start_positions:
                                        current_prot = prot[start:]
                                        if '*' in current_prot:
                                            stop = current_prot.index('*')
                                            current_prot = current_prot[:stop]
                                            orf = current_prot[:stop]
                                            if len(current_prot) > 50:
                                                dna_start = frame + (start * 3)
                                                dna_end = dna_start + (len(current_prot) * 3) + 3

                                                all_orfs.append({
                                                    'frame': frame + 1,
                                                    'protein': str(orf),
                                                    'length': len(orf),
                                                    'dna_start': dna_start + 1,  # 1-based coordinates
                                                    'dna_end': dna_end,
                                                    'gene_family': bed[3].split("_")[0],
                                                    'sequence': record.seq,
                                                    'bed_start': int(bed[1]),
                                                    'bed_end': int(bed[2]),
                                                    'strand': bed[5],
                                                    'bed12_lengths': bed[10],
                                                    'bed12_starts': bed[11],
                                                    "chromosome": bed[0]
                                                })
                            all_orfs.sort(key=lambda x: x['length'], reverse=True)

                            # if gene_family=="HSFY":
                            #     print("FOUND HSFY")
                            #     print(all_orfs)
                            if len(all_orfs) == 0:
                                assert False
                            [print(f"{sample}\t{orf['gene_family']}\t{orf['sequence']}\t{orf['protein']}\t{orf['chromosome']}\t{orf['strand']}\t{orf['bed_start']}\t{orf['bed_end']}\t{orf['bed12_lengths']}\t{orf['bed12_starts']}", file=open(f"alignments/{file_name}_genes.txt", "a")) for orf in all_orfs[:1]]
                            [print(f"{sample}\t{orf['gene_family']}\t{orf['sequence']}\t{orf['protein']}\t{orf['chromosome']}\t{orf['strand']}\t{orf['bed_start']}\t{orf['bed_end']}\t{orf['bed12_lengths']}\t{orf['bed12_starts']}", file=open(all_out_file, "a")) for orf in all_orfs[:1]]
                           
                            # print(f"{gene_family}\t{sequence}\t{prot_seq}\t{strand}\t{start}\t{end}\t{bed12_lengths}\t{bed12_starts}", file=open(f"alignments/{file_name}_genes.txt", "a"))

                 


if __name__ == "__main__":
    main()