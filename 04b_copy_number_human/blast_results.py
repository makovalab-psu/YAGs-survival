from Bio import SeqIO
import subprocess
import time
import tempfile
import os


def blast_orfs(input_table, out_file):

    blast_db = "/Users/kxp5629/proj/Y/data/blast/human_protein"

    with open(input_table, "r") as table, open(out_file,'w') as output:

        
        for line in table:
            # print(line)
            if line.startswith("sample"):
                continue
            line = line.strip().split("\t")
            orf = line[3]


            with tempfile.NamedTemporaryFile("w+", delete=False) as temp_seq, \
                tempfile.NamedTemporaryFile("w+", delete=False) as temp_out:

                try:

                    temp_seq.write(f">{line[3]}\n{line[3]}\n")
                    temp_seq.flush()


                    cmd = [ 
                        'blastp',
                        '-query', temp_seq.name,
                        '-db', blast_db,
                        '-outfmt', '6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle',
                        '-max_target_seqs', '1',
                        '-num_threads', '4',
                        '-evalue', '1e-10',
                        '-out', temp_out.name
                    ]
                
                    process = subprocess.run(cmd, capture_output=True, text=True, check=True)


                    with open(temp_out.name) as f:

                        blast_result = f.readline().strip()

                        res = "no blastp hit"
                        if blast_result:
                            res = blast_result.split("\t")[11]

                        print(f"{'\t'.join(line)}\t{res}")

                except subprocess.CalledProcessError as e:
                    print(f"BLAST error: {e.stderr}")
                    return None
                
                finally:
                    os.unlink(temp_seq.name)
                    os.unlink(temp_out.name)
                    


if __name__ == "__main__":
    blast_orfs("/Users/kxp5629/proj/Y/src/18_get_human_genes/alignments_v3/all_genes.txt", "blast_result.tsv")