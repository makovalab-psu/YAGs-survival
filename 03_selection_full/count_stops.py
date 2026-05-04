import os
import glob

directories = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]

stops = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}

print("| Gene | Internal Stop Codons | Sequences with Internal Stops |")
print("|------|-------------------|----------------------|")

for repo in directories:
    if not os.path.exists(repo):
        continue
    
    # Find original fasta file
    fasta_files = [f for f in (glob.glob(os.path.join(repo, "*.fasta")) + glob.glob(os.path.join(repo, "*.Y-only"))) if ".unique." not in f and ".cln." not in f and ".clean." not in f]
    
    if not fasta_files:
        continue
    
    input_fasta = fasta_files[0]
    
    total_stops = 0
    seqs_with_stops = 0
    
    with open(input_fasta, 'r') as f:
        seq = []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if seq:
                    full_seq = "".join(seq)
                    # Check for INTERNAL stops (all except the last codon)
                    seq_stops = 0
                    for k in range(0, len(full_seq) - 3, 3):
                        codon = full_seq[k:k+3]
                        if len(codon) == 3 and codon in stops:
                            seq_stops += 1
                    
                    total_stops += seq_stops
                    if seq_stops > 0:
                        seqs_with_stops += 1
                seq = []
            else:
                seq.append(line)
        
        # Last sequence
        if seq:
            full_seq = "".join(seq)
            seq_stops = 0
            for k in range(0, len(full_seq) - 3, 3):
                codon = full_seq[k:k+3]
                if len(codon) == 3 and codon in stops:
                    seq_stops += 1
            
            total_stops += seq_stops
            if seq_stops > 0:
                seqs_with_stops += 1

    print(f"| {repo} | {total_stops} | {seqs_with_stops} |")