import sys
from Bio import SeqIO

def clean_stops(input_fasta, output_fasta):
    cleaned_records = []
    # Stop codons
    stops = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}
    
    # Read manually to preserve exact headers if Biopython messes them up?
    # Biopython usually does a good job, but let's be careful.
    # Actually, simple string processing is safer for preserving headers exactly.
    
    with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
        header = None
        seq = []
        
        for line in f_in:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header:
                    # Process previous sequence
                    full_seq = "".join(seq)
                    new_seq = []
                    for k in range(0, len(full_seq), 3):
                        codon = full_seq[k:k+3]
                        if len(codon) == 3 and codon in stops:
                            new_seq.append("---")
                        else:
                            new_seq.append(codon)
                    f_out.write(f"{header}\n{''.join(new_seq)}\n")
                
                header = line
                seq = []
            else:
                seq.append(line)
        
        # Process last sequence
        if header:
            full_seq = "".join(seq)
            new_seq = []
            for k in range(0, len(full_seq), 3):
                codon = full_seq[k:k+3]
                if len(codon) == 3 and codon in stops:
                    new_seq.append("---")
                else:
                    new_seq.append(codon)
            f_out.write(f"{header}\n{''.join(new_seq)}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 clean_stops.py <input_fasta> <output_fasta>")
        sys.exit(1)
    
    clean_stops(sys.argv[1], sys.argv[2])
