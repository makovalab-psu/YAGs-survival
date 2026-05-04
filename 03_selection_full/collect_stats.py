import os
import glob
import re

directories = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]

stats = []

for repo in directories:
    repo_path = repo
    if not os.path.exists(repo_path):
        continue
    
    # Original sequences
    fasta_files = [f for f in (glob.glob(os.path.join(repo_path, "*.fasta")) + glob.glob(os.path.join(repo_path, "*.Y-only"))) if ".unique." not in f]
    if not fasta_files:
        continue
    
    orig_file = fasta_files[0]
    with open(orig_file, 'r') as f:
        orig_count = sum(1 for line in f if line.startswith('>'))

    # Unique sequences and codons from Nexus
    nexus_file = os.path.join(repo_path, f"{repo if 'DAZ' not in repo or 'RRM' not in repo or 'repeats' not in repo else repo.replace('_RRM','.unique').replace('_repeats','.unique')}.unique.nexus")
    # Wait, my previous naming logic in run_deduplicate.py was:
    # output_file = os.path.join(repo_path, f"{repo}.unique.nexus")
    # EXCEPT for DAZ cases where I might have messed up the name slightly in the previous script? 
    # Let me check what I actually created.
    
    # Standardized check for what I actually outputted:
    nexus_files = glob.glob(os.path.join(repo_path, "*.unique.nexus"))
    if not nexus_files:
        # Fallback to check if I used a different extension or name
        nexus_files = glob.glob(os.path.join(repo_path, "*.unique.fasta"))
        if not nexus_files:
            continue
            
    nexus_file = nexus_files[0]
    
    ntax = 0
    nchar = 0
    with open(nexus_file, 'r') as f:
        content = f.read()
        ntax_match = re.search(r"NTAX\s*=\s*(\d+)", content, re.IGNORECASE)
        nchar_match = re.search(r"NCHAR\s*=\s*(\d+)", content, re.IGNORECASE)
        if ntax_match:
            ntax = int(ntax_match.group(1))
        if nchar_match:
            nchar = int(nchar_match.group(1))

    stats.append({
        "Gene": repo,
        "Original Seqs": orig_count,
        "Unique Seqs": ntax,
        "Codons": nchar // 3
    })

# Print Markdown Table
print("| Gene | Original Seqs | Unique Seqs | Codons |")
print("|------|---------------|-------------|--------|")
for s in stats:
    print(f"| {s['Gene']} | {s['Original Seqs']} | {s['Unique Seqs']} | {s['Codons']} |")
