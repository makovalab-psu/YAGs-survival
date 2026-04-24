import os
import subprocess
import glob

# Configuration
HYPHY_SCRIPT = "/Users/sergei/Development/hyphy-analyses/remove-duplicates/remove-duplicates.bf"
BASE_DIR = "/Users/sergei/Dropbox/Work/Collaborations/Y-chromosome/2026/round6_export_20260126"

directories = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]

def run_hyphy_dedup(msa_path, tree_path, output_path):
    cmd = [
        "hyphy", HYPHY_SCRIPT,
        "--msa", os.path.abspath(msa_path),
        "--tree", os.path.abspath(tree_path),
        "--output", os.path.abspath(output_path)
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running HyPhy on {msa_path}")
        print(result.stderr)
        print(result.stdout)
    else:
        print(f"Successfully processed {msa_path}")

def run_python_clean(input_path, output_path):
    cmd = [
        "python3", "clean_stops.py",
        os.path.abspath(input_path),
        os.path.abspath(output_path)
    ]
    print(f"Cleaning (Python): {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error cleaning {input_path}")
        print(result.stderr)
        print(result.stdout)
    else:
        print(f"Successfully cleaned {input_path}")

for repo in directories:
    repo_path = os.path.join(BASE_DIR, repo)
    if not os.path.exists(repo_path):
        print(f"Directory {repo_path} not found, skipping.")
        continue
    
    # Find tree file
    tree_files = glob.glob(os.path.join(repo_path, "*.treefile"))
    if not tree_files:
        print(f"No tree file found in {repo_path}, skipping.")
        continue
    tree_file = tree_files[0]
    
    # Find fasta file
    fasta_files = [f for f in (glob.glob(os.path.join(repo_path, "*.fasta")) + glob.glob(os.path.join(repo_path, "*.Y-only"))) if ".unique." not in f and ".cln." not in f and ".clean." not in f]
    if not fasta_files:
        print(f"No fasta file found in {repo_path}, skipping.")
        continue
    
    msa_file = fasta_files[0]
    
    clean_fasta = os.path.join(repo_path, f"{repo}.clean.fasta")
    final_nexus = os.path.join(repo_path, f"{repo}.unique.nexus")
    
    # Step 1: Clean (remove terminal stops) - using Python script
    run_python_clean(msa_file, clean_fasta)
    
    # Step 2: De-duplicate using the cleaned MSA and the original tree
    if os.path.exists(clean_fasta):
        run_hyphy_dedup(clean_fasta, tree_file, final_nexus)
