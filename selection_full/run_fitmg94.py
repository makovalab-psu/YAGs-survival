import os
import subprocess
import glob

# Configuration
HYPHY_SCRIPT = "/Users/sergei/Development/hyphy-analyses/FitMG94/FitMG94.bf"
BASE_DIR = "/Users/sergei/Dropbox/Work/Collaborations/Y-chromosome/2026/round6_export_20260126"

directories = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]

def run_fitmg94(alignment_path, output_path):
    cmd = [
        "hyphy", HYPHY_SCRIPT,
        "--alignment", os.path.abspath(alignment_path),
        "--lrt", "Yes",
        "--output", os.path.abspath(output_path)
    ]
    print(f"Running FitMG94: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running FitMG94 on {alignment_path}")
        print(result.stderr)
        print(result.stdout)
    else:
        print(f"Successfully processed {alignment_path}")

for repo in directories:
    repo_path = os.path.join(BASE_DIR, repo)
    if not os.path.exists(repo_path):
        continue
    
    alignment_file = os.path.join(repo_path, f"{repo}.unique.nexus")
    if not os.path.exists(alignment_file):
        # Handle cases where the gene name in directory might differ from the file name prefix
        nexus_files = glob.glob(os.path.join(repo_path, "*.unique.nexus"))
        if not nexus_files:
            print(f"No unique nexus file found in {repo_path}, skipping.")
            continue
        alignment_file = nexus_files[0]
        
    output_json = os.path.join(repo_path, "MG94fit.json")
    
    run_fitmg94(alignment_file, output_json)
