import os
import subprocess
import glob
import json

# Configuration
DIRECTORIES = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]
BASE_DIR = "/Users/sergei/Dropbox/Work/Collaborations/Y-chromosome/2026/round6_export_20260126"

def run_absrel(alignment, output):
    cmd = [
        "hyphy", "absrel",
        "--alignment", alignment,
        "--branches", "All",
        "--output", output
    ]
    print(f"Running aBSREL: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FAILED: {alignment}")
        return False
    return True

def get_selected_branches(json_path):
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        test_results = data.get("test results", {})
        positives = test_results.get("positive test results", 0)
        tested = test_results.get("tested", 0)
        
        return positives, tested
    except Exception as e:
        print(f"Error parsing {json_path}: {e}")
        return 0, 0

results = []

for repo in DIRECTORIES:
    repo_path = os.path.join(BASE_DIR, repo)
    if not os.path.exists(repo_path):
        continue
    
    nexus_files = glob.glob(os.path.join(repo_path, "*.unique.nexus"))
    if not nexus_files:
        continue
    alignment = nexus_files[0]
    
    output_json = os.path.join(repo_path, "ABSREL.json")
    
    if run_absrel(alignment, output_json):
        positives, tested = get_selected_branches(output_json)
        results.append({
            "Gene": repo,
            "Selected Branches": positives,
            "Tested Branches": tested
        })

print("\n### aBSREL Results\n")
print("| Gene | Selected Branches | Tested Branches |")
print("|------|-------------------|-----------------|")
for r in results:
    print(f"| {r['Gene']} | {r['Selected Branches']} | {r['Tested Branches']} |")
