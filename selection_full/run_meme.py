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

def run_meme(alignment, output):
    cmd = [
        "hyphy", "meme",
        "--alignment", alignment,
        "--output", output
    ]
    print(f"Running MEME: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FAILED: {alignment}")
        # print(result.stderr)
        return False
    return True

def count_selected_sites(json_path):
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # MEME results structure: "MLE" -> "content" -> "0" (headers) ...
        # We need to find the column for p-value.
        # Header usually includes: "alpha", "beta", "p-value", etc.
        
        mle = data.get("MLE", {})
        headers = mle.get("headers", [])
        content = mle.get("content", [])
        
        # Content is usually a dict keyed by partition index "0", "1", ...
        # We want to aggregate all sites from all partitions.
        all_rows = []
        if isinstance(content, dict):
            for part_key in content:
                all_rows.extend(content[part_key])
        elif isinstance(content, list):
            all_rows = content
            
        try:
            pval_idx = -1
            for i, h in enumerate(headers):
                # Header can be ["ColName", "Desc"] or just "ColName"
                col_name = ""
                if isinstance(h, list) and len(h) > 0:
                    col_name = h[0]
                elif isinstance(h, str):
                    col_name = h
                
                if "p-value" in col_name:
                    pval_idx = i
                    break
            
            # Count sites with p <= 0.05
            count = 0
            for row in all_rows:
                if pval_idx != -1 and len(row) > pval_idx:
                    p = row[pval_idx]
                    # HyPhy sometimes puts null or errors? Assume float.
                    if isinstance(p, (int, float)) and p <= 0.05:
                        count += 1
                elif len(row) > 6:
                     # Fallback to index 6 if standard
                     p = row[6]
                     if isinstance(p, (int, float)) and p <= 0.05:
                        count += 1
                        
            return count, len(all_rows)
            
        except Exception as e:
            print(f"Error processing content: {e}")
            return 0, 0

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
    
    output_json = os.path.join(repo_path, "MEME.json")
    
    if not os.path.exists(output_json):
        if not run_meme(alignment, output_json):
            continue
    
    selected, total = count_selected_sites(output_json)
    results.append({
        "Gene": repo,
        "Selected Sites": selected,
        "Total Sites": total
    })

print("\n### MEME Results (p < 0.1)\n")
print("| Gene | Selected Sites | Total Sites |")
print("|------|----------------|-------------|")
for r in results:
    print(f"| {r['Gene']} | {r['Selected Sites']} | {r['Total Sites']} |")
