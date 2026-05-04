import os
import subprocess
import glob
import json
from scipy.stats import chi2

# Configuration
HYPHY_SCRIPT = "/Users/sergei/Development/hyphy-analyses/FitMG94/FitMG94.bf"
BASE_DIR = "/Users/sergei/Dropbox/Work/Collaborations/Y-chromosome/2026/round6_export_20260126"

directories = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]

results = []

def run_mg94_local(alignment_path, output_path):
    cmd = [
        "hyphy", HYPHY_SCRIPT,
        "--alignment", os.path.abspath(alignment_path),
        "--type", "local",
        "--output", os.path.abspath(output_path)
    ]
    print(f"Running FitMG94 Local: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running FitMG94 Local on {alignment_path}")
        # print(result.stderr)
        return False
    return True

def get_fit_stats(json_path, model_name):
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # Structure depends on model name
        # For global: "fits" -> "Global MG94xREV" (usually)
        # For local: "fits" -> "Local MG94xREV" (likely)
        
        fits = data.get("fits", {})
        
        # Try to find the fit by keyword or just take the first one?
        # Global run usually has "Global MG94xREV"
        # Local run usually has "Local MG94xREV"
        
        fit_data = fits.get(model_name)
        if not fit_data:
            # Fallback: look for any key containing "MG94"
            for k in fits.keys():
                if "MG94" in k:
                    fit_data = fits[k]
                    break
        
        if not fit_data:
            return None, None

        log_l = fit_data.get("Log Likelihood")
        params = fit_data.get("estimated parameters")
        return log_l, params
        
    except Exception as e:
        print(f"Error parsing {json_path}: {e}")
        return None, None

print("| Gene | Global LogL | Global Params | Local LogL | Local Params | LRT | df | p-value |")
print("|------|-------------|---------------|------------|--------------|-----|----|---------|")

for repo in directories:
    repo_path = os.path.join(BASE_DIR, repo)
    if not os.path.exists(repo_path):
        continue
    
    # Find alignment
    nexus_files = glob.glob(os.path.join(repo_path, "*.unique.nexus"))
    if not nexus_files:
        continue
    alignment_file = nexus_files[0]
    
    global_json = os.path.join(repo_path, "MG94fit.json")
    local_json = os.path.join(repo_path, "MGlocal.json")
    
    # Run Local Analysis
    if not run_mg94_local(alignment_file, local_json):
        continue
        
    # Get Stats
    # Global model name usually "Global MG94xREV"
    # Local model name usually "Local MG94xREV"
    
    g_logl, g_params = get_fit_stats(global_json, "Global MG94xREV")
    l_logl, l_params = get_fit_stats(local_json, "Local MG94xREV")
    
    if g_logl is not None and l_logl is not None:
        lrt = 2 * (l_logl - g_logl)
        df = l_params - g_params
        
        p_val = 1.0
        if df > 0:
            p_val = chi2.sf(lrt, df)
        elif df == 0:
            p_val = 1.0 # No difference in parameters implies models are effectively same or nested boundary
        else:
            p_val = "N/A (df<0)" # Should not happen if Local nests Global
            
        p_str = f"{p_val:.4g}" if isinstance(p_val, float) else p_val
        
        print(f"| {repo} | {g_logl:.2f} | {g_params} | {l_logl:.2f} | {l_params} | {lrt:.2f} | {df} | {p_str} |")
        
        results.append({
            "Gene": repo,
            "LRT": lrt,
            "df": df,
            "p-value": p_str
        })
    else:
        print(f"| {repo} | N/A | N/A | N/A | N/A | N/A | N/A | N/A |")

