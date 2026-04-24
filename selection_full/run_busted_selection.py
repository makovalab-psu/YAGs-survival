import os
import subprocess
import glob
import json
import shutil

# Configuration
DIRECTORIES = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]
BASE_DIR = "/Users/sergei/Dropbox/Work/Collaborations/Y-chromosome/2026/round6_export_20260126"

def get_aic_anomaly_pvalue(json_path):
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # Get c-AIC
        fit = data.get("fits", {}).get("Unconstrained model", {})
        aic = fit.get("AIC-c")
        
        # Check for anomaly
        rate_dist = fit.get("Rate Distributions", {})
        omega_dist = rate_dist.get("non-synonymous/synonymous rate ratio")
        
        anomaly_score = 0.0
        if isinstance(omega_dist, list):
             for entry in omega_dist:
                 if len(entry) >= 2:
                     rate = entry[0]
                     prop = entry[1]
                     if rate >= 99.0:
                         anomaly_score += prop

        # Get p-value
        test_results = data.get("test results", {})
        p_value = test_results.get("p-value", "N/A")
        
        return aic, anomaly_score, p_value
    except Exception as e:
        return None, 0.0, "N/A"

# ... (omitted parts same as before) ...

    # Run initial model
    temp_json = os.path.join(repo_path, f"BUSTED_temp.json")
    final_json = os.path.join(repo_path, "BUSTED_best.json")
    
    # We rely on previous run existing, or re-run?
    # Since we want to extract p-value from the BEST model found, we should read BUSTED_best.json
    
    if os.path.exists(final_json):
        aic, anomaly, p_val = get_aic_anomaly_pvalue(final_json)
        
        # Format p-value
        if isinstance(p_val, float):
            p_str = f"{p_val:.4g}"
        else:
            p_str = str(p_val)
            
        results_table.append({
            "Gene": repo,
            # We need to re-find the best N/M or just assume the previous run worked?
            # Actually, I should just write a small script to parse the EXISTING best json files.
            # No need to re-run the whole selection.
            "p-value": p_str
        })

# Instead of re-running the whole selection, let's make a new script "collect_busted_pvalues.py"

def run_busted(alignment, output, rates, syn_rates, multiple_hits="None"):

    cmd = [

        "hyphy", "busted",

        "--alignment", alignment,

        "--rates", str(rates),

        "--error-sink", "Yes",

        "--starting-points", "1",

        "--multiple-hits", multiple_hits,

        "--output", output

    ]

    

    if syn_rates == 1:

        cmd.extend(["--srv", "No"])

    else:

        cmd.extend(["--srv", "Yes", "--syn-rates", str(syn_rates)])

        

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:

        # print(f"FAILED: {' '.join(cmd)}")

        return False

    return True



results_table = []



for repo in DIRECTORIES:

    repo_path = os.path.join(BASE_DIR, repo)

    if not os.path.exists(repo_path):

        continue

    

    # Find alignment

    nexus_files = glob.glob(os.path.join(repo_path, "*.unique.nexus"))

    if not nexus_files:

        continue

    alignment = nexus_files[0]

    

    print(f"\nProcessing {repo}...")

    

    # Hill Climbing starting at N=1, M=1

    current_n = 1

    current_m = 1

    

    best_aic = float('inf')

    best_model = (1, 1)

    best_anomaly = 0.0

    best_mh = "None"

    

    temp_json = os.path.join(repo_path, f"BUSTED_temp.json")

    final_json = os.path.join(repo_path, "BUSTED_best.json")

    

    # Initial Run

    if run_busted(alignment, temp_json, current_n, current_m):

        aic, anomaly, _ = get_aic_anomaly_pvalue(temp_json)

        if aic is not None:

            best_aic = aic

            best_anomaly = anomaly

            shutil.copy(temp_json, final_json)

            print(f"  Initial (N={current_n}, M={current_m}): AIC-c={aic:.2f}")

        else:

            print("  Failed to extract AIC from initial run.")

            continue

    else:

        print("  Initial run failed.")

        continue

    

    # Rate Class Optimization

    while True:

        candidates = [

            (current_n + 1, current_m),

            (current_n, current_m + 1)

        ]

        

        improved = False

        step_best_aic = best_aic

        step_best_model = best_model

        step_best_anomaly = best_anomaly

        

        for (n, m) in candidates:

            # Safety limit

            if n > 6 or m > 6:

                continue

                

            print(f"  Testing (N={n}, M={m})...")

            if run_busted(alignment, temp_json, n, m):

                aic, anomaly, _ = get_aic_anomaly_pvalue(temp_json)

                if aic is not None:

                    print(f"    AIC-c={aic:.2f}")

                    if aic < step_best_aic - 0.01: # Use a small epsilon for improvement

                        step_best_aic = aic

                        step_best_model = (n, m)

                        step_best_anomaly = anomaly

                        shutil.copy(temp_json, final_json)

                        improved = True

            else:

                print("    Run failed.")

        

        if improved:

            print(f"  Improved! Moving to {step_best_model} with AIC-c={step_best_aic:.2f}")

            current_n, current_m = step_best_model

            best_aic = step_best_aic

            best_anomaly = step_best_anomaly

            best_model = step_best_model

        else:

            print(f"  No rate class improvement. Stopping at {best_model}.")

            break

            

    # Multiple Hits Optimization

    print(f"  Testing Multiple Hits (Double+Triple) on top of {best_model}...")

    if run_busted(alignment, temp_json, best_model[0], best_model[1], multiple_hits="Double+Triple"):

        aic, anomaly, _ = get_aic_anomaly_pvalue(temp_json)

        if aic is not None:

            print(f"    MH AIC-c={aic:.2f}")

            if aic < best_aic - 0.01:

                print(f"  Multiple Hits improved AIC! ({best_aic:.2f} -> {aic:.2f})")

                best_aic = aic

                best_anomaly = anomaly

                best_mh = "Double+Triple"

                shutil.copy(temp_json, final_json)

            else:

                print("  Multiple Hits did not improve AIC.")

        else:

            print("    Failed to extract AIC from MH run.")

    else:

        print("    MH Run failed.")



    # Get final p-value from best json

    _, _, best_p = get_aic_anomaly_pvalue(final_json)

    if isinstance(best_p, float):

        best_p = f"{best_p:.4g}"



    # Record results

    results_table.append({

        "Gene": repo,

        "Rates (N)": best_model[0],

        "Syn-Rates (M)": best_model[1],

        "MH": best_mh,

        "c-AIC": f"{best_aic:.2f}",

        "Anomaly Score": f"{best_anomaly:.4f}",

        "Anomaly Detected": "YES" if best_anomaly > 0.001 else "No",

        "p-value": str(best_p)

    })

    

    if os.path.exists(temp_json):

        os.remove(temp_json)



# Print Summary

print("\n### BUSTED Model Selection Results\n")

print("| Gene | Rates (N) | Syn-Rates (M) | MH | c-AIC | Anomaly Score | Anomaly? | p-value |")

print("|------|-----------|---------------|----|-------|---------------|----------|---------|")

for r in results_table:

    print(f"| {r['Gene']} | {r['Rates (N)']} | {r['Syn-Rates (M)']} | {r['MH']} | {r['c-AIC']} | {r['Anomaly Score']} | {r['Anomaly Detected']} | {r['p-value']} |")
