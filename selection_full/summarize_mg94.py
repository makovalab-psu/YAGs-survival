import os
import json

directories = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]

results = []

for repo in directories:
    json_path = os.path.join(repo, "MG94fit.json")
    if not os.path.exists(json_path):
        continue
        
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        fit_mg94 = data.get("fits", {}).get("Standard MG94", {})
        
        # dN/dS from Rate Distributions
        rate_dist = fit_mg94.get("Rate Distributions", {})
        omega = rate_dist.get("non-synonymous/synonymous rate ratio", "N/A")
        
        # Confidence Intervals
        ci = fit_mg94.get("Confidence Intervals", {}).get("non-synonymous/synonymous rate ratio", {})
        omega_ci = "N/A"
        if isinstance(ci, dict):
             lb = ci.get("LB")
             ub = ci.get("UB")
             if lb is not None and ub is not None:
                 omega_ci = f"[{lb:.3f}, {ub:.3f}]"

        # P-value from "test results"
        test_res = data.get("test results", {}).get("non-synonymous/synonymous rate ratio", {})
        p_val = test_res.get("p-value", "N/A")
        if isinstance(p_val, float):
            p_str = f"{p_val:.4g}"
        else:
            p_str = str(p_val)

        results.append({
            "Gene": repo,
            "dN/dS": f"{omega:.3f}" if isinstance(omega, (int, float)) else omega,
            "dN/dS CI": omega_ci,
            "p-value (dN/dS != 1)": p_str
        })
        
    except Exception as e:
        print(f"Error parsing {repo}: {e}")

# Generate Table
print("| Gene | dN/dS | 95% CI | p-value (dN/dS != 1) |")
print("|------|-------|--------|----------------------|")
for r in results:
    print(f"| {r['Gene']} | {r['dN/dS']} | {r['dN/dS CI']} | {r['p-value (dN/dS != 1)']} |")
