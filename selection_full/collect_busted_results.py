import os
import json

directories = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]

results = []

for repo in directories:
    json_path = os.path.join(repo, "BUSTED_best.json")
    if not os.path.exists(json_path):
        continue
        
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # Fits
        fit = data.get("fits", {}).get("Unconstrained model", {})
        aic = fit.get("AIC-c", "N/A")
        
        # Rate Distributions
        rate_dist = fit.get("Rate Distributions", {})
        
        # N (Rates) from "Test"
        test_dist = rate_dist.get("Test", {})
        # Count keys
        n_total = len(test_dist.keys())
        # Assuming error sink is always 1 class, Biological N = n_total - 1
        # But let's verify if error sink is actually present.
        
        # Anomaly Check & N correction
        anomaly_score = 0.0
        for k, v in test_dist.items():
            omega = v.get("omega", 0)
            prop = v.get("proportion", 0)
            if omega >= 99.0:
                anomaly_score += prop
        
        # We can report N as the parameter N passed to BUSTED.
        # Since we ran with --error-sink Yes, we expect N+1 classes.
        # So N = n_total - 1.
        N = n_total - 1 if n_total > 0 else 0
            
        # M (Syn-Rates) from "Synonymous site-to-site rates"
        syn_dist = rate_dist.get("Synonymous site-to-site rates", {})
        M = len(syn_dist.keys())
        if M == 0:
            M = 1 # Default constant rate

        # P-value
        test_results = data.get("test results", {})
        p_val = test_results.get("p-value", "N/A")
        if isinstance(p_val, float):
            p_str = f"{p_val:.4g}"
        else:
            p_str = str(p_val)

        results.append({
            "Gene": repo,
            "Rates (N)": N,
            "Syn-Rates (M)": M,
            "c-AIC": f"{aic:.2f}" if isinstance(aic, float) else aic,
            "Anomaly Score": f"{anomaly_score:.4f}",
            "Anomaly?": "YES" if anomaly_score > 0.001 else "No",
            "p-value": p_str
        })
        
    except Exception as e:
        print(f"Error parsing {repo}: {e}")

print("| Gene | Rates (N) | Syn-Rates (M) | c-AIC | Anomaly Score | Anomaly? | p-value |")
print("|------|-----------|---------------|-------|---------------|----------|---------|")
for r in results:
    print(f"| {r['Gene']} | {r['Rates (N)']} | {r['Syn-Rates (M)']} | {r['c-AIC']} | {r['Anomaly Score']} | {r['Anomaly?']} | {r['p-value']} |")