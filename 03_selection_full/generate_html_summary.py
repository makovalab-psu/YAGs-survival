import os
import glob
import json
import re
from scipy.stats import chi2

# Configuration
DIRECTORIES = [
    "BPY2", "CDY", "CDY_CDYL", "DAZ_DAZL_repeats", "DAZ_DAZL_RRM", 
    "DAZ_repeats", "DAZ_RRM", "HSFY", "RBMY", "RBMY_RBMX", "TSPY", "VCY"
]
BASE_DIR = "/Users/sergei/Dropbox/Work/Collaborations/Y-chromosome/2026/round6_export_20260126"

def get_alignment_stats(repo_path):
    # Original sequences
    orig_files = [f for f in (glob.glob(os.path.join(repo_path, "*.fasta")) + glob.glob(os.path.join(repo_path, "*.Y-only"))) if ".unique." not in f and ".cln." not in f and ".clean." not in f]
    orig_count = "N/A"
    if orig_files:
        with open(orig_files[0], 'r') as f:
            orig_count = sum(1 for line in f if line.startswith('>'))
            
    # Unique sequences and codons
    nexus_files = glob.glob(os.path.join(repo_path, "*.unique.nexus"))
    unique_count = "N/A"
    codons = "N/A"
    if nexus_files:
        with open(nexus_files[0], 'r') as f:
            content = f.read()
            ntax_match = re.search(r"NTAX\s*=\s*(\d+)", content, re.IGNORECASE)
            nchar_match = re.search(r"NCHAR\s*=\s*(\d+)", content, re.IGNORECASE)
            if ntax_match: unique_count = int(ntax_match.group(1))
            if nchar_match: codons = int(nchar_match.group(1)) // 3
            
    return orig_count, unique_count, codons

def get_stop_stats(repo_path):
    # Find clean fasta used for stops check or original
    # We want INTERNAL stops.
    stops = {"TAA", "TAG", "TGA", "taa", "tag", "tga"}
    fastas = [f for f in (glob.glob(os.path.join(repo_path, "*.fasta")) + glob.glob(os.path.join(repo_path, "*.Y-only"))) if ".unique." not in f and ".cln." not in f and ".clean." not in f]
    
    total_stops = 0
    seqs_with_stops = 0
    
    if fastas:
        with open(fastas[0], 'r') as f:
            seq = []
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if seq:
                        full_seq = "".join(seq)
                        seq_stops = 0
                        for k in range(0, len(full_seq) - 3, 3): # Internal check
                            codon = full_seq[k:k+3]
                            if len(codon) == 3 and codon in stops:
                                seq_stops += 1
                        total_stops += seq_stops
                        if seq_stops > 0: seqs_with_stops += 1
                    seq = []
                else:
                    seq.append(line)
            if seq: # Last seq
                full_seq = "".join(seq)
                seq_stops = 0
                for k in range(0, len(full_seq) - 3, 3):
                    codon = full_seq[k:k+3]
                    if len(codon) == 3 and codon in stops:
                        seq_stops += 1
                total_stops += seq_stops
                if seq_stops > 0: seqs_with_stops += 1
                
    return total_stops, seqs_with_stops

def get_fitmg94_results(repo_path):
    json_path = os.path.join(repo_path, "MG94fit.json")
    omega = "N/A"
    p_val = "N/A"
    g_logl = None
    g_params = None
    
    if os.path.exists(json_path):
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # Global Fit
            fit_global = data.get("fits", {}).get("Global MG94xREV", {}) or data.get("fits", {}).get("Standard MG94", {})
            g_logl = fit_global.get("Log Likelihood")
            g_params = fit_global.get("estimated parameters")
            
            # Omega
            rate_dist = fit_global.get("Rate Distributions", {}) or fit_global.get("fitted parameters", {})
            val = rate_dist.get("non-synonymous/synonymous rate ratio")
            if val is None:
                 # Try finding list format
                 val_list = rate_dist.get("non-synonymous/synonymous rate ratio")
                 if isinstance(val_list, list) and len(val_list)>0:
                     val = val_list[0][0]
            
            if isinstance(val, (int, float)):
                omega = f"{val:.3f}"
            
            # P-value
            test_res = data.get("test results", {}).get("non-synonymous/synonymous rate ratio", {}) or data.get("test results", {}).get("dN/dS != 1", {})
            p = test_res.get("p-value")
            if isinstance(p, (int, float)):
                p_val = f"{p:.4g}"
                
        except:
            pass
            
    return omega, p_val, g_logl, g_params

def get_heterogeneity_results(repo_path, g_logl, g_params):
    json_path = os.path.join(repo_path, "MGlocal.json")
    p_val = "N/A"
    
    if os.path.exists(json_path) and g_logl is not None:
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            fit_local = data.get("fits", {}).get("Local MG94xREV", {})
            l_logl = fit_local.get("Log Likelihood")
            l_params = fit_local.get("estimated parameters")
            
            if l_logl is not None and l_params is not None:
                lrt = 2 * (l_logl - g_logl)
                df = l_params - g_params
                if df > 0:
                    p = chi2.sf(lrt, df)
                    p_val = f"{p:.4g}"
                elif df == 0:
                    p_val = "1.0"
        except:
            pass
            
    return p_val

def get_busted_results(repo_path):
    json_path = os.path.join(repo_path, "BUSTED_best.json")
    p_val = "N/A"
    anomaly = "No"
    mh = "None"
    
    if os.path.exists(json_path):
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # P-value
            test_res = data.get("test results", {})
            p = test_res.get("p-value")
            if isinstance(p, (int, float)):
                p_val = f"{p:.4g}"
            
            # Anomaly
            fit = data.get("fits", {}).get("Unconstrained model", {})
            test_dist = fit.get("Rate Distributions", {}).get("Test", {})
            score = 0.0
            for k, v in test_dist.items():
                if v.get("omega", 0) >= 99.0:
                    score += v.get("proportion", 0)
            if score > 0.001:
                anomaly = "YES"
                
            # Multiple Hits
            settings = data.get("analysis", {}).get("settings", {})
            mh = settings.get("multiple-hit", "None")
            
        except:
            pass
            
    return p_val, anomaly, mh

def get_absrel_results(repo_path):
    json_path = os.path.join(repo_path, "ABSREL.json")
    branches = "N/A"
    if os.path.exists(json_path):
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            branches = data.get("test results", {}).get("positive test results", 0)
        except:
            pass
    return branches

def get_meme_results(repo_path):
    json_path = os.path.join(repo_path, "MEME.json")
    sites = "N/A"
    if os.path.exists(json_path):
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # Extract p-value column index
            mle = data.get("MLE", {})
            headers = mle.get("headers", [])
            content = mle.get("content", {})
            
            pval_idx = 6 # Default
            for i, h in enumerate(headers):
                col_name = h[0] if isinstance(h, list) else h
                if "p-value" in col_name:
                    pval_idx = i
                    break
            
            count = 0
            # Handle partition dict or flat list
            rows = []
            if isinstance(content, dict):
                for k in content: rows.extend(content[k])
            elif isinstance(content, list):
                rows = content
                
            for row in rows:
                if len(row) > pval_idx:
                    p = row[pval_idx]
                    if isinstance(p, (int, float)) and p <= 0.05:
                        count += 1
            sites = count
        except:
            pass
    return sites

# ----------------- Generate HTML -----------------

# ...

html = """
<!DOCTYPE html>
<html>
<head>
<style>
table {
  font-family: Arial, sans-serif;
  border-collapse: collapse;
  width: 100%;
}
td, th {
  border: 1px solid #dddddd;
  text-align: left;
  padding: 8px;
}
tr:nth-child(even) {
  background-color: #f2f2f2;
}
th {
  background-color: #4CAF50;
  color: white;
}
</style>
</head>
<body>

<h2>Y-Chromosome Selection Analysis Results</h2>

<table>
  <tr>
    <th>Gene</th>
    <th>Codons</th>
    <th>Orig/Unique Seqs</th>
    <th>Internal Stops</th>
    <th>Global dN/dS</th>
    <th>FitMG94 p-val</th>
    <th>Heterogeneity p-val</th>
    <th>BUSTED p-val</th>
    <th>BUSTED Anomaly?</th>
    <th>BUSTED MH</th>
    <th>aBSREL Branches</th>
    <th>MEME Sites</th>
  </tr>
"""

for repo in DIRECTORIES:
    repo_path = os.path.join(BASE_DIR, repo)
    if not os.path.exists(repo_path): continue
    
    # Collect Data
    orig, unique, codons = get_alignment_stats(repo_path)
    stops, seqs_stops = get_stop_stats(repo_path)
    omega, mg94_pval, g_logl, g_params = get_fitmg94_results(repo_path)
    het_pval = get_heterogeneity_results(repo_path, g_logl, g_params)
    busted_pval, anomaly, mh = get_busted_results(repo_path)
    absrel_branches = get_absrel_results(repo_path)
    meme_sites = get_meme_results(repo_path)
    
    # Format Row
    html += f"""
  <tr>
    <td>{repo}</td>
    <td>{codons}</td>
    <td>{orig}/{unique}</td>
    <td>{stops}</td>
    <td>{omega}</td>
    <td>{mg94_pval}</td>
    <td>{het_pval}</td>
    <td>{busted_pval}</td>
    <td>{anomaly}</td>
    <td>{mh}</td>
    <td>{absrel_branches}</td>
    <td>{meme_sites}</td>
  </tr>
"""

html += """
</table>

</body>
</html>
"""

with open("summary_table.html", "w") as f:
    f.write(html)

print("Successfully generated summary_table.html")
