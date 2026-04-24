"""
Step 08 — Collect results from all HyPhy analyses into a summary table.
Reads existing JSON outputs; does not re-run any analyses.

Output: printed markdown table + output/summary.tsv
"""

import os
import json
from pathlib import Path
from scipy.stats import chi2

os.chdir(Path(__file__).parent)
from pipeline_config import OUTPUT_DIR, iter_domains


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def mg94_stats(d: Path) -> dict:
    g = load_json(d / "MG94fit.json")
    l = load_json(d / "MGlocal.json")

    g_fits = g.get("fits", {})
    l_fits = l.get("fits", {})

    def get_fit(fits, key):
        return fits.get(key) or next((v for k, v in fits.items() if "MG94" in k), {})

    gf = get_fit(g_fits, "Global MG94xREV")
    lf = get_fit(l_fits, "Local MG94xREV")

    # dN/dS from global + 95% CI
    omega = (gf.get("Rate Distributions", {})
               .get("non-synonymous/synonymous rate ratio", "N/A"))
    omega_str = f"{omega:.3f}" if isinstance(omega, (int, float)) else "N/A"

    ci = (gf.get("Confidence Intervals", {})
            .get("non-synonymous/synonymous rate ratio", {}))
    lb, ub = ci.get("LB"), ci.get("UB")
    omega_ci = f"[{lb:.3f}, {ub:.3f}]" if lb is not None and ub is not None else "N/A"

    # LRT dN/dS != 1 (from --lrt Yes in global run)
    omega_p = g.get("test results", {}).get(
        "non-synonymous/synonymous rate ratio", {}
    ).get("p-value")
    omega_p_str = f"{omega_p:.4g}" if isinstance(omega_p, (int, float)) else "N/A"

    # LRT global vs local
    g_ll = gf.get("Log Likelihood")
    l_ll = lf.get("Log Likelihood")
    g_np = gf.get("estimated parameters")
    l_np = lf.get("estimated parameters")
    lrt_p = "N/A"
    if g_ll and l_ll and g_np and l_np:
        lrt = 2 * (l_ll - g_ll)
        df  = l_np - g_np
        lrt_p = f"{chi2.sf(lrt, df):.4g}" if df > 0 else "1.0"

    return {"omega": omega_str, "omega_ci": omega_ci, "omega_p": omega_p_str, "lrt_p": lrt_p}


def busted_stats(d: Path) -> dict:
    data = load_json(d / "BUSTED_best.json")
    p = data.get("test results", {}).get("p-value", "N/A")
    fit = data.get("fits", {}).get("Unconstrained model", {})
    omega_dist = fit.get("Rate Distributions", {}).get("non-synonymous/synonymous rate ratio", [])
    anomaly = sum(e[1] for e in omega_dist if isinstance(e, list) and len(e) >= 2 and e[0] >= 99.0)
    return {
        "busted_p": f"{p:.4g}" if isinstance(p, float) else str(p),
        "anomaly": f"{anomaly:.4f}" if anomaly else "0",
    }


def absrel_stats(d: Path) -> dict:
    data = load_json(d / "ABSREL.json")
    tr = data.get("test results", {})
    pos = tr.get("positive test results", "N/A")
    tested = tr.get("tested", "N/A")
    return {"absrel_pos": str(pos), "absrel_tested": str(tested)}


def meme_stats(d: Path, p_thr: float = 0.05) -> dict:
    data = load_json(d / "MEME.json")
    mle = data.get("MLE", {})
    headers = mle.get("headers", [])
    content = mle.get("content", {})
    pval_idx = next(
        (i for i, h in enumerate(headers)
         if (h[0] if isinstance(h, list) else h) == "p-value"), 6
    )
    all_rows = []
    if isinstance(content, dict):
        for part in content.values():
            all_rows.extend(part)
    else:
        all_rows = content
    sel = sum(
        1 for row in all_rows
        if len(row) > pval_idx and isinstance(row[pval_idx], (int, float))
        and row[pval_idx] <= p_thr
    )
    return {"meme_sel": str(sel) if all_rows else "N/A", "meme_total": str(len(all_rows))}


def main() -> None:
    rows = []
    for d in iter_domains():
        row = {"domain": d.name}
        row.update(mg94_stats(d))
        row.update(busted_stats(d))
        row.update(absrel_stats(d))
        row.update(meme_stats(d))
        rows.append(row)

    headers = ["domain", "omega", "omega_ci", "omega_p", "lrt_p", "busted_p", "anomaly",
               "absrel_pos", "absrel_tested", "meme_sel", "meme_total"]
    labels  = ["Domain", "dN/dS", "dN/dS 95% CI", "p (dN/dS≠1)", "LRT p (heterogeneity)",
               "BUSTED p", "BUSTED anomaly", "aBSREL positives", "aBSREL tested",
               "MEME sites (p≤0.05)", "MEME total sites"]

    # Markdown
    print("| " + " | ".join(labels) + " |")
    print("|" + "|".join(["---"] * len(labels)) + "|")
    for r in rows:
        print("| " + " | ".join(r.get(h, "N/A") for h in headers) + " |")

    # TSV
    tsv_path = OUTPUT_DIR / "summary.tsv"
    with open(tsv_path, "w") as f:
        f.write("\t".join(labels) + "\n")
        for r in rows:
            f.write("\t".join(r.get(h, "N/A") for h in headers) + "\n")
    print(f"\nSaved: {tsv_path}")


if __name__ == "__main__":
    main()
