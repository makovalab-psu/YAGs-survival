"""
Step 04 — Fit MG94 model (global and local omega) with HyPhy.
  Global: single dN/dS across all branches + LRT test (dN/dS != 1)
  Local:  branch-specific dN/dS; LRT vs global tests rate heterogeneity

Input:  output/{family}_{domain}/{ANALYSIS_ALIGNMENT} + tree.treefile
Output: output/{family}_{domain}/MG94fit.json
        output/{family}_{domain}/MGlocal.json
"""

import os
import json
import subprocess
from pathlib import Path
from scipy.stats import chi2

os.chdir(Path(__file__).parent)
from pipeline_config import HYPHY, FITMG94, ANALYSIS_ALIGNMENT, iter_domains


def run_fitmg94(alignment: Path, tree: Path, output: Path, mode: str = "global") -> bool:
    cmd = [
        HYPHY, FITMG94,
        "--alignment", str(alignment),
        "--tree", str(tree),
        "--type", mode,
        "--rooted", "Yes",
        "--output", str(output),
    ]
    if mode == "global":
        cmd += ["--lrt", "Yes"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    HyPhy FAILED ({mode}):\n{result.stderr[-300:]}")
        return False
    return True


def parse_mg94(json_path: Path, model_key: str) -> tuple[float | None, int | None]:
    """Return (log_likelihood, n_params)."""
    try:
        data = json.loads(json_path.read_text())
        fits = data.get("fits", {})
        fit = fits.get(model_key) or next(
            (v for k, v in fits.items() if "MG94" in k), None
        )
        if fit is None:
            return None, None
        return fit.get("Log Likelihood"), fit.get("estimated parameters")
    except Exception as e:
        print(f"    parse error {json_path.name}: {e}")
        return None, None


def main() -> None:
    results = []
    for d in iter_domains():
        aln = d / ANALYSIS_ALIGNMENT
        tree = d / "tree.treefile"
        if not aln.exists() or not tree.exists():
            print(f"  SKIP {d.name}: missing alignment or tree")
            continue

        global_json = d / "MG94fit.json"
        local_json  = d / "MGlocal.json"

        if not global_json.exists():
            print(f"  {d.name}: running FitMG94 global...")
            run_fitmg94(aln, tree, global_json, "global")
        if not local_json.exists():
            print(f"  {d.name}: running FitMG94 local...")
            run_fitmg94(aln, tree, local_json, "local")

        g_logl, g_params = parse_mg94(global_json, "Global MG94xREV")
        l_logl, l_params = parse_mg94(local_json,  "Local MG94xREV")

        if g_logl and l_logl:
            lrt = 2 * (l_logl - g_logl)
            df  = (l_params or 0) - (g_params or 0)
            p   = chi2.sf(lrt, df) if df > 0 else 1.0
            results.append({"domain": d.name, "g_logl": g_logl, "l_logl": l_logl,
                             "lrt": lrt, "df": df, "p": p})
            print(f"  {d.name}: LRT={lrt:.2f} df={df} p={p:.4g}")

    print("\n| Domain | Global LogL | Local LogL | LRT | df | p (heterogeneity) |")
    print("|--------|-------------|------------|-----|----|-------------------|")
    for r in results:
        print(f"| {r['domain']} | {r['g_logl']:.2f} | {r['l_logl']:.2f} "
              f"| {r['lrt']:.2f} | {r['df']} | {r['p']:.4g} |")


if __name__ == "__main__":
    main()
