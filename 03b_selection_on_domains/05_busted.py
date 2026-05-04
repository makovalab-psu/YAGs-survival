"""
Step 05 — BUSTED: test for episodic diversifying selection gene-wide.
Uses AIC-based hill-climbing to select the best number of rate classes (N, M)
and tests for improvement with multiple-hit models.

Input:  output/{family}_{domain}/{ANALYSIS_ALIGNMENT} + tree.treefile
Output: output/{family}_{domain}/BUSTED_best.json
"""

import os
import json
import shutil
import subprocess
from pathlib import Path

os.chdir(Path(__file__).parent)
from pipeline_config import HYPHY, ANALYSIS_ALIGNMENT, iter_domains


def run_busted(
    alignment: Path, tree: Path, output: Path,
    rates: int, syn_rates: int, multiple_hits: str = "None"
) -> bool:
    cmd = [
        HYPHY, "busted",
        "--alignment", str(alignment),
        "--tree", str(tree),
        "--rates", str(rates),
        "--error-sink", "Yes",
        "--starting-points", "1",
        "--multiple-hits", multiple_hits,
        "--output", str(output),
    ]
    cmd += ["--srv", "No"] if syn_rates == 1 else ["--srv", "Yes", "--syn-rates", str(syn_rates)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def parse_busted(json_path: Path) -> tuple[float | None, float, float | None]:
    """Return (aic_c, anomaly_score, p_value)."""
    try:
        data = json.loads(json_path.read_text())
        fit = data.get("fits", {}).get("Unconstrained model", {})
        aic = fit.get("AIC-c")
        omega_dist = fit.get("Rate Distributions", {}).get("non-synonymous/synonymous rate ratio", [])
        anomaly = sum(e[1] for e in omega_dist if isinstance(e, list) and len(e) >= 2 and e[0] >= 99.0)
        p_val = data.get("test results", {}).get("p-value")
        return aic, anomaly, p_val
    except Exception as e:
        print(f"    parse error {json_path.name}: {e}")
        return None, 0.0, None


def optimize_busted(d: Path, aln: Path, tree: Path) -> None:
    final_json = d / "BUSTED_best.json"
    if final_json.exists():
        print(f"  SKIP {d.name}: BUSTED_best.json exists")
        return

    tmp = d / "BUSTED_temp.json"
    best_aic, best_n, best_m, best_mh = float("inf"), 1, 1, "None"

    print(f"  {d.name}: initial run N=1 M=1...")
    if not run_busted(aln, tree, tmp, 1, 1):
        print(f"  {d.name}: initial run FAILED")
        return

    aic, anomaly, _ = parse_busted(tmp)
    if aic is None:
        print(f"  {d.name}: could not parse initial AIC")
        return

    best_aic = aic
    shutil.copy(tmp, final_json)
    print(f"    AIC-c={aic:.2f} anomaly={anomaly:.4f}")

    # Hill-climb on rate classes
    n, m = 1, 1
    while True:
        improved = False
        for (cn, cm) in [(n + 1, m), (n, m + 1)]:
            if cn > 6 or cm > 6:
                continue
            print(f"    testing N={cn} M={cm}...")
            if run_busted(aln, tree, tmp, cn, cm):
                aic, anomaly, _ = parse_busted(tmp)
                if aic is not None:
                    print(f"      AIC-c={aic:.2f}")
                    if aic < best_aic - 0.01:
                        best_aic, best_n, best_m = aic, cn, cm
                        shutil.copy(tmp, final_json)
                        improved = True
        if improved:
            n, m = best_n, best_m
        else:
            break

    # Test multiple hits
    print(f"    testing Double+Triple hits on N={best_n} M={best_m}...")
    if run_busted(aln, tree, tmp, best_n, best_m, "Double+Triple"):
        aic, anomaly, _ = parse_busted(tmp)
        if aic is not None and aic < best_aic - 0.01:
            best_aic, best_mh = aic, "Double+Triple"
            shutil.copy(tmp, final_json)
            print(f"      MH improved AIC-c → {aic:.2f}")

    if tmp.exists():
        tmp.unlink()

    _, anomaly, p_val = parse_busted(final_json)
    p_str = f"{p_val:.4g}" if isinstance(p_val, float) else str(p_val)
    print(f"  {d.name}: best N={best_n} M={best_m} MH={best_mh} "
          f"AIC-c={best_aic:.2f} anomaly={anomaly:.4f} p={p_str}")


def main() -> None:
    for d in iter_domains():
        aln  = d / ANALYSIS_ALIGNMENT
        tree = d / "tree.treefile"
        if not aln.exists() or not tree.exists():
            print(f"  SKIP {d.name}: missing alignment or tree")
            continue
        optimize_busted(d, aln, tree)


if __name__ == "__main__":
    main()
