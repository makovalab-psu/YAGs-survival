"""
Step 07 — MEME: mixed effects model of evolution for site-specific episodic selection.
Identifies individual codon sites under episodic positive selection.

Input:  output/{family}_{domain}/{ANALYSIS_ALIGNMENT} + tree.treefile
Output: output/{family}_{domain}/MEME.json
"""

import os
import json
import subprocess
from pathlib import Path

os.chdir(Path(__file__).parent)
from pipeline_config import HYPHY, ANALYSIS_ALIGNMENT, iter_domains


def run_meme(alignment: Path, tree: Path, output: Path) -> bool:
    cmd = [
        HYPHY, "meme",
        "--alignment", str(alignment),
        "--tree", str(tree),
        "--output", str(output),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    MEME FAILED:\n{result.stderr[-300:]}")
        return False
    return True


def parse_meme(json_path: Path, p_threshold: float = 0.05) -> tuple[int, int]:
    """Return (selected_sites, total_sites) at given p-value threshold."""
    try:
        data = json.loads(json_path.read_text())
        mle = data.get("MLE", {})
        headers = mle.get("headers", [])
        content = mle.get("content", {})

        # Find p-value column index
        pval_idx = next(
            (i for i, h in enumerate(headers)
             if (h[0] if isinstance(h, list) else h) == "p-value"),
            6  # fallback
        )

        all_rows = []
        if isinstance(content, dict):
            for part in content.values():
                all_rows.extend(part)
        else:
            all_rows = content

        selected = sum(
            1 for row in all_rows
            if len(row) > pval_idx and isinstance(row[pval_idx], (int, float))
            and row[pval_idx] <= p_threshold
        )
        return selected, len(all_rows)
    except Exception as e:
        print(f"    parse error {json_path.name}: {e}")
        return 0, 0


def main() -> None:
    results = []
    for d in iter_domains():
        aln  = d / ANALYSIS_ALIGNMENT
        tree = d / "tree.treefile"
        out  = d / "MEME.json"
        if not aln.exists() or not tree.exists():
            print(f"  SKIP {d.name}: missing alignment or tree")
            continue
        if out.exists():
            print(f"  SKIP {d.name}: MEME.json exists")
        else:
            print(f"  {d.name}: running MEME...")
            if not run_meme(aln, tree, out):
                continue
        sel, total = parse_meme(out)
        results.append({"domain": d.name, "sel": sel, "total": total})
        print(f"  {d.name}: {sel}/{total} sites under selection (p≤0.05)")

    print("\n| Domain | Selected Sites (p≤0.05) | Total Sites |")
    print("|--------|------------------------|-------------|")
    for r in results:
        print(f"| {r['domain']} | {r['sel']} | {r['total']} |")


if __name__ == "__main__":
    main()
