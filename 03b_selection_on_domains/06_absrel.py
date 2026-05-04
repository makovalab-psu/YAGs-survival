"""
Step 06 — aBSREL: adaptive branch-specific episodic diversifying selection.
Tests each branch for evidence of episodic positive selection.

Input:  output/{family}_{domain}/{ANALYSIS_ALIGNMENT} + tree.treefile
Output: output/{family}_{domain}/ABSREL.json
"""

import os
import json
import subprocess
from pathlib import Path

os.chdir(Path(__file__).parent)
from pipeline_config import HYPHY, ANALYSIS_ALIGNMENT, iter_domains


def run_absrel(alignment: Path, tree: Path, output: Path) -> bool:
    cmd = [
        HYPHY, "absrel",
        "--alignment", str(alignment),
        "--tree", str(tree),
        "--branches", "All",
        "--output", str(output),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    aBSREL FAILED:\n{result.stderr[-300:]}")
        return False
    return True


def parse_absrel(json_path: Path) -> tuple[int, int]:
    """Return (positively_selected_branches, branches_tested)."""
    try:
        data = json.loads(json_path.read_text())
        tr = data.get("test results", {})
        return tr.get("positive test results", 0), tr.get("tested", 0)
    except Exception as e:
        print(f"    parse error {json_path.name}: {e}")
        return 0, 0


def main() -> None:
    results = []
    for d in iter_domains():
        aln  = d / ANALYSIS_ALIGNMENT
        tree = d / "tree.treefile"
        out  = d / "ABSREL.json"
        if not aln.exists() or not tree.exists():
            print(f"  SKIP {d.name}: missing alignment or tree")
            continue
        if out.exists():
            print(f"  SKIP {d.name}: ABSREL.json exists")
        else:
            print(f"  {d.name}: running aBSREL...")
            if not run_absrel(aln, tree, out):
                continue
        pos, tested = parse_absrel(out)
        results.append({"domain": d.name, "pos": pos, "tested": tested})
        print(f"  {d.name}: {pos}/{tested} branches under selection")

    print("\n| Domain | Selected Branches | Tested Branches |")
    print("|--------|-------------------|-----------------|")
    for r in results:
        print(f"| {r['domain']} | {r['pos']} | {r['tested']} |")


if __name__ == "__main__":
    main()
