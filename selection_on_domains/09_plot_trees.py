"""Step 09 — Plot phylogenetic trees (calls 09_plot_trees.R)."""

import os
import subprocess
from pathlib import Path

os.chdir(Path(__file__).parent)


def main() -> None:
    import pipeline_config
    result = subprocess.run(
        ["Rscript", "09_plot_trees.R", str(pipeline_config.OUTPUT_DIR)],
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError("09_plot_trees.R failed")


if __name__ == "__main__":
    main()
