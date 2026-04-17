#!/usr/bin/env python3
"""Run mutant_notation.py and sift_table.py for all gene folders (skipping venv).

Usage: python sift_table_all.py
Output:
  <gene>/mutant_notation.txt  -- consensus FASTA + mutation lines
  <gene>/sift_table.tsv       -- SIFT results in gapped coordinates (if sift files present)
"""

import subprocess
import sys
from pathlib import Path

SKIP = {'venv'}
BASE = Path(__file__).parent
MUTANT_SCRIPT = BASE / 'mutant_notation.py'
SIFT_SCRIPT = BASE / 'sift_table.py'


def run(script: Path, gene: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(script), gene],
        capture_output=True, text=True, cwd=BASE
    )


def main() -> None:
    gene_dirs = sorted(d for d in BASE.iterdir() if d.is_dir() and d.name not in SKIP)

    combined_rows: list[str] = []
    header = "Gene\tMSA_pos\tOrig_AA\tChanged_AA\tSIFT_prediction\tScore\tConfidence"

    for gene_dir in gene_dirs:
        gene = gene_dir.name

        # mutant_notation.py requires *.clean.fasta and mutated_positions.txt
        has_fasta = bool(list(gene_dir.glob('*.clean.fasta')))
        has_positions = (gene_dir / 'mutated_positions.txt').exists()
        if has_fasta and has_positions:
            result = run(MUTANT_SCRIPT, gene)
            if result.returncode != 0:
                sys.stderr.write(f"mutant_notation error [{gene}]:\n{result.stderr}\n")
            else:
                (gene_dir / 'mutant_notation.txt').write_text(result.stdout)
                sys.stderr.write(result.stderr)
                print(f"Written: {gene_dir / 'mutant_notation.txt'}")
        else:
            sys.stderr.write(f"Skipping mutant_notation [{gene}]: missing fasta or mutated_positions.txt\n")

        # sift_table.py requires sift_input and sift_result
        if (gene_dir / 'sift_input').exists() and (gene_dir / 'sift_result').exists():
            result = run(SIFT_SCRIPT, gene)
            if result.returncode != 0:
                sys.stderr.write(f"sift_table error [{gene}]:\n{result.stderr}\n")
            else:
                (gene_dir / 'sift_table.tsv').write_text(result.stdout)
                sys.stderr.write(result.stderr)
                print(f"Written: {gene_dir / 'sift_table.tsv'}")
                # collect data rows (skip header line) for combined table
                for line in result.stdout.splitlines()[1:]:
                    if line.strip():
                        combined_rows.append(f"{gene}\t{line}")
        else:
            sys.stderr.write(f"Skipping sift_table [{gene}]: missing sift_input or sift_result\n")

    if combined_rows:
        combined_file = BASE / 'sift_table_all.tsv'
        combined_file.write_text('\n'.join([header] + combined_rows) + '\n')
        print(f"Written: {combined_file}")


if __name__ == "__main__":
    main()
