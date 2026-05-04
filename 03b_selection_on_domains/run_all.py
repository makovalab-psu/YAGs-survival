"""
Master script — runs the full selection pipeline in order.
Each step is skipped if its output already exists (idempotent re-runs).

Steps:
  00  align_domains.py   — MAFFT domain alignments (run manually with gene family arg)
  00b import_alignments  — copy pre-built alignments from test_alignments/ into output/
  01  clean_stops        — replace in-frame stop codons with gaps
  02  deduplicate        — [SKIPPED] identical domain sequences are real gene copies;
                           HyPhy 2.5+ handles zero-length branches without issues
  03  build_trees        — IQtree2 with outgroup rooting
  04  fitmg94            — MG94 global + local dN/dS
  05  busted             — BUSTED episodic selection (AIC-optimized)
  06  absrel             — aBSREL branch-specific selection
  07  meme               — MEME site-specific selection
  08  summarize          — collect all results into summary table

Usage:
    python run_all.py                           # run all steps, all domains
    python run_all.py 01 03                     # run only steps 01 and 03, all domains
    python run_all.py 05 06 07 --only DAZ       # run steps 05-07 for DAZ-DAZL_RBD only
    python run_all.py --only CDY 03 04          # steps and --only can be in any order
    python run_all.py 01 03 --dir output_full   # run against a different output directory
"""

import sys
import os
import importlib
from pathlib import Path

os.chdir(Path(__file__).parent)

STEPS = {
    "00a": ("00_setup_from_alignments", "00_setup_from_alignments.py"),
    "00b": ("00_import_alignments", "00_import_alignments.py"),
    "01": ("01_clean_stops",   "01_clean_stops.py"),
    "02": ("02_deduplicate",   "02_deduplicate.py"),
    "03": ("03_build_trees",   "03_build_trees.py"),
    "04": ("04_fitmg94",       "04_fitmg94.py"),
    "05": ("05_busted",        "05_busted.py"),
    "06": ("06_absrel",        "06_absrel.py"),
    "07": ("07_meme",          "07_meme.py"),
    "08": ("08_summarize",     "08_summarize.py"),
    "09": ("09_plot_trees",    "09_plot_trees.py"),
}

# Parse --only and --dir from args
args = sys.argv[1:]
domain_filter = None
output_dir = None

for flag in ("--only", "--dir"):
    if flag in args:
        idx = args.index(flag)
        val = args[idx + 1]
        args = args[:idx] + args[idx + 2:]
        if flag == "--only":
            domain_filter = val
        else:
            output_dir = val

requested = args if args else list(STEPS.keys())

import pipeline_config

# Override output directory if specified
if output_dir:
    pipeline_config.OUTPUT_DIR = Path(output_dir)
    print(f"Output dir: '{output_dir}'")

# Inject domain filter into pipeline_config so iter_domains() respects it
if domain_filter:
    _orig_iter = pipeline_config.iter_domains
    pipeline_config.iter_domains = lambda: [
        d for d in _orig_iter() if domain_filter.lower() in d.name.lower()
    ]
    print(f"Domain filter: '{domain_filter}'")

for step_id in requested:
    if step_id not in STEPS:
        print(f"Unknown step '{step_id}'. Valid: {list(STEPS.keys())}")
        sys.exit(1)
    module_name, script = STEPS[step_id]
    print(f"\n{'='*60}")
    print(f"STEP {step_id}: {script}")
    print("="*60)
    mod = importlib.import_module(module_name)
    mod.main()
