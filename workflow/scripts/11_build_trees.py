#!/usr/bin/env python3
"""
Step 11: Build phylogenetic trees using FastTree (LG + gamma).

For each final filtered alignment, builds an ML phylogeny.
Requires >= 4 sequences per dataset.

Output: results/11_trees/ — one Newick tree per dataset
"""

import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path, get_tool, limit_datasets

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])
params = cfg['params']

INPUT_DIR = outdir / '10_final_alignments'
OUTPUT_DIR = outdir / '11_trees'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FASTTREE = get_tool(cfg, 'fasttree')
MODEL_OPTS = params['fasttree_model'].split()
TIMEOUT = params['fasttree_timeout']
MIN_SEQS = params['min_seqs_for_tree']

print("=" * 80)
print("STEP 11: BUILD PHYLOGENETIC TREES")
print("=" * 80)
print(f"  Model: {' '.join(MODEL_OPTS)}")
print(f"  Timeout: {TIMEOUT}s")
print(f"  Min sequences: {MIN_SEQS}")

alignment_files = sorted(INPUT_DIR.glob('*.fasta'))
alignment_files = limit_datasets(alignment_files, cfg, label='alignments')
print(f"\nFound {len(alignment_files)} alignments")

stats = {'success': 0, 'skipped': 0, 'failed': 0, 'few_seqs': 0}
total_time = 0

for i, al_file in enumerate(alignment_files, 1):
    ds_name = al_file.stem
    tree_file = OUTPUT_DIR / f"{ds_name}.tree"

    if tree_file.exists():
        stats['skipped'] += 1
        continue

    # Count sequences
    n_seqs = sum(1 for line in open(al_file) if line.startswith('>'))
    if n_seqs < MIN_SEQS:
        stats['few_seqs'] += 1
        continue

    # Run FastTree
    cmd = [FASTTREE] + MODEL_OPTS + [str(al_file)]
    t0 = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        dt = time.time() - t0
        total_time += dt

        if result.returncode == 0 and result.stdout.strip():
            with open(tree_file, 'w') as f:
                f.write(result.stdout)
            stats['success'] += 1
            if stats['success'] <= 20 or stats['success'] % 500 == 0:
                print(f"  [{i}/{len(alignment_files)}] {ds_name}: "
                      f"{n_seqs} seqs, {dt:.1f}s")
        else:
            stats['failed'] += 1
    except subprocess.TimeoutExpired:
        stats['failed'] += 1
    except Exception:
        stats['failed'] += 1

print(f"\nDone.")
print(f"  Success:   {stats['success']}")
print(f"  Skipped:   {stats['skipped']} (existing)")
print(f"  Too few:   {stats['few_seqs']}")
print(f"  Failed:    {stats['failed']}")
print(f"  Total time: {total_time/60:.1f} min")
print(f"  Output: {OUTPUT_DIR}")
