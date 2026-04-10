#!/usr/bin/env python3
"""
Step 10: Filter datasets for database inclusion (biosynthetic relevance).

Filtering criteria:
  - Must have >= 1 sequence from AS_db or MIBiG (ensures BGC linkage)

NOTE: Unlike the original pipeline, eukaryotic datasets are NOT excluded.
All datasets passing the BGC-linkage criterion are kept regardless of
taxonomic composition.

Output: results/10_final_alignments/
        results/10_final_annotations/
"""

import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])
params = cfg['params']

IN_ALIGN = outdir / '09_deduplicated_alignments'
IN_ANNOT = outdir / '09_deduplicated_annotations'
OUT_ALIGN = outdir / '10_final_alignments'
OUT_ANNOT = outdir / '10_final_annotations'

OUT_ALIGN.mkdir(parents=True, exist_ok=True)
OUT_ANNOT.mkdir(parents=True, exist_ok=True)

REQUIRE_BGC = params.get('require_bgc_linkage', True)

print("=" * 80)
print("STEP 10: FILTER DATASETS FOR DATABASE INCLUSION")
print("=" * 80)
print(f"  Require BGC linkage: {REQUIRE_BGC}")
print(f"  Eukaryotic filtering: DISABLED (all kingdoms included)")


def has_bgc_linkage(annot_file):
    """Check if dataset has at least one AS_db or MIBiG entry."""
    try:
        df = pd.read_csv(annot_file, sep='\t', dtype=str)
    except Exception:
        return False
    if df.empty:
        return False
    # Check all_sources column (from deduplication) or Source
    if 'all_sources' in df.columns:
        return df['all_sources'].apply(
            lambda x: 'AS_db' in str(x) or 'MIBiG' in str(x)
        ).any()
    if 'Source' in df.columns:
        return df['Source'].isin(['AS_db', 'MIBiG']).any()
    return False


annot_files = sorted(IN_ANNOT.glob('*_annotations.tsv'))
print(f"\nProcessing {len(annot_files)} datasets...")

kept = excluded = missing = 0
kept_list = []
excluded_list = []

for annot_file in annot_files:
    ds_name = annot_file.stem.replace('_annotations', '')
    al_file = IN_ALIGN / f"{ds_name}.fasta"

    if not al_file.exists():
        missing += 1
        continue

    keep = True
    reason = 'passed'

    if REQUIRE_BGC and not has_bgc_linkage(annot_file):
        keep = False
        reason = 'no_AS_or_MIBiG'

    if keep:
        shutil.copy2(al_file, OUT_ALIGN / al_file.name)
        shutil.copy2(annot_file, OUT_ANNOT / annot_file.name)
        kept += 1
        kept_list.append(ds_name)
    else:
        excluded += 1
        excluded_list.append((ds_name, reason))

    if (kept + excluded) <= 10 or (kept + excluded) % 1000 == 0:
        status = 'KEPT' if keep else f'EXCL ({reason})'
        print(f"  [{kept + excluded}/{len(annot_files)}] {ds_name}: {status}")

# Save lists
with open(OUT_ANNOT / 'kept_datasets.txt', 'w') as f:
    for ds in kept_list:
        f.write(ds + '\n')

with open(OUT_ANNOT / 'excluded_datasets.txt', 'w') as f:
    f.write("dataset\treason\n")
    for ds, reason in excluded_list:
        f.write(f"{ds}\t{reason}\n")

print(f"\nDone.")
print(f"  Kept:     {kept} ({kept/(kept+excluded)*100:.1f}%)" if kept + excluded > 0 else "")
print(f"  Excluded: {excluded}")
print(f"  Missing alignment: {missing}")
print(f"  Output: {OUT_ALIGN}")
