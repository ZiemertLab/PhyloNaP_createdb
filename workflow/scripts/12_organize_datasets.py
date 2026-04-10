#!/usr/bin/env python3
"""
Step 12: Organize final datasets into PT-numbered folders.

Combines original steps 11 (unaligned extraction) and 17 (PT ID assignment +
folder renaming + sequence merging).

For each final filtered dataset:
  1. Extracts unaligned sequences from the original FASTA files
  2. Assigns a PT ID (PT######) based on sequential numbering
  3. Creates a PT folder with all dataset artifacts

Output structure per dataset:
  results/12_datasets/PT######/
    PT######.fasta              Unaligned sequences
    PT######_al.fa              Full alignment (trimmed)
    PT######.nwk                Tree (to be rooted in step 13)
    PT######_annotations.tsv    Annotation table

Also produces:
  results/12_datasets/pt_id_mapping.tsv
  results/12_datasets/all_sequences_merged.fasta
"""

import shutil
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])

ALIGNMENT_DIR = outdir / '10_final_alignments'
ANNOTATION_DIR = outdir / '10_final_annotations'
TREE_DIR = outdir / '11_trees'
CLUSTER_FASTA_DIR = outdir / '02_filtered_clusters'
OUTPUT_DIR = outdir / '12_datasets'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("STEP 12: ORGANIZE DATASETS (unaligned extraction + PT IDs)")
print("=" * 80)

# ── Step 1: Build index of all original (unaligned) sequences ──────────────
print("\n--- Step 1: Indexing original sequences ---")
seq_index = {}
n_indexed = 0
for fasta_file in sorted(CLUSTER_FASTA_DIR.glob('*.fasta')):
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id not in seq_index:
            seq_index[record.id] = record
            n_indexed += 1

# Also check the initial merged FASTA if cluster files don't cover everything
initial_fasta = outdir / '01_cluster_results' / 'combined_initial_data.fasta'
if initial_fasta.exists():
    for record in SeqIO.parse(initial_fasta, 'fasta'):
        if record.id not in seq_index:
            seq_index[record.id] = record
            n_indexed += 1

print(f"  Indexed {n_indexed:,} sequences")

# ── Step 2: Gather all datasets ───────────────────────────────────────────
alignment_files = sorted(ALIGNMENT_DIR.glob('*.fasta'))
print(f"\n--- Step 2: Processing {len(alignment_files)} datasets ---")

mapping_rows = []
all_sequences_out = OUTPUT_DIR / 'all_sequences_merged.fasta'
merged_fh = open(all_sequences_out, 'w')

stats = {'ok': 0, 'no_tree': 0, 'missing_seqs': 0}

for idx, al_file in enumerate(alignment_files, 1):
    ds_name = al_file.stem
    pt_id = f"PT{idx:06d}"

    # Create PT folder
    pt_dir = OUTPUT_DIR / pt_id
    pt_dir.mkdir(exist_ok=True)

    # 3a. Copy trimmed alignment
    shutil.copy2(al_file, pt_dir / f'{pt_id}_al.fa')

    # 3b. Copy annotation
    annot_src = ANNOTATION_DIR / f'{ds_name}_annotations.tsv'
    if annot_src.exists():
        shutil.copy2(annot_src, pt_dir / f'{pt_id}_annotations.tsv')

    # 3c. Copy tree (unrooted)
    tree_src = TREE_DIR / f'{ds_name}.tree'
    if tree_src.exists():
        shutil.copy2(tree_src, pt_dir / f'{pt_id}_unrooted.nwk')
    else:
        stats['no_tree'] += 1

    # 3d. Extract unaligned sequences
    seq_ids = [r.id for r in SeqIO.parse(al_file, 'fasta')]
    unaligned = []
    missing = 0
    for sid in seq_ids:
        if sid in seq_index:
            unaligned.append(seq_index[sid])
        else:
            missing += 1

    if missing > 0:
        stats['missing_seqs'] += missing

    if unaligned:
        SeqIO.write(unaligned, pt_dir / f'{pt_id}.fasta', 'fasta')
        # Write to merged FASTA with dataset prefix
        for rec in unaligned:
            merged_fh.write(f'>{pt_id}|{rec.id}\n{str(rec.seq)}\n')

    mapping_rows.append((pt_id, ds_name))

    stats['ok'] += 1
    if idx <= 20 or idx % 500 == 0:
        print(f"  [{idx}/{len(alignment_files)}] {ds_name} -> {pt_id} "
              f"(n={len(seq_ids)})")

merged_fh.close()

# ── Save mapping ──────────────────────────────────────────────────────────
mapping_file = OUTPUT_DIR / 'pt_id_mapping.tsv'
with open(mapping_file, 'w') as f:
    f.write('pt_id\toriginal_name\n')
    for pt_id, ds_name in mapping_rows:
        f.write(f'{pt_id}\t{ds_name}\n')

print(f"\nDone.")
print(f"  Datasets organized: {stats['ok']}")
print(f"  Missing trees:      {stats['no_tree']}")
print(f"  Missing sequences:  {stats['missing_seqs']}")
print(f"  Mapping:            {mapping_file}")
print(f"  Merged FASTA:       {all_sequences_out}")
print(f"  Output:             {OUTPUT_DIR}")
