#!/usr/bin/env python3
"""
Step 02: Filter MMseqs2 clusters and extract per-cluster FASTA files.

Filtering criteria:
  - Cluster size > min_cluster_size (default 10)
  - At least min_bgc_sequences BGC-derived sequences (AS* or BGC* IDs)

Output: results/02_filtered_clusters/ — one FASTA per passing cluster.
"""

import sys
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path

# ── Config ────────────────────────────────────────────────────────────────
cfg = load_config()
outdir = resolve_path(cfg['output_dir'])
params = cfg['params']

CLUSTER_TSV = outdir / '01_cluster_results' / 'cluster_results_cluster.tsv'
ALL_SEQS_FASTA = outdir / '01_cluster_results' / 'cluster_results_all_seqs.fasta'
OUTPUT_DIR = outdir / '02_filtered_clusters'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_SIZE = params['min_cluster_size']
MIN_BGC = params['min_bgc_sequences']

print("=" * 80)
print("STEP 02: FILTER CLUSTERS AND EXTRACT FASTA")
print("=" * 80)
print(f"  Min cluster size:   {MIN_SIZE}")
print(f"  Min BGC sequences:  {MIN_BGC}")

# ── Parse clusters ────────────────────────────────────────────────────────
clusters = defaultdict(list)
with open(CLUSTER_TSV) as f:
    for line in f:
        rep, member = line.strip().split('\t')
        clusters[rep].append(member)

print(f"\nTotal clusters: {len(clusters):,}")

# ── Filter ────────────────────────────────────────────────────────────────
passing = {}
for rep, members in clusters.items():
    if len(members) <= MIN_SIZE:
        continue
    bgc_count = sum(1 for m in members if m.startswith('AS') or m.startswith('BGC'))
    if bgc_count >= MIN_BGC:
        passing[rep] = set(members)

print(f"Passing clusters: {len(passing):,}")

# ── Index sequences ──────────────────────────────────────────────────────
print("\nIndexing sequences...")
all_needed = set()
for members in passing.values():
    all_needed.update(members)

seq_index = {}
for record in SeqIO.parse(ALL_SEQS_FASTA, 'fasta'):
    if record.id in all_needed:
        seq_index[record.id] = record

print(f"  Indexed {len(seq_index):,} / {len(all_needed):,} needed sequences")

# ── Write per-cluster FASTA ──────────────────────────────────────────────
print("\nWriting per-cluster FASTA files...")
written = 0
for i, (rep, members) in enumerate(sorted(passing.items()), 1):
    records = [seq_index[m] for m in members if m in seq_index]
    if len(records) <= MIN_SIZE:
        continue
    outfile = OUTPUT_DIR / f"clust_{i}.fasta"
    SeqIO.write(records, outfile, 'fasta')
    written += 1
    if written <= 10 or written % 500 == 0:
        print(f"  [{written}] clust_{i}: {len(records)} sequences")

# ── Save cluster membership ──────────────────────────────────────────────
mapping_file = OUTPUT_DIR / 'cluster_membership.tsv'
with open(mapping_file, 'w') as f:
    f.write("cluster_id\trep_id\tmember_id\n")
    for i, (rep, members) in enumerate(sorted(passing.items()), 1):
        for m in sorted(members):
            f.write(f"clust_{i}\t{rep}\t{m}\n")

print(f"\nDone. {written} cluster FASTA files written to {OUTPUT_DIR}")
print(f"Cluster membership: {mapping_file}")
