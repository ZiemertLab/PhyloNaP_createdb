#!/usr/bin/env python3
"""
Step 06: Create sub-clusters, align, trim, and QC.

Splits multi-cluster datasets into per-cluster FASTA files based on TreeCluster
results. Single-cluster datasets are copied as-is. Sub-clusters below
min_subcluster_size are discarded.

Then runs AliStat on each to compute Ca/Cr quality metrics.

Output: results/06_subclusters/ — aligned and trimmed sub-cluster FASTAs
        results/06_subclusters/alistat_stats.tsv — quality metrics
"""

import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path, get_tool, limit_datasets

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])
params = cfg['params']

TREECLUSTER_DIR = outdir / '05_treecluster'
TRIMMED_DIR = outdir / '03_trimmed_alignments'
FASTA_DIR = outdir / '02_filtered_clusters'
OUTPUT_DIR = outdir / '06_subclusters'
ALIGNED_DIR = OUTPUT_DIR / 'aligned'
TRIMMED_OUT = OUTPUT_DIR / 'trimmed'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
ALIGNED_DIR.mkdir(exist_ok=True)
TRIMMED_OUT.mkdir(exist_ok=True)

MIN_SIZE = params['min_subcluster_size']
MAFFT = get_tool(cfg, 'mafft')
TRIMAL = get_tool(cfg, 'trimal')
ALISTAT = get_tool(cfg, 'alistat')
MAX_ITER = params['mafft_maxiterate']

print("=" * 80)
print("STEP 06: CREATE SUBCLUSTERS, ALIGN, TRIM, QC")
print("=" * 80)
print(f"  Min subcluster size: {MIN_SIZE}")

# ── Parse TreeCluster results ─────────────────────────────────────────────
tc_files = sorted(TREECLUSTER_DIR.glob('*_clusters.tsv'))
tc_files = limit_datasets(tc_files, cfg, label='TreeCluster files')
print(f"\nFound {len(tc_files)} TreeCluster result files")

all_subclusters = []  # list of (base_name, subcluster_id, seq_ids)

for tc_file in tc_files:
    base_name = tc_file.stem.replace('_clusters', '')
    clusters = defaultdict(list)
    with open(tc_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                seq_id, cluster_id = parts[0], parts[1]
                if cluster_id == '-1':
                    cluster_id = '0'  # singleton → main cluster
                clusters[cluster_id].append(seq_id)

    if len(clusters) <= 1:
        # Single cluster — use original sequences
        all_ids = []
        for ids in clusters.values():
            all_ids.extend(ids)
        if len(all_ids) >= MIN_SIZE:
            all_subclusters.append((base_name, None, all_ids))
    else:
        # Multiple clusters — split
        for cid, seq_ids in sorted(clusters.items()):
            if len(seq_ids) >= MIN_SIZE:
                all_subclusters.append((base_name, cid, seq_ids))

print(f"Total subclusters (>= {MIN_SIZE} seqs): {len(all_subclusters)}")

# ── Extract subcluster FASTAs from original FASTAs ────────────────────────
print("\nExtracting subcluster FASTAs...")
written = 0

for base_name, cid, seq_ids in all_subclusters:
    if cid is not None:
        sc_name = f"{base_name}_c{cid}"
    else:
        sc_name = base_name

    outfile = OUTPUT_DIR / f"{sc_name}.fasta"
    if outfile.exists():
        written += 1
        continue

    # Find source FASTA
    src = FASTA_DIR / f"{base_name}.fasta"
    if not src.exists():
        continue

    id_set = set(seq_ids)
    records = [r for r in SeqIO.parse(src, 'fasta') if r.id in id_set]
    if len(records) >= MIN_SIZE:
        SeqIO.write(records, outfile, 'fasta')
        written += 1

print(f"  Written: {written} subcluster FASTAs")

# ── Align and trim each subcluster ────────────────────────────────────────
print("\nAligning and trimming subclusters...")
processed = 0

for fasta in sorted(OUTPUT_DIR.glob('*.fasta')):
    sc_name = fasta.stem
    aligned = ALIGNED_DIR / f"{sc_name}_aligned.fasta"
    trimmed = TRIMMED_OUT / f"{sc_name}_trimmed.fasta"

    if trimmed.exists():
        processed += 1
        continue

    # Clean non-standard AAs
    cleaned_lines = []
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                cleaned_lines.append(line)
            else:
                import re
                cleaned_lines.append(re.sub(r'[JjBbZzOoUu]', 'X', line))
    cleaned_path = OUTPUT_DIR / f"{sc_name}_clean.fasta"
    with open(cleaned_path, 'w') as f:
        f.writelines(cleaned_lines)

    # Align
    try:
        with open(aligned, 'w') as al_fh:
            subprocess.run(
                [MAFFT, '--auto', '--maxiterate', str(MAX_ITER), '--quiet', str(cleaned_path)],
                stdout=al_fh, stderr=subprocess.PIPE, text=True, check=True, timeout=1800,
            )
    except Exception as e:
        print(f"  WARN: MAFFT failed for {sc_name}: {e}")
        cleaned_path.unlink(missing_ok=True)
        continue

    # Trim
    try:
        subprocess.run(
            [TRIMAL, '-in', str(aligned), '-automated1', '-out', str(trimmed)],
            capture_output=True, check=True, timeout=600
        )
    except Exception:
        try:
            subprocess.run(
                [TRIMAL, '-in', str(aligned), '-gt', '0.5', '-out', str(trimmed)],
                capture_output=True, check=True, timeout=600
            )
        except Exception:
            import shutil
            shutil.copy2(aligned, trimmed)

    cleaned_path.unlink(missing_ok=True)
    processed += 1

    if processed <= 10 or processed % 500 == 0:
        print(f"  [{processed}/{len(all_subclusters)}] {sc_name}")

# ── Run AliStat ──────────────────────────────────────────────────────────
print("\nRunning AliStat quality assessment...")
stats_file = OUTPUT_DIR / 'alistat_stats.tsv'
stats_rows = []

for trimmed in sorted(TRIMMED_OUT.glob('*_trimmed.fasta')):
    sc_name = trimmed.stem.replace('_trimmed', '')
    try:
        result = subprocess.run(
            [ALISTAT, '-i', str(trimmed), '-t', '6'],
            capture_output=True, text=True, timeout=120
        )
        # Parse AliStat output for Ca and Cr values
        ca = cr_min = n_seqs = al_len = ''
        for line in result.stdout.split('\n'):
            if 'Ca:' in line:
                ca = line.split(':')[-1].strip()
            elif 'Cr(min):' in line or 'Cr_min:' in line:
                cr_min = line.split(':')[-1].strip()
            elif 'Number of sequences' in line:
                n_seqs = line.split(':')[-1].strip()
            elif 'Alignment length' in line:
                al_len = line.split(':')[-1].strip()
        stats_rows.append(f"{sc_name}\t{n_seqs}\t{al_len}\t{ca}\t{cr_min}")
    except Exception:
        pass

with open(stats_file, 'w') as f:
    f.write("dataset\tn_sequences\talignment_length\tCa\tCr_min\n")
    for row in stats_rows:
        f.write(row + '\n')

print(f"\nDone. {len(stats_rows)} datasets assessed.")
print(f"  Subclusters: {OUTPUT_DIR}")
print(f"  Trimmed:     {TRIMMED_OUT}")
print(f"  Stats:       {stats_file}")
