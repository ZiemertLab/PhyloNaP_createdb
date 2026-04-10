#!/usr/bin/env python3
"""
Step 09: Deduplicate identical sequences within each dataset.

For each filtered alignment:
  1. Group sequences by gap-stripped amino acid string
  2. Keep representative with highest-priority source (MITE>MIBiG>SwissProt>AS_db)
  3. Aggregate annotations across group members (pipe-separated)
  4. Track n_deduplicated, deduplicated_ids, all_sources

Output: results/09_deduplicated_alignments/
        results/09_deduplicated_annotations/
"""

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path, SOURCE_PRIORITY

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])

ALIGNMENT_DIR = outdir / '07_filtered_datasets'
ANNOTATION_DIR = outdir / '08_annotations'
OUT_ALIGN_DIR = outdir / '09_deduplicated_alignments'
OUT_ANNOT_DIR = outdir / '09_deduplicated_annotations'

OUT_ALIGN_DIR.mkdir(parents=True, exist_ok=True)
OUT_ANNOT_DIR.mkdir(parents=True, exist_ok=True)

TEST_MODE = '--test' in sys.argv

print("=" * 80)
print("STEP 09: DEDUPLICATE IDENTICAL SEQUENCES")
print("=" * 80)


def aggregate_values(values):
    """Pipe-join unique non-empty values."""
    seen = set()
    unique = []
    for v in values:
        s = str(v).strip() if v and pd.notna(v) else ''
        if s and s.lower() != 'nan' and s not in seen:
            seen.add(s)
            unique.append(s)
    return '|'.join(unique) if unique else ''


alignment_files = sorted(ALIGNMENT_DIR.glob('*.fasta'))
alignment_files = [f for f in alignment_files if not f.name.startswith(('kept_', 'excl'))]
print(f"Found {len(alignment_files)} alignment files")

if TEST_MODE:
    alignment_files = alignment_files[:5]

total_in = total_out = total_removed = 0
processed = 0

for al_file in alignment_files:
    ds_name = al_file.stem
    annot_file = ANNOTATION_DIR / f"{ds_name}_annotations.tsv"
    if not annot_file.exists():
        continue

    sequences = list(SeqIO.parse(al_file, 'fasta'))
    df_annot = pd.read_csv(annot_file, sep='\t', dtype=str).fillna('')

    # Group by gap-stripped sequence
    seq_groups = defaultdict(list)
    for rec in sequences:
        key = str(rec.seq).replace('-', '').upper()
        seq_groups[key].append(rec.id)

    dedup_seqs = []
    dedup_annots = []

    for seq_str, seq_ids in seq_groups.items():
        # Get annotations
        annots = df_annot[df_annot.iloc[:, 0].isin(seq_ids)].copy()
        if annots.empty:
            # Keep first sequence anyway
            best_id = seq_ids[0]
        else:
            annots['_priority'] = annots.get('Source', pd.Series(dtype=str)).map(
                lambda x: SOURCE_PRIORITY.get(str(x).strip(), 0)
            )
            annots = annots.sort_values('_priority', ascending=False)
            best_id = annots.iloc[0].iloc[0]  # ID column

        # Keep best sequence
        best_seq = next((s for s in sequences if s.id == best_id), None)
        if best_seq:
            dedup_seqs.append(best_seq)

        # Aggregate annotations
        if len(annots) > 0:
            agg = {'n_deduplicated': len(seq_ids), 'deduplicated_ids': '|'.join(seq_ids)}
            agg['all_sources'] = aggregate_values(annots.get('Source', []))
            # Keep best row's values, aggregate others
            best_row = annots.iloc[0].drop('_priority', errors='ignore').to_dict()
            for col in annots.columns:
                if col == '_priority':
                    continue
                if col in ('n_deduplicated', 'deduplicated_ids', 'all_sources'):
                    continue
                if len(seq_ids) > 1:
                    agg[col] = aggregate_values(annots[col])
                else:
                    agg[col] = best_row.get(col, '')
            dedup_annots.append(agg)
        else:
            dedup_annots.append({
                df_annot.columns[0]: best_id,
                'n_deduplicated': len(seq_ids),
                'deduplicated_ids': '|'.join(seq_ids),
                'all_sources': '',
            })

    # Save
    out_al = OUT_ALIGN_DIR / f"{ds_name}.fasta"
    SeqIO.write(dedup_seqs, out_al, 'fasta')

    df_dedup = pd.DataFrame(dedup_annots)
    out_an = OUT_ANNOT_DIR / f"{ds_name}_annotations.tsv"
    df_dedup.to_csv(out_an, sep='\t', index=False)

    n_removed = len(sequences) - len(dedup_seqs)
    total_in += len(sequences)
    total_out += len(dedup_seqs)
    total_removed += n_removed
    processed += 1

    if processed <= 10 or processed % 500 == 0:
        print(f"  [{processed}/{len(alignment_files)}] {ds_name}: "
              f"{len(sequences)} -> {len(dedup_seqs)} (-{n_removed})")

print(f"\nDone.")
print(f"  Datasets: {processed}")
print(f"  Input sequences:  {total_in:,}")
print(f"  Output sequences: {total_out:,}")
print(f"  Removed: {total_removed:,} ({total_removed/total_in*100:.1f}%)" if total_in else "")
print(f"  Output: {OUT_ALIGN_DIR}")
