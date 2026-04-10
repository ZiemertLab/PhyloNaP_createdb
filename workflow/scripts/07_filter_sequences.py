#!/usr/bin/env python3
"""
Step 07: Cascade filtration — filter bad sequences, re-align, re-trim, QC gate.

For each subcluster:
  1. Run AliStat to get per-sequence Cr values
  2. Remove sequences with Cr < cr_threshold
  3. Re-align (MAFFT) and re-trim (trimAl) the filtered set
  4. Apply quality gate:
     - N sequences >= min_seqs_after_filter
     - Ca > ca_threshold
     - Retention ratio >= retention_threshold
     - Mean original sequence length > min_mean_seq_length

Output: results/07_filtered_datasets/  — datasets passing all QC
        results/07_filtered_datasets/filter_report.tsv
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path, get_tool

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])
params = cfg['params']

TRIMMED_DIR = outdir / '06_subclusters' / 'trimmed'
RAW_FASTA_DIR = outdir / '06_subclusters'
OUTPUT_DIR = outdir / '07_filtered_datasets'
REFILTERED_DIR = OUTPUT_DIR / 'realigned'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
REFILTERED_DIR.mkdir(exist_ok=True)

CR_THRESH = params['cr_threshold']
CA_THRESH = params['ca_threshold']
RET_THRESH = params['retention_threshold']
MIN_LEN = params['min_mean_seq_length']
MIN_SEQS = params['min_seqs_after_filter']

MAFFT = get_tool(cfg, 'mafft')
TRIMAL = get_tool(cfg, 'trimal')
ALISTAT = get_tool(cfg, 'alistat')
MAX_ITER = params['mafft_maxiterate']

print("=" * 80)
print("STEP 07: CASCADE FILTRATION — FILTER BAD SEQUENCES")
print("=" * 80)
print(f"  Cr threshold:     {CR_THRESH}")
print(f"  Ca threshold:     {CA_THRESH}")
print(f"  Retention thresh: {RET_THRESH}")
print(f"  Min mean length:  {MIN_LEN}")
print(f"  Min sequences:    {MIN_SEQS}")


def run_alistat_per_seq(fasta_path):
    """Run AliStat and parse per-sequence Cr values, Ca, and alignment length."""
    try:
        result = subprocess.run(
            [ALISTAT, '-i', str(fasta_path), '-t', '6', '-s'],
            capture_output=True, text=True, timeout=120
        )
        ca = 0.0
        al_len = 0
        seq_cr = {}
        in_per_seq = False

        for line in result.stdout.split('\n'):
            line = line.strip()
            if 'Ca:' in line:
                try:
                    ca = float(line.split(':')[-1].strip())
                except ValueError:
                    pass
            elif 'Alignment length' in line:
                try:
                    al_len = int(line.split(':')[-1].strip())
                except ValueError:
                    pass
            elif line.startswith('Seq_ID') or line.startswith('seq_id'):
                in_per_seq = True
                continue
            elif in_per_seq and '\t' in line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    try:
                        seq_cr[parts[0]] = float(parts[1])
                    except ValueError:
                        pass
        return ca, al_len, seq_cr
    except Exception:
        return 0.0, 0, {}


def compute_mean_raw_length(seq_ids, raw_fasta_dir, base_name):
    """Compute mean original (unaligned) sequence length."""
    # Try to find original raw FASTA
    raw = raw_fasta_dir / f"{base_name}.fasta"
    if not raw.exists():
        return 999  # don't penalize if we can't find raw

    id_set = set(seq_ids)
    lengths = []
    for r in SeqIO.parse(raw, 'fasta'):
        if r.id in id_set:
            lengths.append(len(str(r.seq).replace('-', '')))
    return sum(lengths) / len(lengths) if lengths else 0


# ── Process each trimmed alignment ────────────────────────────────────────
trimmed_files = sorted(TRIMMED_DIR.glob('*_trimmed.fasta'))
print(f"\nProcessing {len(trimmed_files)} trimmed alignments...")

report_rows = []
kept = 0
excluded = 0

for i, trimmed in enumerate(trimmed_files, 1):
    sc_name = trimmed.stem.replace('_trimmed', '')
    base_name = sc_name.split('_c')[0] if '_c' in sc_name else sc_name

    # 1. Get per-sequence Cr
    ca, al_len, seq_cr = run_alistat_per_seq(trimmed)

    # 2. Filter bad sequences
    good_ids = {sid for sid, cr in seq_cr.items() if cr >= CR_THRESH}
    all_ids = set(seq_cr.keys())
    if not all_ids:
        # Fallback: keep all
        all_ids = {r.id for r in SeqIO.parse(trimmed, 'fasta')}
        good_ids = all_ids

    removed = len(all_ids) - len(good_ids)

    # 3. Check if enough sequences remain
    if len(good_ids) < MIN_SEQS:
        report_rows.append({
            'dataset': sc_name, 'status': 'excluded',
            'reason': f'too_few_after_filter ({len(good_ids)})',
            'n_original': len(all_ids), 'n_filtered': len(good_ids),
            'n_removed': removed, 'Ca': ca
        })
        excluded += 1
        continue

    # 4. If sequences were removed, re-align and re-trim
    if removed > 0:
        # Extract filtered sequences (unaligned) from raw FASTA
        raw = RAW_FASTA_DIR / f"{sc_name}.fasta"
        if not raw.exists():
            raw = RAW_FASTA_DIR / f"{base_name}.fasta"

        if raw.exists():
            records = [r for r in SeqIO.parse(raw, 'fasta') if r.id in good_ids]
        else:
            # Fall back to trimmed alignment, strip gaps
            records = []
            for r in SeqIO.parse(trimmed, 'fasta'):
                if r.id in good_ids:
                    from Bio.SeqRecord import SeqRecord
                    from Bio.Seq import Seq
                    records.append(SeqRecord(
                        Seq(str(r.seq).replace('-', '')), id=r.id, description=''
                    ))

        if len(records) < MIN_SEQS:
            report_rows.append({
                'dataset': sc_name, 'status': 'excluded',
                'reason': f'too_few_records ({len(records)})',
                'n_original': len(all_ids), 'n_filtered': len(records),
                'n_removed': removed, 'Ca': ca
            })
            excluded += 1
            continue

        # Write filtered FASTA
        filt_fasta = REFILTERED_DIR / f"{sc_name}_filtered.fasta"
        SeqIO.write(records, filt_fasta, 'fasta')

        # Re-align
        realigned = REFILTERED_DIR / f"{sc_name}_aligned.fasta"
        try:
            with open(realigned, 'w') as out:
                subprocess.run(
                    [MAFFT, '--auto', '--maxiterate', str(MAX_ITER), '--quiet', str(filt_fasta)],
                    stdout=out, stderr=subprocess.DEVNULL, check=True, timeout=1800
                )
        except Exception:
            report_rows.append({
                'dataset': sc_name, 'status': 'excluded', 'reason': 'realign_failed',
                'n_original': len(all_ids), 'n_filtered': len(records),
                'n_removed': removed, 'Ca': ca
            })
            excluded += 1
            continue

        # Re-trim
        retrimmed = REFILTERED_DIR / f"{sc_name}_trimmed.fasta"
        try:
            subprocess.run(
                [TRIMAL, '-in', str(realigned), '-automated1', '-out', str(retrimmed)],
                capture_output=True, check=True, timeout=600
            )
        except Exception:
            try:
                subprocess.run(
                    [TRIMAL, '-in', str(realigned), '-gt', '0.5', '-out', str(retrimmed)],
                    capture_output=True, check=True, timeout=600
                )
            except Exception:
                import shutil
                shutil.copy2(realigned, retrimmed)

        final_trimmed = retrimmed
    else:
        final_trimmed = trimmed

    # 5. Quality gate on final alignment
    new_ca, new_al_len, _ = run_alistat_per_seq(final_trimmed)

    # Count sequences in final
    n_final = sum(1 for _ in SeqIO.parse(final_trimmed, 'fasta'))

    # Retention ratio
    n_original_seqs = sum(1 for _ in SeqIO.parse(trimmed, 'fasta')) if removed > 0 else len(all_ids)
    retention = n_final / n_original_seqs if n_original_seqs > 0 else 1.0

    # Mean raw length
    mean_len = compute_mean_raw_length(good_ids, RAW_FASTA_DIR, base_name)

    passes = True
    fail_reasons = []
    if n_final < MIN_SEQS:
        passes = False; fail_reasons.append(f'n_seqs={n_final}')
    if new_ca < CA_THRESH:
        passes = False; fail_reasons.append(f'Ca={new_ca:.3f}')
    if retention < RET_THRESH:
        passes = False; fail_reasons.append(f'retention={retention:.3f}')
    if mean_len < MIN_LEN:
        passes = False; fail_reasons.append(f'mean_len={mean_len:.0f}')

    if passes:
        # Copy to output
        import shutil
        shutil.copy2(final_trimmed, OUTPUT_DIR / f"{sc_name}.fasta")
        kept += 1
        report_rows.append({
            'dataset': sc_name, 'status': 'kept', 'reason': '',
            'n_original': len(all_ids), 'n_filtered': n_final,
            'n_removed': removed, 'Ca': new_ca
        })
    else:
        excluded += 1
        report_rows.append({
            'dataset': sc_name, 'status': 'excluded',
            'reason': '; '.join(fail_reasons),
            'n_original': len(all_ids), 'n_filtered': n_final,
            'n_removed': removed, 'Ca': new_ca
        })

    if i <= 10 or i % 500 == 0:
        status = 'KEPT' if passes else 'EXCL'
        print(f"  [{i}/{len(trimmed_files)}] {sc_name}: {status} (n={n_final}, Ca={new_ca:.3f})")

# ── Save report ──────────────────────────────────────────────────────────
report_file = OUTPUT_DIR / 'filter_report.tsv'
pd.DataFrame(report_rows).to_csv(report_file, sep='\t', index=False)

# Also save list of kept/excluded datasets
kept_list = OUTPUT_DIR / 'kept_datasets.txt'
excl_list = OUTPUT_DIR / 'excluded_datasets.txt'
with open(kept_list, 'w') as f:
    for r in report_rows:
        if r['status'] == 'kept':
            f.write(r['dataset'] + '\n')
with open(excl_list, 'w') as f:
    f.write("dataset\treason\n")
    for r in report_rows:
        if r['status'] == 'excluded':
            f.write(f"{r['dataset']}\t{r['reason']}\n")

print(f"\nDone.")
print(f"  Kept:     {kept}")
print(f"  Excluded: {excluded}")
print(f"  Report:   {report_file}")
print(f"  Output:   {OUTPUT_DIR}")
