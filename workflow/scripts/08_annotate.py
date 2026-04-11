#!/usr/bin/env python3
"""
Step 08: Comprehensive annotation of filtered datasets.

Merges metadata from multiple sources with priority-based resolution:
  MITE > MIBiG > SwissProt > AS_db

Also integrates EggNOG (KEGG, COG, PFAMs, EC) and superfamily HMM hits.

This script consolidates the original 7d.rebuild_annotations.py,
14.organise_and_enrich_annotations.py, and 14d.rebuild_auto_annotations.py
into a single pass.

Annotation sources:
  1. MITE — enzyme name, tailoring, cofactors, taxonomy (summary CSV + JSONs)
  2. MIBiG — BGC ID, biosynthetic class, product, taxonomy
  3. SwissProt — UniProt ID, gene name, EC number, PDB IDs, Rhea reactions
  4. AntiSMASH — BGC type, genome ID, taxonomy
  5. EggNOG — KEGG reactions, COG category, PFAMs, EC numbers
  6. Superfamily — HMM-based classification
  7. PanBGC — GCF/OG identifiers (optional)

MITE_ID is only set for PEs whose original_id directly matches a MITE-linked
protein. Dataset_MITE_IDs records cluster-level MITE links (informational).

Output: results/08_annotations/ — one TSV per dataset
"""

import gc
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils import (
    load_config, resolve_path, val, clean_species_name,
    SOURCE_PRIORITY, SOURCE_ORDER, ANNOTATION_COLUMNS, limit_datasets
)

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])
annot_cfg = cfg['annotations']
id_cfg = cfg['id_mapping']

FILTERED_DIR = outdir / '07_filtered_datasets'
OUTPUT_DIR = outdir / '08_annotations'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Data source paths
MITE_SUMMARY = resolve_path(annot_cfg['mite_summary_csv'])
MITE_JSON_DIR = resolve_path(annot_cfg['mite_json_dir'])
MIBIG_FILE = resolve_path(annot_cfg['mibig_annot_tsv'])
AS_FILE = resolve_path(annot_cfg['antismash_annot_tsv'])
SWISS_FILE = resolve_path(annot_cfg['swissprot_annot_tsv'])
EGGNOG_FILE = resolve_path(annot_cfg['eggnog_annotations'])
SF_FILE = resolve_path(annot_cfg['superfamily_hits_tsv'])
PE_MAPPING = resolve_path(id_cfg['pe_mapping_tsv'])

# Optional
PANBGC_FILE = resolve_path(annot_cfg['panbgc_matches_tsv']) if annot_cfg.get('panbgc_matches_tsv') else None
CLUSTER_MAPPING = resolve_path(id_cfg['cluster_mapping_tsv']) if id_cfg.get('cluster_mapping_tsv') else None

print("=" * 80)
print("STEP 08: COMPREHENSIVE ANNOTATION")
print("=" * 80)


# ===== LOAD DATA SOURCES =====

# 1. PE ID mapping
print("1. Loading PE ID mapping...")
new_to_original = {}
original_to_new = {}
with open(PE_MAPPING) as f:
    next(f)  # header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            orig, new_id = parts[0], parts[1]
            new_to_original[new_id] = orig
            original_to_new[orig] = new_id
print(f"   {len(new_to_original):,} mappings")

# 2. MITE summary
print("2. Loading MITE summary...")
df_mite = pd.read_csv(MITE_SUMMARY)
mite_annotations = {}
for rec in df_mite.to_dict('records'):
    mid = val(rec.get('accession', ''))
    if not mid:
        continue
    mite_annotations[mid] = {
        'Source': 'MITE', 'MITE_ID': mid,
        'gene_name': val(rec.get('name', '')),
        'enzyme_name': val(rec.get('name', '')),
        'tailoring': val(rec.get('tailoring', '')),
        'description': val(rec.get('description', '')),
        'reaction_description': val(rec.get('reaction_description', '')),
        'cofactors_organic': val(rec.get('cofactors_organic', '')),
        'cofactors_inorganic': val(rec.get('cofactors_inorganic', '')),
        'organism': clean_species_name(rec.get('organism', '')),
        'domain': val(rec.get('domain', '')),
        'phylum': val(rec.get('phylum', '')),
        'class': val(rec.get('class', '')),
        'order': val(rec.get('order', '')),
        'family': val(rec.get('family', '')),
    }
print(f"   {len(mite_annotations):,} MITE records")
del df_mite; gc.collect()

# 3. MITE JSON ID mappings
print("3. Loading MITE JSON ID mappings...")
mite_id_lookup = defaultdict(set)
mite_protein_ids = defaultdict(set)
json_files = sorted(MITE_JSON_DIR.glob('MITE*.json'))

for jf in json_files:
    try:
        with open(jf) as f:
            data = json.load(f)
        mid = data.get('accession', '')
        if not mid or 'enzyme' not in data:
            continue
        db = data['enzyme'].get('databaseIds', {})
        uniprot = db.get('uniprot', '')
        genpept = db.get('genpept', '')
        mibig = db.get('mibig', '')
        # Index only protein-level IDs (uniprot, genpept) — NOT bare MIBiG IDs
        # which are cluster-level and would cause false MITE matches for
        # all proteins in the same BGC.
        for pid in [uniprot, genpept]:
            if pid:
                mite_id_lookup[pid].add(mid)
                mite_protein_ids[mid].add(pid)
        if mibig:
            mite_protein_ids[mid].add(mibig)  # track for info, not for lookup
        # Composite keys (BGC_protein, MITE_protein) for exact original_id matching
        for combo in [f"{mibig}_{uniprot}", f"{mibig}_{genpept}",
                      f"{mid}_{genpept}", f"{mid}_{uniprot}"]:
            if combo and '_' in combo and not combo.startswith('_'):
                mite_id_lookup[combo].add(mid)
                mite_protein_ids[mid].add(combo)
    except Exception:
        pass

mite_id_lookup = {k: list(v) for k, v in mite_id_lookup.items()}
print(f"   {len(mite_id_lookup):,} lookup keys, {len(mite_protein_ids):,} MITE accessions")

# 4. MIBiG
print("4. Loading MIBiG annotations...")
df_mibig = pd.read_csv(MIBIG_FILE, sep='\t', low_memory=False)
mibig_annotations = defaultdict(list)
for rec in df_mibig.to_dict('records'):
    seq_id = val(rec.get('ID', ''))
    if seq_id:
        mibig_annotations[seq_id].append({
            'Source': 'MIBiG', 'MIBiG_ID': val(rec.get('BGC', '')),
            'BGC_type': val(rec.get('biosyn_class', '')),
            'gene_name': val(rec.get('name', '')),
            'enzyme_name': val(rec.get('name', '')),
            'tailoring': val(rec.get('function_category', '')),
            'product': val(rec.get('product', '')),
            'organism': clean_species_name(rec.get('organizm_name', '')),
            'Species': clean_species_name(rec.get('Species', '')),
        })
print(f"   {len(mibig_annotations):,} MIBiG entries")
del df_mibig; gc.collect()

# 5. SwissProt
print("5. Loading SwissProt annotations...")
df_swiss = pd.read_csv(SWISS_FILE, sep='\t', low_memory=False)
swiss_annotations = defaultdict(list)
for rec in df_swiss.to_dict('records'):
    seq_id = val(rec.get('ID', ''))
    if seq_id:
        swiss_annotations[seq_id].append({
            'Source': 'SwissProt', 'Uniprot_ID': seq_id,
            'Entry_Name': val(rec.get('Entry Name', '')),
            'gene_name': val(rec.get('Gene', '')),
            'ProteinExistence': val(rec.get('ProteinExistence', '')),
            'EC_Number': val(rec.get('EC_Number', '')),
            'PDB_IDs': val(rec.get('PDB_IDs', '')),
            'Species': clean_species_name(rec.get('Species', '')),
            'Rhea': val(rec.get('Rhea', '')),
            'Superkingdom': val(rec.get('Superkingdom', '')),
            'BGC_type': val(rec.get('biosyn_class', '')),
            'Enzyme_function': val(rec.get('Enzyme_function', '')),
        })
print(f"   {len(swiss_annotations):,} SwissProt entries")
del df_swiss; gc.collect()

# 6. AntiSMASH
print("6. Loading AntiSMASH annotations...")
df_as = pd.read_csv(AS_FILE, sep='\t', low_memory=False)
as_annotations = defaultdict(list)
for rec in df_as.to_dict('records'):
    seq_id = val(rec.get('AS_id', ''))
    if seq_id:
        as_annotations[seq_id].append({
            'Source': 'AS_db', 'BGC_type': val(rec.get('BGC_type', '')),
            'genome_ID': val(rec.get('genome_ID', '')),
            'Species': clean_species_name(rec.get('strain', '')),
            'order': val(rec.get('Order name', '')),
            'class': val(rec.get('Class name', '')),
            'phylum': val(rec.get('Phylum name', '')),
            'Superkingdom': val(rec.get('Superkingdom name', '')),
        })
print(f"   {len(as_annotations):,} AS entries")
del df_as; gc.collect()

# 7. EggNOG
print("7. Loading EggNOG annotations...")
eggnog_annotations = {}
with open(EGGNOG_FILE) as fh:
    header = None
    for line in fh:
        if line.startswith('##'):
            continue
        line = line.rstrip('\n')
        if header is None:
            header = line.lstrip('#').split('\t')
            continue
        parts = line.split('\t')
        if len(parts) < len(header):
            continue
        row = dict(zip(header, parts))
        pe_id = re.sub(r'^P?T\d+_', '', row.get('query', ''))
        if pe_id:
            eggnog_annotations[pe_id] = {
                'KEGG_Reaction': val(row.get('KEGG_Reaction', '')).replace('-', ''),
                'KEGG_rclass': val(row.get('KEGG_rclass', '')).replace('-', ''),
                'COG_category': val(row.get('COG_category', '')).replace('-', ''),
                'PFAMs': val(row.get('PFAMs', '')).replace('-', ''),
                'eggnog_EC': val(row.get('EC', '')).replace('-', ''),
            }
print(f"   {len(eggnog_annotations):,} EggNOG records")

# 8. Superfamily
print("8. Loading superfamily annotations...")
superfamily_annotations = {}
df_sf = pd.read_csv(SF_FILE, sep='\t')
for rec in df_sf.to_dict('records'):
    pe_id = val(rec.get('PE_ID', ''))
    if pe_id:
        superfamily_annotations[pe_id] = {
            'Superfamily': val(rec.get('best_superfamily', '')),
            'All_superfamilies': val(rec.get('all_superfamilies', '')),
        }
print(f"   {len(superfamily_annotations):,} superfamily records")
del df_sf; gc.collect()

# 9. PanBGC (optional)
panbgc_data = {}
if PANBGC_FILE and PANBGC_FILE.exists():
    print("9. Loading PanBGC matches...")
    df_pb = pd.read_csv(PANBGC_FILE, sep='\t')
    for rec in df_pb.to_dict('records'):
        pe_id = val(rec.get('PE_ID', ''))
        mappings = val(rec.get('panBGC_mappings', ''))
        status = val(rec.get('match_status', ''))
        if pe_id and mappings and status == 'mapped':
            # Take only the first mapping if multiple (separated by ;)
            first_mapping = mappings.split(';')[0].strip()
            panbgc_data[pe_id] = first_mapping
    print(f"   {len(panbgc_data):,} PanBGC matches")
    del df_pb; gc.collect()
else:
    print("9. PanBGC: not provided (skipping)")


# ===== ANNOTATION FUNCTIONS =====

def resolve_mite_ids(original_id):
    """Find MITE accessions linked to an original_id.

    Only returns MITE IDs for direct protein-level matches:
      - Full original_id in mite_id_lookup (covers protein accessions
        and composite keys like BGC_protein or MITE_protein)
      - MITE accession prefix in underscore-separated IDs

    Does NOT split original_id by '_' to check individual parts against
    mite_id_lookup, because that would incorrectly match all proteins
    sharing a BGC with a MITE enzyme (the broadcast bug).
    """
    found = set()
    # 1. Direct full-ID lookup (covers uniprot, genpept, and composite keys)
    if original_id in mite_id_lookup:
        found.update(mite_id_lookup[original_id])
    # 2. If the ID itself is or starts with a MITE accession
    if original_id.startswith('MITE'):
        parts = original_id.split('_', 1)
        if parts[0] in mite_annotations:
            found.add(parts[0])
    # 3. For composite IDs, only extract embedded MITE accessions —
    #    do NOT do mite_id_lookup on individual parts.
    elif '_' in original_id:
        for part in original_id.split('_'):
            if part.startswith('MITE') and part in mite_annotations:
                found.add(part)
    return list(found)


def find_annotations_for_id(original_id):
    """Find all annotations for an original ID from all sources."""
    annot = {'mite': [], 'mibig': [], 'swiss': [], 'as_db': []}
    mite_ids = resolve_mite_ids(original_id)
    for mid in mite_ids:
        if mid in mite_annotations:
            annot['mite'].append(mite_annotations[mid])
    if original_id in mibig_annotations:
        annot['mibig'].extend(mibig_annotations[original_id])
    if original_id in swiss_annotations:
        annot['swiss'].extend(swiss_annotations[original_id])
    if '_' in original_id:
        for part in original_id.split('_'):
            if part in swiss_annotations and part != original_id:
                annot['swiss'].extend(swiss_annotations[part])
    if original_id in as_annotations:
        annot['as_db'].extend(as_annotations[original_id])
    return annot


def merge_annotations(annot_dict):
    """Merge with priority MITE > MIBiG > SwissProt > AS_db."""
    merged = {col: '' for col in ANNOTATION_COLUMNS}
    if annot_dict['mite']:
        merged['Source'] = 'MITE'
    elif annot_dict['mibig']:
        merged['Source'] = 'MIBiG'
    elif annot_dict['swiss']:
        merged['Source'] = 'SwissProt'
    elif annot_dict['as_db']:
        merged['Source'] = 'AS_db'
    else:
        return merged
    # Fill: highest priority first (values stick)
    for source_key in SOURCE_ORDER:
        for entry in annot_dict[source_key]:
            for key, value in entry.items():
                if key == 'Source':
                    continue
                v = val(value)
                if v and key in merged and not merged[key]:
                    merged[key] = v
    return merged


def annotate_pe(pe_id, all_pe_ids=None):
    """Build annotation row for a PE. Checks all dedup-group members."""
    if all_pe_ids is None:
        all_pe_ids = [pe_id]
    combined = {'mite': [], 'mibig': [], 'swiss': [], 'as_db': []}
    for pid in all_pe_ids:
        orig = new_to_original.get(pid, pid)
        annot = find_annotations_for_id(orig)
        for src in combined:
            combined[src].extend(annot[src])
    merged = merge_annotations(combined)
    merged['ID'] = pe_id
    # EggNOG
    for pid in all_pe_ids:
        egg = eggnog_annotations.get(pid, {})
        for col in ['KEGG_Reaction', 'KEGG_rclass', 'COG_category', 'PFAMs', 'eggnog_EC']:
            if not merged[col] and col in egg:
                merged[col] = egg[col]
    # Superfamily
    for pid in all_pe_ids:
        sf = superfamily_annotations.get(pid, {})
        for col in ['Superfamily', 'All_superfamilies']:
            if not merged[col] and col in sf:
                merged[col] = sf[col]
    # Cluster: MIBiG_ID first, then PanBGC mapping as fallback
    if not merged.get('Cluster'):
        if merged.get('MIBiG_ID'):
            merged['Cluster'] = merged['MIBiG_ID']
        else:
            for pid in all_pe_ids:
                pb = panbgc_data.get(pid, '')
                if pb:
                    merged['Cluster'] = pb
                    break
    return merged


# ===== PROCESS DATASETS =====

print("\n" + "=" * 80)
print("ANNOTATING DATASETS")
print("=" * 80)

dataset_files = sorted(FILTERED_DIR.glob('*.fasta'))
# Exclude non-dataset files
dataset_files = [f for f in dataset_files if not f.name.startswith('kept_') and not f.name.startswith('excl')]
dataset_files = limit_datasets(dataset_files, cfg, label='datasets')
print(f"Found {len(dataset_files):,} datasets")

processed = 0
total_seqs = 0

for ds_file in dataset_files:
    ds_name = ds_file.stem
    output_file = OUTPUT_DIR / f"{ds_name}_annotations.tsv"

    # Read sequences
    sequences = list(SeqIO.parse(ds_file, 'fasta'))
    if not sequences:
        continue

    # Annotate each sequence
    rows = []
    dataset_mite_ids = set()
    for seq in sequences:
        row = annotate_pe(seq.id, [seq.id])
        rows.append(row)
        # Collect MITE IDs at dataset level (for informational field)
        if row.get('MITE_ID'):
            dataset_mite_ids.add(row['MITE_ID'])

    # Set Dataset_MITE_IDs for all rows (informational, cluster-level)
    ds_mite_str = ';'.join(sorted(dataset_mite_ids)) if dataset_mite_ids else ''
    for row in rows:
        row['Dataset_MITE_IDs'] = ds_mite_str

    df = pd.DataFrame(rows, columns=ANNOTATION_COLUMNS)
    df.to_csv(output_file, sep='\t', index=False)

    processed += 1
    total_seqs += len(sequences)
    if processed <= 10 or processed % 500 == 0:
        n_ann = int((df['Source'] != '').sum())
        print(f"  [{processed}/{len(dataset_files)}] {ds_name}: "
              f"{len(sequences)} seqs, {n_ann} annotated")

print(f"\nDone. {processed} datasets annotated ({total_seqs:,} sequences).")
print(f"Output: {OUTPUT_DIR}")
