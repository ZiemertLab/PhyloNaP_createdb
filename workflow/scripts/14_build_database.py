#!/usr/bin/env python3
"""
Step 14: Build the PhyloNaP database (JSON + SQLite).

Reads organized PT datasets from step 12, groups them by superfamily
annotation, and produces:
  - results/14_database/db_structure.json
  - results/14_database/phylonap.db  (SQLite)

Optionally merges curated datasets from curated.json.
Also generates genome-location mapping tables.
"""

import json
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path, val, limit_datasets

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])

DATASETS_DIR = outdir / '12_datasets'
ANNOT_DIR = outdir / '10_final_annotations'
ALISTAT_FILE = outdir / '07_filtered_datasets' / 'filter_report.tsv'
MAPPING_FILE = DATASETS_DIR / 'pt_id_mapping.tsv'
OUTPUT_DIR = outdir / '14_database'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_JSON = OUTPUT_DIR / 'db_structure.json'
OUTPUT_DB = OUTPUT_DIR / 'phylonap.db'

# Optional curated data
CURATED_JSON = None
curated_json_path = cfg.get('curated_json', '')
if curated_json_path:
    p = resolve_path(curated_json_path)
    if p.exists():
        CURATED_JSON = p

print("=" * 80)
print("STEP 14: BUILD PHYLONAP DATABASE")
print("=" * 80)
print(f"  Datasets:   {DATASETS_DIR}")
print(f"  Curated:    {CURATED_JSON or 'none'}")
print(f"  Output:     {OUTPUT_DIR}")

# ── Metadata columns ─────────────────────────────────────────────────────
METADATA_COLUMNS = [
    'Source', 'MITE_ID', 'Uniprot_ID', 'MIBiG_ID',
    'gene_name', 'enzyme_name', 'tailoring',
    'description', 'reaction_description',
    'cofactors_organic', 'cofactors_inorganic',
    'product', 'Rhea', 'Enzyme_function',
    'Entry_Name', 'ProteinExistence', 'EC_Number', 'PDB_IDs',
    'organism', 'Species', 'domain', 'Superkingdom',
    'phylum', 'class', 'order', 'family',
    'BGC_type', 'genome_ID',
    'Superfamily', 'All_superfamilies',
    'KEGG_Reaction', 'KEGG_rclass', 'COG_category', 'PFAMs', 'eggnog_EC',
]


def parse_annotation(annot_path):
    """Parse annotation file and return dataset statistics."""
    sources = []
    names = []
    col_has_data = {c: False for c in METADATA_COLUMNS}
    sf_labels = set()
    ds_pfams = set()
    ds_cog = set()
    ds_ec = set()

    with open(annot_path) as f:
        header = next(f).strip().split('\t')
        col_idx = {c: header.index(c) for c in METADATA_COLUMNS if c in header}
        sf_idx = col_idx.get('Superfamily')
        pfams_idx = col_idx.get('PFAMs')
        cog_idx = col_idx.get('COG_category')
        ec_idx = col_idx.get('EC_Number')
        eggnog_ec_idx = col_idx.get('eggnog_EC')

        for line in f:
            parts = line.rstrip('\n').split('\t')
            while len(parts) < len(header):
                parts.append('')
            row = dict(zip(header, parts))

            src = row.get('Source', '').strip()
            sources.append(src)

            for field in ('enzyme_name', 'gene_name', 'Enzyme_function'):
                v = row.get(field, '').strip()
                if v and v not in ('', 'nan'):
                    names.append(v)
                    break

            for col, idx in col_idx.items():
                if not col_has_data[col]:
                    v = parts[idx].strip() if idx < len(parts) else ''
                    if v and v != 'nan':
                        col_has_data[col] = True

            if sf_idx is not None:
                sv = parts[sf_idx].strip() if sf_idx < len(parts) else ''
                if sv and sv != 'nan':
                    sf_labels.add(sv)

            if pfams_idx is not None:
                pv = parts[pfams_idx].strip() if pfams_idx < len(parts) else ''
                if pv and pv != 'nan':
                    for v in pv.split(','):
                        v = v.strip()
                        if v:
                            ds_pfams.add(v)
            if cog_idx is not None:
                cv = parts[cog_idx].strip() if cog_idx < len(parts) else ''
                if cv and cv != 'nan':
                    for ch in cv:
                        if ch.strip():
                            ds_cog.add(ch)
            for eidx in (ec_idx, eggnog_ec_idx):
                if eidx is not None:
                    ev = parts[eidx].strip() if eidx < len(parts) else ''
                    if ev and ev != 'nan':
                        for v in ev.split(','):
                            v = v.strip()
                            if v:
                                ds_ec.add(v)

    n_total = len(sources)
    n_char = sum(1 for s in sources if s in ('MITE', 'SwissProt'))
    n_npval = sum(1 for s in sources if s == 'MIBiG')
    n_nppred = sum(1 for s in sources if s == 'AS_db')
    description = Counter(names).most_common(1)[0][0] if names else ''
    populated_cols = [c for c in METADATA_COLUMNS if col_has_data.get(c, False)]

    return {
        'n_total': n_total, 'n_char': n_char, 'n_npval': n_npval,
        'n_nppred': n_nppred, 'description': description,
        'populated_cols': populated_cols, 'sf_labels': sf_labels,
        'pfams': ds_pfams, 'cog': ds_cog, 'ec': ds_ec,
    }


# ── Load alistat data ────────────────────────────────────────────────────
alistat = {}
if ALISTAT_FILE.exists():
    with open(ALISTAT_FILE) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= len(header):
                row = dict(zip(header, parts))
                alistat[row.get('dataset', '')] = row
    print(f"  Loaded alistat for {len(alistat)} datasets")

# ── Load PT mapping ──────────────────────────────────────────────────────
pt_to_original = {}
if MAPPING_FILE.exists():
    with open(MAPPING_FILE) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                pt_to_original[parts[0]] = parts[1]

# ── Build dataset entries ────────────────────────────────────────────────
print("\nBuilding dataset entries...")
folders = sorted([p for p in DATASETS_DIR.iterdir()
                  if p.is_dir() and p.name.startswith('PT')])
folders = limit_datasets(folders, cfg, label='PT folders')
print(f"  {len(folders)} PT folders")

sf_groups = defaultdict(list)
stats = {'total': 0, 'no_annot': 0, 'no_tree': 0, 'errors': 0}

for i, folder in enumerate(folders, 1):
    pt_id = folder.name
    stats['total'] += 1

    try:
        annot_file = folder / f'{pt_id}_annotations.tsv'
        tree_file = folder / f'{pt_id}.nwk'
        trimmed_fa = folder / f'{pt_id}_al.fa'
        sequences_fa = folder / f'{pt_id}.fasta'

        tree_path = f'datasets/{pt_id}/{pt_id}.nwk' if tree_file.exists() else ''
        if not tree_file.exists():
            stats['no_tree'] += 1

        # Alistat from original name
        original_name = pt_to_original.get(pt_id, '')
        ali = alistat.get(original_name, {})

        # Parse annotation
        if annot_file.exists():
            ann = parse_annotation(annot_file)
        else:
            stats['no_annot'] += 1
            ann = {
                'n_total': 0, 'n_char': 0, 'n_npval': 0, 'n_nppred': 0,
                'description': '', 'populated_cols': [], 'sf_labels': set(),
                'pfams': set(), 'cog': set(), 'ec': set(),
            }

        # Superfamily combo
        combo = ' | '.join(sorted(ann['sf_labels'])) if ann['sf_labels'] else 'Unknown'

        entry = {
            "name": pt_id,
            "id": pt_id,
            "description": ann['description'],
            "tree": tree_path,
            "tree_model": "LG",
            "tree_rooting_method": "",
            "metadata": f'datasets/{pt_id}/{pt_id}_annotations.tsv',
            "metadata_columns": ann['populated_cols'],
            "N_proteins": ann['n_total'],
            "N_characterized": ann['n_char'],
            "N_np_val": ann['n_npval'],
            "N_np_pred": ann['n_nppred'],
            "trimmed_length": int(float(ali.get('trimmed_length', 0))) if ali.get('trimmed_length') else 0,
            "retention": float(ali.get('retention', 0)) if ali.get('retention') else 0.0,
            "Ca": float(ali.get('Ca', 0)) if ali.get('Ca') else 0.0,
            "alignment": f'datasets/{pt_id}/{pt_id}_al.fa',
            "sequences": f'datasets/{pt_id}/{pt_id}.fasta',
            "source": "automatic",
            "data_type": "protein",
            "reviewed": "yes",
            "PFAMs": sorted(ann['pfams']),
            "COG_category": sorted(ann['cog']),
            "EC": sorted(ann['ec']),
        }

        sf_groups[combo].append(entry)

        if i <= 5 or i % 5000 == 0:
            print(f"  [{i}/{len(folders)}] {pt_id}  N={ann['n_total']}  SF={combo[:60]}")

    except Exception as e:
        stats['errors'] += 1
        print(f"  [{i}] ERROR {pt_id}: {e}")

# ── Update rooting method from summary ──────────────────────────────────
rooting_summary = outdir / '13_rooting_summary.tsv'
if rooting_summary.exists():
    rooting_methods = {}
    with open(rooting_summary) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                row = dict(zip(header, parts))
                rooting_methods[row.get('dataset_id', '')] = row.get('method', '')

    for sf_datasets in sf_groups.values():
        for entry in sf_datasets:
            m = rooting_methods.get(entry['id'], '')
            entry['tree_rooting_method'] = m

# ── Build JSON ───────────────────────────────────────────────────────────
print(f"\nSuperfamily groups: {len(sf_groups)}")

# Load curated superfamilies if available
curated_sfs = []
if CURATED_JSON:
    with open(CURATED_JSON) as f:
        curated = json.load(f)
    curated_sfs = curated.get('superfamilies', [])
    n_cur = sum(len(sf['datasets']) for sf in curated_sfs)
    print(f"  Curated: {len(curated_sfs)} superfamilies, {n_cur} datasets")

auto_sfs = []
for sf_idx, combo_name in enumerate(sorted(sf_groups.keys()), 1):
    datasets = sf_groups[combo_name]
    datasets.sort(key=lambda d: d['id'])
    hmm_list = [] if combo_name == 'Unknown' else [
        s.strip() for s in combo_name.split(' | ')]
    auto_sfs.append({
        "name": f"Superfamily_{sf_idx}",
        "hmm_name": hmm_list,
        "datasets": datasets,
    })

combined = {"superfamilies": curated_sfs + auto_sfs}
n_auto = sum(len(sf['datasets']) for sf in auto_sfs)
print(f"  Automatic: {len(auto_sfs)} superfamilies, {n_auto} datasets")

with open(OUTPUT_JSON, 'w') as f:
    json.dump(combined, f, indent=2)
print(f"  Saved: {OUTPUT_JSON}")

# ── Build SQLite ─────────────────────────────────────────────────────────
print("\nBuilding SQLite database...")

CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS datasets (
    id                    INTEGER PRIMARY KEY AUTOINCREMENT,
    dataset_id            TEXT UNIQUE NOT NULL,
    superfamily_name      TEXT,
    superfamily_hmm_names TEXT,
    name                  TEXT,
    description           TEXT,
    tree                  TEXT,
    tree_model            TEXT,
    tree_rooting_method   TEXT,
    metadata              TEXT,
    metadata_columns      TEXT,
    alignment             TEXT,
    sequences             TEXT,
    source                TEXT,
    data_type             TEXT,
    reviewed              TEXT,
    N_proteins            INTEGER,
    N_characterized       INTEGER,
    N_np_val              INTEGER,
    N_np_pred             INTEGER,
    trimmed_length        REAL,
    retention             REAL,
    Ca                    REAL,
    PFAMs                 TEXT,
    COG_category          TEXT,
    EC                    TEXT,
    citation_authors      TEXT,
    citation_doi          TEXT,
    created_at            TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
"""

INSERT_ROW = """
INSERT OR REPLACE INTO datasets
    (dataset_id, superfamily_name, superfamily_hmm_names,
     name, description, tree, tree_model, tree_rooting_method,
     metadata, metadata_columns, alignment, sequences,
     source, data_type, reviewed,
     N_proteins, N_characterized, N_np_val, N_np_pred,
     trimmed_length, retention, Ca,
     PFAMs, COG_category, EC,
     citation_authors, citation_doi)
VALUES
    (:dataset_id, :superfamily_name, :superfamily_hmm_names,
     :name, :description, :tree, :tree_model, :tree_rooting_method,
     :metadata, :metadata_columns, :alignment, :sequences,
     :source, :data_type, :reviewed,
     :N_proteins, :N_characterized, :N_np_val, :N_np_pred,
     :trimmed_length, :retention, :Ca,
     :PFAMs, :COG_category, :EC,
     :citation_authors, :citation_doi);
"""

if OUTPUT_DB.exists():
    OUTPUT_DB.unlink()

con = sqlite3.connect(OUTPUT_DB)
cur = con.cursor()
cur.executescript(CREATE_TABLE)

n_inserted = 0
for sf in combined['superfamilies']:
    sf_name = sf.get('name', '')
    sf_hmm = json.dumps(sf.get('hmm_name', [])) if sf.get('hmm_name') else None

    for ds in sf.get('datasets', []):
        cite = ds.get('cite', {}) or {}
        raw_name = cite.get('name')
        citation_authors = json.dumps(raw_name) if raw_name else None
        citation_doi = cite.get('doi') or None

        row = {
            'dataset_id': ds.get('id'),
            'superfamily_name': sf_name,
            'superfamily_hmm_names': sf_hmm,
            'name': ds.get('name'),
            'description': ds.get('description'),
            'tree': ds.get('tree'),
            'tree_model': ds.get('tree_model'),
            'tree_rooting_method': ds.get('tree_rooting_method', ''),
            'metadata': ds.get('metadata'),
            'metadata_columns': json.dumps(ds.get('metadata_columns')),
            'alignment': ds.get('alignment'),
            'sequences': ds.get('sequences'),
            'source': ds.get('source'),
            'data_type': ds.get('data_type'),
            'reviewed': ds.get('reviewed'),
            'N_proteins': ds.get('N_proteins'),
            'N_characterized': ds.get('N_characterized'),
            'N_np_val': ds.get('N_np_val'),
            'N_np_pred': ds.get('N_np_pred'),
            'trimmed_length': ds.get('trimmed_length'),
            'retention': ds.get('retention'),
            'Ca': ds.get('Ca'),
            'PFAMs': json.dumps(ds.get('PFAMs', [])),
            'COG_category': json.dumps(ds.get('COG_category', [])),
            'EC': json.dumps(ds.get('EC', [])),
            'citation_authors': citation_authors,
            'citation_doi': citation_doi,
        }
        if row['dataset_id']:
            cur.execute(INSERT_ROW, row)
            n_inserted += 1

con.commit()

total = cur.execute("SELECT COUNT(*) FROM datasets").fetchone()[0]
auto = cur.execute("SELECT COUNT(*) FROM datasets WHERE source='automatic'").fetchone()[0]
cura = cur.execute("SELECT COUNT(*) FROM datasets WHERE source='curated'").fetchone()[0]
print(f"  Rows: {n_inserted}  (auto={auto}, curated={cura})")
print(f"  Saved: {OUTPUT_DB}")
con.close()

print(f"\nDone.")
print(f"  Total datasets:      {stats['total']}")
print(f"  Missing annotations: {stats['no_annot']}")
print(f"  Missing trees:       {stats['no_tree']}")
print(f"  Errors:              {stats['errors']}")
