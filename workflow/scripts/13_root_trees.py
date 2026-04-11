#!/usr/bin/env python3
"""
Step 13: Root phylogenetic trees using taxonomy-based outgroup detection,
         with MAD rooting fallback via an R helper script.

Algorithm (bipartition-based taxonomy rooting):
  1. Enumerate ALL internal edges (bipartitions) in each tree.
  2. For each edge, evaluate BOTH sides at each taxonomy level
     (Superkingdom -> phylum -> class -> order -> Genus -> Species)
     in two modes:

     Mode A - Pure outgroup (single minority taxon):
       - >=2 leaves on the proposed outgroup side
       - 100% taxonomic purity among annotated outgroup leaves
       - >=70% annotation coverage within the outgroup
       - Outgroup taxon differs from the dominant taxon overall
       - Tiered stem length: Superkingdom/phylum: any positive stem;
         class/order/Genus/Species: stem > median(branch lengths)
       - Parent-level scoping (outgroup within dominant parent taxon)

     Mode B - Dominant monophyly:
       - Outgroup has ZERO annotated members of the dominant taxon
       - Ingroup is 100% pure for the dominant taxon (among annotated)
       - Both sides have >=70% annotation coverage
       - Same tiered stem length check

  3. Candidate ranking: higher taxonomy level > larger outgroup >
     pure over monophyly > longer stem.

  4. Guards: dominant taxa cascade, sub-level scoping,
     low-level annotation coverage, single-genus species guard.

  5. Fallback: MAD rooting via 13_root_trees_mad.R.

Output:
  results/12_datasets/PT######/PT######.nwk            (rooted tree)
  results/12_datasets/PT######/PT######_unrooted.nwk   (preserved)
  results/13_rooting_summary.tsv
"""

import csv
import re
import statistics
import subprocess
import sys
from collections import Counter
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade as BPClade

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_config, resolve_path, get_tool, limit_datasets

cfg = load_config()
outdir = resolve_path(cfg['output_dir'])

DATASETS_DIR = outdir / '12_datasets'
SUMMARY_FILE = outdir / '13_rooting_summary.tsv'
MAD_SCRIPT = Path(__file__).parent / '13_root_trees_mad.R'

# Taxonomy levels and thresholds
TAX_COLUMN_NAMES = ['Superkingdom', 'phylum', 'class', 'order', 'Species']
TAX_LEVELS = ['Superkingdom', 'phylum', 'class', 'order', 'Genus', 'Species']
LENIENT_STEM_LEVELS = {'Superkingdom', 'phylum'}
HIGH_ANNOT_LEVELS = {'order', 'Genus', 'Species'}

# Read rooting config
params = cfg.get('params', {})
USE_TAXONOMY = params.get('taxonomy_rooting', True)
USE_MAD_FALLBACK = params.get('mad_rooting_fallback', True)
MIN_OUTGROUP_SIZE = params.get('taxonomy_min_outgroup', 2)
MIN_ANNOTATION_COVERAGE = params.get('taxonomy_min_coverage', 0.70)
MIN_GLOBAL_ANNOT_COVERAGE = 0.50
SINGLE_GENUS_STEM_FACTOR = 20

print("=" * 80)
print("STEP 13: ROOT PHYLOGENETIC TREES")
print("=" * 80)

# ── Helpers ──────────────────────────────────────────────────────────────────

def load_annotations(annot_file):
    """Read annotation TSV -> dict: seq_id -> {level: value}."""
    annot = {}
    with open(annot_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        col_indices = {}
        for col_name in TAX_COLUMN_NAMES:
            try:
                col_indices[col_name] = header.index(col_name)
            except ValueError:
                col_indices[col_name] = None

        for row in reader:
            if not row:
                continue
            seq_id = row[0]
            tax = {}
            for col_name in TAX_COLUMN_NAMES:
                idx = col_indices.get(col_name)
                if idx is not None and idx < len(row):
                    val = row[idx].strip()
                    tax[col_name] = val if val else None
                else:
                    tax[col_name] = None

            species = tax.get('Species')
            if species:
                words = species.split()
                genus = words[0].strip('[]')
                tax['Genus'] = genus if genus else None
                if len(words) < 2 or words[-1] in ('sp.', 'sp'):
                    tax['Species'] = None
            else:
                tax['Genus'] = None
            annot[seq_id] = tax
    return annot


def get_leaf_names(clade):
    return {t.name for t in clade.get_terminals() if t.name}


def compute_median_branch_length(tree):
    bls = [c.branch_length for c in tree.find_clades()
           if c.branch_length is not None and c.branch_length > 0]
    return statistics.median(bls) if bls else 0.0


def determine_dominant_taxa(annot, all_leaves):
    dominant = {}
    scope = all_leaves
    for level in TAX_LEVELS:
        counts = Counter()
        for leaf in scope:
            if leaf in annot and annot[leaf].get(level) is not None:
                counts[annot[leaf][level]] += 1
        if counts:
            dominant[level] = counts.most_common(1)[0][0]
            annotated_count = sum(counts.values())
            coverage = annotated_count / len(scope) if scope else 0
            if coverage >= 0.50:
                scope = {leaf for leaf in scope
                         if leaf in annot
                         and annot[leaf].get(level) == dominant[level]}
        else:
            dominant[level] = None
            break
    return dominant


def evaluate_outgroup(outgroup_leaves, ingroup_leaves,
                      annot, level, dominant_ingroup, median_bl, stem):
    n_out = len(outgroup_leaves)
    if n_out < MIN_OUTGROUP_SIZE or len(ingroup_leaves) < 1:
        return None

    annotated = [l for l in outgroup_leaves
                 if l in annot and annot[l].get(level) is not None]
    coverage = len(annotated) / n_out
    if coverage < MIN_ANNOTATION_COVERAGE:
        return None

    taxa = Counter(annot[l][level] for l in annotated)
    if len(taxa) != 1:
        return None
    outgroup_taxon = list(taxa.keys())[0]
    if outgroup_taxon == dominant_ingroup:
        return None

    if level in LENIENT_STEM_LEVELS:
        if stem <= 0:
            return None
    else:
        if stem <= median_bl:
            return None

    return {
        'valid': True, 'level': level, 'outgroup_taxon': outgroup_taxon,
        'outgroup_size': n_out, 'annotated_pct': coverage, 'stem_length': stem,
    }


def evaluate_dominant_monophyly(outgroup_leaves, ingroup_leaves,
                                annot, level, dominant_taxon, median_bl, stem):
    n_out = len(outgroup_leaves)
    n_in = len(ingroup_leaves)
    if n_out < MIN_OUTGROUP_SIZE or n_in < MIN_OUTGROUP_SIZE:
        return None
    if n_out >= n_in:
        return None

    out_annotated = [l for l in outgroup_leaves
                     if l in annot and annot[l].get(level) is not None]
    out_coverage = len(out_annotated) / n_out
    if out_coverage < MIN_ANNOTATION_COVERAGE:
        return None
    if sum(1 for l in out_annotated if annot[l][level] == dominant_taxon) > 0:
        return None

    in_annotated = [l for l in ingroup_leaves
                    if l in annot and annot[l].get(level) is not None]
    in_coverage = len(in_annotated) / n_in
    if in_coverage < MIN_ANNOTATION_COVERAGE:
        return None
    if sum(1 for l in in_annotated if annot[l][level] == dominant_taxon) != len(in_annotated):
        return None

    if level in LENIENT_STEM_LEVELS:
        if stem <= 0:
            return None
    else:
        if stem <= median_bl:
            return None

    return {
        'valid': True, 'level': level, 'outgroup_taxon': f'non-{dominant_taxon}',
        'outgroup_size': n_out, 'annotated_pct': out_coverage, 'stem_length': stem,
    }


def fix_rooted_tree(tree, best):
    """Fix trifurcations and balance root branch lengths."""
    # Fix trifurcation
    if len(tree.root.clades) > 2:
        outgroup_leaf_names = get_leaf_names(best['root_clade'])
        outgroup_child = None
        ingroup_children = []
        for child in tree.root.clades:
            child_leaves = {t.name for t in child.get_terminals()}
            if child_leaves <= outgroup_leaf_names and outgroup_child is None:
                outgroup_child = child
            else:
                ingroup_children.append(child)
        if outgroup_child is None:
            outgroup_child = tree.root.clades[0]
            ingroup_children = tree.root.clades[1:]
        if len(ingroup_children) > 1:
            new_ingroup = BPClade(branch_length=0.0, clades=ingroup_children)
            tree.root.clades = [outgroup_child, new_ingroup]

    # Balance root branch lengths
    if len(tree.root.clades) == 2:
        c1, c2 = tree.root.clades
        bl1 = c1.branch_length if c1.branch_length else 0.0
        bl2 = c2.branch_length if c2.branch_length else 0.0
        total = bl1 + bl2
        if total > 0:
            c1.branch_length = total / 2
            c2.branch_length = total / 2
        else:
            c1.branch_length = 1e-6
            c2.branch_length = 1e-6

    # Clean root node
    tree.root.branch_length = None
    tree.root.confidence = None
    tree.root.name = None


def find_outgroup_and_root(tree_file, annot_file):
    """Full taxonomy rooting pipeline for one dataset."""
    result_info = {
        'dataset_id': '', 'method': 'skipped', 'level': '',
        'outgroup_taxon': '', 'outgroup_size': 0,
        'outgroup_annotated_pct': 0.0, 'ingroup_size': 0,
        'stem_length': 0.0, 'median_branch': 0.0,
        'total_leaves': 0, 'reason': '',
    }

    try:
        tree = Phylo.read(tree_file, 'newick')
    except Exception as e:
        result_info['reason'] = f'tree parse error: {e}'
        return None, result_info

    all_leaves = get_leaf_names(tree.root)
    result_info['total_leaves'] = len(all_leaves)
    if len(all_leaves) < 4:
        result_info['reason'] = f'too few leaves ({len(all_leaves)})'
        return None, result_info

    try:
        annot = load_annotations(annot_file)
    except Exception as e:
        result_info['reason'] = f'annotation parse error: {e}'
        return None, result_info

    median_bl = compute_median_branch_length(tree)
    result_info['median_branch'] = median_bl
    if median_bl == 0:
        result_info['reason'] = 'all branch lengths zero'
        return None, result_info

    dominant_taxa = determine_dominant_taxa(annot, all_leaves)

    # Cache leaf sets for internal clades
    clade_leaf_cache = {}
    for clade in tree.find_clades(order='postorder'):
        if not clade.is_terminal():
            clade_leaf_cache[id(clade)] = get_leaf_names(clade)

    # Global annotation coverage per level
    global_annot_coverage = {}
    for level in TAX_LEVELS:
        ann = [l for l in all_leaves
               if l in annot and annot[l].get(level) is not None]
        global_annot_coverage[level] = len(ann) / len(all_leaves)

    # Dominant taxon share per level
    dominant_share = {}
    for level in TAX_LEVELS:
        dom = dominant_taxa.get(level)
        if dom is None:
            dominant_share[level] = 0.0
        else:
            ann_at = [l for l in all_leaves
                      if l in annot and annot[l].get(level) is not None]
            dominant_share[level] = (
                sum(1 for l in ann_at if annot[l][level] == dom)
                / len(ann_at) if ann_at else 0.0
            )

    # Genus diversity check
    genus_counts = Counter()
    for l in all_leaves:
        if l in annot and annot[l].get('Genus') is not None:
            genus_counts[annot[l]['Genus']] += 1
    single_genus = (len(genus_counts) <= 1)

    # Higher-level uniformity check
    HIGHER_LEVELS = ['Superkingdom', 'phylum', 'class', 'order']
    higher_levels_uniform = True
    for hl in HIGHER_LEVELS:
        ann_at = [l for l in all_leaves
                  if l in annot and annot[l].get(hl) is not None]
        if ann_at:
            hl_coverage = len(ann_at) / len(all_leaves)
            if hl_coverage < MIN_GLOBAL_ANNOT_COVERAGE:
                continue
            taxa_at = set(annot[l][hl] for l in ann_at)
            if len(taxa_at) > 1:
                higher_levels_uniform = False
                break

    # Evaluate all bipartitions
    candidates = []
    for clade in tree.find_clades(order='postorder'):
        if clade.is_terminal():
            continue

        clade_leaves = clade_leaf_cache[id(clade)]
        complement_leaves = all_leaves - clade_leaves
        stem = clade.branch_length if clade.branch_length is not None else 0.0

        for outgroup_leaves, ingroup_leaves in [
            (clade_leaves, complement_leaves),
            (complement_leaves, clade_leaves),
        ]:
            for level in TAX_LEVELS:
                dom = dominant_taxa.get(level)
                if dom is None:
                    continue

                if level in ('Genus', 'Species') and not higher_levels_uniform:
                    continue
                if level in HIGH_ANNOT_LEVELS:
                    if global_annot_coverage.get(level, 0) < MIN_GLOBAL_ANNOT_COVERAGE:
                        continue
                if level in ('Genus', 'Species'):
                    if dominant_share.get(level, 0) < 0.33:
                        continue

                species_blocked = (
                    level == 'Species' and single_genus
                    and stem <= median_bl * SINGLE_GENUS_STEM_FACTOR
                )

                # Mode A: Pure outgroup
                if not species_blocked:
                    parent_ok = True
                    level_idx = TAX_LEVELS.index(level)
                    if level_idx > 0:
                        parent_level = TAX_LEVELS[level_idx - 1]
                        parent_dom = dominant_taxa.get(parent_level)
                        if parent_dom is None:
                            parent_ok = False
                        else:
                            annotated_parent = [
                                l for l in outgroup_leaves
                                if l in annot and annot[l].get(parent_level) is not None
                            ]
                            if annotated_parent:
                                parent_taxa = set(annot[l][parent_level]
                                                  for l in annotated_parent)
                                if parent_taxa != {parent_dom}:
                                    parent_ok = False

                    if parent_ok:
                        result = evaluate_outgroup(
                            outgroup_leaves, ingroup_leaves,
                            annot, level, dom, median_bl, stem)
                        if result is not None:
                            result['root_clade'] = clade
                            result['eval_mode'] = 'pure'
                            candidates.append(result)

                # Mode B: Dominant monophyly (blocked at Species level)
                monophyly_blocked = level == 'Species' or species_blocked
                if not monophyly_blocked:
                    result_mono = evaluate_dominant_monophyly(
                        outgroup_leaves, ingroup_leaves,
                        annot, level, dom, median_bl, stem)
                    if result_mono is not None:
                        result_mono['root_clade'] = clade
                        result_mono['eval_mode'] = 'monophyly'
                        candidates.append(result_mono)

    if not candidates:
        reasons = []
        for level in TAX_LEVELS:
            dom = dominant_taxa.get(level)
            if dom is None:
                reasons.append(f'{level}: no annotations')
                break
            scope = all_leaves
            for prev_idx in range(TAX_LEVELS.index(level)):
                prev_level = TAX_LEVELS[prev_idx]
                prev_dom = dominant_taxa.get(prev_level)
                if prev_dom:
                    scope = {l for l in scope if l in annot
                             and annot[l].get(prev_level) == prev_dom}
            counts = Counter()
            for l in scope:
                if l in annot and annot[l].get(level) is not None:
                    counts[annot[l][level]] += 1
            if len(counts) < 2:
                reasons.append(f'{level}: single taxon ({dom})')
            else:
                minorities = {k: v for k, v in counts.items() if k != dom}
                reasons.append(f'{level}: no valid clade ({dict(minorities)})')
        result_info['reason'] = '; '.join(reasons)
        return None, result_info

    # Pick best candidate
    def sort_key(c):
        level_priority = TAX_LEVELS.index(c['level'])
        mode_priority = 0 if c.get('eval_mode') == 'pure' else 1
        return (level_priority, -c['outgroup_size'], mode_priority,
                -c['stem_length'])

    candidates.sort(key=sort_key)
    best = candidates[0]

    try:
        tree.root_with_outgroup(best['root_clade'])
    except Exception as e:
        result_info['reason'] = f'rooting failed: {e}'
        return None, result_info

    fix_rooted_tree(tree, best)

    result_info.update({
        'method': 'taxonomy', 'level': best['level'],
        'outgroup_taxon': best['outgroup_taxon'],
        'outgroup_size': best['outgroup_size'],
        'outgroup_annotated_pct': best['annotated_pct'],
        'ingroup_size': len(all_leaves) - best['outgroup_size'],
        'stem_length': best['stem_length'], 'reason': 'success',
    })

    output = StringIO()
    Phylo.write(tree, output, 'newick')
    nwk_str = output.getvalue().strip()
    nwk_str = re.sub(r'\)[^(;)]*;$', ');', nwk_str)

    return nwk_str, result_info


# ── MAD rooting fallback ─────────────────────────────────────────────────

def try_mad_rooting(tree_file, rooted_file):
    """Attempt MAD rooting via R script. Returns True if successful."""
    if not MAD_SCRIPT.exists():
        return False
    try:
        result = subprocess.run(
            ['Rscript', str(MAD_SCRIPT), str(tree_file), str(rooted_file)],
            capture_output=True, text=True, timeout=300)
        return result.returncode == 0 and rooted_file.exists()
    except Exception:
        return False


# ── Main ─────────────────────────────────────────────────────────────────────

dataset_dirs = sorted([
    d for d in DATASETS_DIR.iterdir()
    if d.is_dir() and d.name.startswith('PT')
])
dataset_dirs = limit_datasets(dataset_dirs, cfg, label='dataset dirs')

print(f"\nFound {len(dataset_dirs)} datasets")
print(f"  Taxonomy levels: {' -> '.join(TAX_LEVELS)}")
print(f"  Min outgroup size: {MIN_OUTGROUP_SIZE}")
print(f"  Min annotation coverage: {MIN_ANNOTATION_COVERAGE:.0%}")

stats = Counter()
summary_fields = [
    'dataset_id', 'method', 'level', 'outgroup_taxon', 'outgroup_size',
    'outgroup_annotated_pct', 'ingroup_size', 'stem_length',
    'median_branch', 'total_leaves', 'reason',
]

with open(SUMMARY_FILE, 'w') as sf:
    sf.write('\t'.join(summary_fields) + '\n')

    for i, ds_dir in enumerate(dataset_dirs, 1):
        ds_id = ds_dir.name
        tree_file = ds_dir / f'{ds_id}_unrooted.nwk'
        annot_file = ds_dir / f'{ds_id}_annotations.tsv'
        rooted_file = ds_dir / f'{ds_id}.nwk'

        if not tree_file.exists():
            stats['missing_tree'] += 1
            continue
        if not annot_file.exists():
            stats['missing_annot'] += 1
            continue

        # Try taxonomy rooting first
        rooted_nwk = None
        info = {'dataset_id': ds_id, 'method': 'skipped', 'level': '',
                'outgroup_taxon': '', 'outgroup_size': 0,
                'outgroup_annotated_pct': 0.0, 'ingroup_size': 0,
                'stem_length': 0.0, 'median_branch': 0.0,
                'total_leaves': 0, 'reason': ''}

        if USE_TAXONOMY:
            rooted_nwk, info = find_outgroup_and_root(tree_file, annot_file)
            info['dataset_id'] = ds_id

        if rooted_nwk is not None:
            with open(rooted_file, 'w') as f:
                f.write(rooted_nwk + '\n')
            stats['taxonomy'] += 1
        elif USE_MAD_FALLBACK:
            # Fallback: MAD rooting only when taxonomy rooting unavailable
            if try_mad_rooting(tree_file, rooted_file):
                info['method'] = 'MAD'
                info['reason'] = 'MAD fallback (taxonomy unavailable)'
                stats['mad'] += 1
            else:
                stats['skipped'] += 1
        else:
            stats['skipped'] += 1

        if i <= 20 or i % 1000 == 0:
            method = info['method']
            if method == 'taxonomy':
                print(f"  [{i}/{len(dataset_dirs)}] {ds_id}: "
                      f"ROOTED by {info['level']} "
                      f"(outgroup={info['outgroup_taxon']}, n={info['outgroup_size']})")
            elif method == 'MAD':
                print(f"  [{i}/{len(dataset_dirs)}] {ds_id}: ROOTED by MAD")
            else:
                print(f"  [{i}/{len(dataset_dirs)}] {ds_id}: "
                      f"SKIP ({info['reason'][:80]})")

        row_vals = [str(info.get(f, '')) for f in summary_fields]
        sf.write('\t'.join(row_vals) + '\n')

        if i % 500 == 0:
            sf.flush()

print(f"\nDone.")
print(f"  Taxonomy rooted: {stats['taxonomy']}")
print(f"  MAD rooted:      {stats['mad']}")
print(f"  Skipped:         {stats['skipped']}")
print(f"  Missing tree:    {stats.get('missing_tree', 0)}")
print(f"  Missing annot:   {stats.get('missing_annot', 0)}")
total = stats['taxonomy'] + stats['mad'] + stats['skipped']
if total > 0:
    pct = (stats['taxonomy'] + stats['mad']) / total * 100
    print(f"  Success rate:    {pct:.1f}%")
print(f"  Summary:         {SUMMARY_FILE}")
