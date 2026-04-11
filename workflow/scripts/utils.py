#!/usr/bin/env python3
"""
PhyloNaP Pipeline — Shared utilities and configuration loader.

All scripts import paths and parameters from config/config.yaml via this module.
"""

import os
import sys
import yaml
from pathlib import Path


def load_config(config_path=None):
    """Load pipeline configuration from YAML file."""
    if config_path is None:
        # Try environment variable, then default location
        # Pipeline root is three levels up: scripts/ -> workflow/ -> pipeline_root/
        config_path = os.environ.get(
            'PHYLONAP_CONFIG',
            Path(__file__).parent.parent.parent / 'config' / 'config.yaml'
        )
    config_path = Path(config_path)
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
        print("Set PHYLONAP_CONFIG or place config.yaml in config/", file=sys.stderr)
        sys.exit(1)
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    return cfg


def resolve_path(path_str, base_dir=None):
    """Resolve a path from config (absolute or relative to pipeline root)."""
    p = Path(path_str)
    if p.is_absolute():
        return p
    if base_dir is None:
        # Pipeline root is three levels up: scripts/ -> workflow/ -> pipeline_root/
        base_dir = Path(__file__).parent.parent.parent
    return (base_dir / p).resolve()


def get_tool(cfg, tool_name):
    """Get tool path from config, falling back to $PATH."""
    custom = cfg.get('tools', {}).get(tool_name, '')
    return custom if custom else tool_name


def is_test_mode(cfg):
    """Check if pipeline is running in test mode.
    
    Respects both config.yaml and the PHYLONAP_TEST_MODE env var
    (set by Snakemake when using --config test_mode=true).
    """
    env = os.environ.get('PHYLONAP_TEST_MODE', '')
    if env.lower() in ('1', 'true', 'yes'):
        return True
    return bool(cfg.get('test_mode', False))


def test_limit(cfg):
    """Return number of datasets to process in test mode (default 10)."""
    return int(cfg.get('test_n_datasets', 10))


def limit_datasets(file_list, cfg, label="datasets"):
    """In test mode, limit a sorted file list to N items and print a notice."""
    if not is_test_mode(cfg):
        return file_list
    n = test_limit(cfg)
    if len(file_list) > n:
        print(f"  ** TEST MODE: limiting {label} from {len(file_list)} to {n} **")
        return file_list[:n]
    return file_list


def val(v):
    """Clean string value: return '' if NaN/None/empty."""
    if v is None:
        return ''
    import pandas as pd
    if isinstance(v, float) and pd.isna(v):
        return ''
    s = str(v).strip()
    return '' if s.lower() == 'nan' else s


def clean_species_name(species_value):
    """Keep only Genus species from a species string."""
    if not species_value:
        return ''
    import pandas as pd
    if isinstance(species_value, float) and pd.isna(species_value):
        return ''
    s = str(species_value).strip()
    if not s:
        return ''
    words = s.split()
    return ' '.join(words[:2]) if len(words) >= 2 else s


# Source priority: higher = better
SOURCE_PRIORITY = {'MITE': 4, 'MIBiG': 3, 'SwissProt': 2, 'AS_db': 1, '': 0}
SOURCE_ORDER = ['mite', 'mibig', 'swiss', 'as_db']  # iteration order (highest first)

# Annotation output columns
ANNOTATION_COLUMNS = [
    'ID', 'Source', 'MITE_ID', 'Dataset_MITE_IDs', 'Uniprot_ID', 'MIBiG_ID',
    'Cluster',
    'gene_name', 'enzyme_name', 'tailoring', 'description', 'reaction_description',
    'cofactors_organic', 'cofactors_inorganic', 'product', 'Rhea', 'Enzyme_function',
    'Entry_Name', 'ProteinExistence', 'EC_Number', 'PDB_IDs',
    'organism', 'Species', 'domain', 'Superkingdom', 'phylum', 'class', 'order', 'family',
    'BGC_type', 'genome_ID', 'Superfamily', 'All_superfamilies',
    'KEGG_Reaction', 'KEGG_rclass', 'COG_category', 'PFAMs', 'eggnog_EC',
]

# Metadata columns tracked in the database
METADATA_COLUMNS = [
    'Source', 'MITE_ID', 'Dataset_MITE_IDs', 'Uniprot_ID', 'MIBiG_ID', 'mibig_id',
    'Cluster',
    'gene_name', 'enzyme_name', 'tailoring', 'description', 'reaction_description',
    'cofactors_organic', 'cofactors_inorganic', 'product', 'Rhea', 'Enzyme_function',
    'Entry_Name', 'ProteinExistence', 'EC_Number', 'PDB_IDs',
    'organism', 'Species', 'domain', 'Superkingdom', 'phylum', 'class', 'order', 'family',
    'BGC_type', 'genome_ID', 'Superfamily', 'All_superfamilies',
    'KEGG_Reaction', 'KEGG_rclass', 'COG_category', 'PFAMs', 'eggnog_EC',
]
