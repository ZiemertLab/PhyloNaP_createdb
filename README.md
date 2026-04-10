# PhyloNaP Pipeline

Automated construction of the **PhyloNaP** enzyme phylogeny database — a comprehensive resource linking natural product biosynthetic enzymes (tailoring enzymes) with phylogenetic, functional, and genomic context.

## Overview

The pipeline processes protein sequences from multiple annotation sources (MITE, MIBiG, SwissProt, AntiSMASH), clusters them, builds phylogenies, and produces a searchable database with superfamily-grouped datasets.

```
Input FASTA files
      │
      ▼
┌─────────────────────┐
│ 01. MMseqs2 cluster │
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 02. Filter clusters │  (>10 seqs, ≥3 BGC-linked)
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 03. MAFFT + trimAl  │
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 04. Fast pre-trees  │
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 05. TreeCluster     │  Clade-based subclustering
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 06. Split + re-QC   │  Subclusters → align → trim → AliStat
└────────┬────────────┘
         ▼
┌─────────────────────────────┐
│ 07. Cascade filtration      │  Cr ≥ 0.6, Ca > 0.7,
│     (remove bad sequences,  │  retention ≥ 0.75,
│      re-align, quality gate)│  mean length > 70, n ≥ 10
└────────┬────────────────────┘
         ▼
┌─────────────────────┐
│ 08. Annotate        │  MITE > MIBiG > SwissProt > AS_db
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 09. Deduplicate     │  Remove identical sequences
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 10. Final filter    │  ≥1 BGC-linked sequence required
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 11. Build trees     │  FastTree LG + gamma
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 12. Organize        │  PT IDs, merged FASTA
└────────┬────────────┘
         ▼
┌─────────────────────┐
│ 13. Root trees      │  Taxonomy outgroup + MAD fallback
└────────┬────────────┘
         ▼
┌──────────────────────────┐
│ 14. Build database       │  JSON + SQLite, superfamily grouping
└────────┬─────────────────┘
         ▼
┌──────────────────────────┐
│ 15. Reference DB         │  MMseqs2 dereplicated search DB
└──────────────────────────┘
```

## Quick Start

### 1. Create the Conda environment

```bash
conda env create -f environment.yml
conda activate phylonap_createdb
```

### 2. Configure

Edit `config/config.yaml` to set your input data paths and parameters.

Required input files:
- **FASTA sequences** from AntiSMASH (AS_db), MIBiG, SwissProt, and/or MITE
- **Annotation tables** (TSV) for each data source
- **ID mapping** file linking cluster IDs to original accessions

### 3. Run

```bash
# Full pipeline
snakemake --cores 8

# Dry run (preview steps)
snakemake --cores 1 -n

# Run up to a specific step
snakemake --cores 8 --until annotate

# Run a single step
snakemake --cores 4 build_trees
```

## Configuration

All parameters are centralized in [`config/config.yaml`](config/config.yaml):

| Section | Description |
|---------|-------------|
| `output_dir` | Base output directory |
| `input_data` | Paths to FASTA input files |
| `annotations` | Paths to annotation/metadata TSVs |
| `id_mapping` | Sequence ID mapping file |
| `tools` | Paths to external tools (mmseqs, mafft, etc.) |
| `params` | All pipeline thresholds and parameters |

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mmseqs_sensitivity` | 7.5 | MMseqs2 clustering sensitivity |
| `mmseqs_min_aln_len` | 80 | Minimum alignment length for clustering |
| `min_cluster_size` | 10 | Minimum sequences per cluster |
| `min_bgc_sequences` | 3 | Minimum BGC-derived sequences per cluster |
| `treecluster_threshold` | 4.0 | TreeCluster clade distance threshold |
| `min_cr` | 0.6 | Minimum per-column residue score |
| `min_ca` | 0.7 | Minimum alignment column score |
| `min_retention` | 0.75 | Minimum trimAl retention fraction |
| `min_mean_length` | 70 | Minimum mean sequence length after trimming |
| `fasttree_model` | `-lg -gamma` | FastTree substitution model |

## Output Structure

```
results/
├── 01_cluster_results/     MMseqs2 clustering output
├── 02_filtered_clusters/   Per-cluster FASTA files
├── 03_trimmed_alignments/  MAFFT + trimAl output
├── 04_fast_trees/          Pre-trees for subclustering
├── 05_treecluster/         TreeCluster assignments
├── 06_subclusters/         Split subclusters + QC stats
├── 07_filtered_datasets/   Quality-filtered datasets
├── 08_annotations/         Merged annotation tables
├── 09_deduplicated_*/      Deduplicated alignments + annotations
├── 10_final_*/             Final filtered datasets
├── 11_trees/               ML phylogenies
├── 12_datasets/            Organized PT-numbered folders
│   ├── PT000001/
│   │   ├── PT000001.fasta
│   │   ├── PT000001_al.fa
│   │   ├── PT000001.nwk
│   │   └── PT000001_annotations.tsv
│   ├── pt_id_mapping.tsv
│   └── all_sequences_merged.fasta
├── 13_rooting_summary.tsv
├── 14_database/
│   ├── db_structure.json
│   └── phylonap.db
└── 15_reference_db/
    ├── dereplicated.fasta
    └── mmseqs_db/
```

## Dependencies

Managed via Conda (`environment.yml`):

**Bioinformatics tools:**
- MMseqs2 (sequence clustering + reference DB)
- MAFFT (multiple sequence alignment)
- trimAl (alignment trimming)
- FastTree (phylogenetic tree inference)
- AliStat (alignment quality statistics)
- TreeCluster (clade-based tree cutting)

**Python packages:**
- BioPython, pandas, numpy, PyYAML

**R packages** (for MAD rooting):
- ape, phangorn, Rcpp, inline

## Annotation Priority

Sequences are annotated using a priority cascade:

1. **MITE** — Experimentally characterized enzyme reactions
2. **MIBiG** — Validated BGC-associated enzymes
3. **SwissProt** — Curated protein annotations
4. **AS_db** — AntiSMASH-predicted BGC enzymes

When a sequence appears in multiple sources, higher-priority annotations take precedence. MITE_IDs are assigned only by direct sequence match (not propagated across clusters).

## Tree Rooting Algorithm

Trees are rooted using a bipartition-based taxonomy outgroup detection:

1. All internal edges are evaluated as potential rooting positions
2. Both sides of each edge are tested at each taxonomy level (Superkingdom → phylum → class → order → Genus → Species)
3. **Mode A** (pure outgroup): minority taxon forms a pure clade
4. **Mode B** (dominant monophyly): rooting makes the dominant taxon monophyletic
5. Guards: annotation coverage thresholds, stem length checks, single-genus species protection
6. Fallback: MAD (Minimal Ancestor Deviation) rooting via R

## Cascade Filtration

The quality filtration preserves phylogenetic signal through iterative refinement:

1. **Remove bad sequences**: column residue score (Cr) < 0.6
2. **Re-align** remaining sequences (MAFFT)
3. **Re-trim** (trimAl automated1)
4. **Quality gate**: Ca > 0.7, retention ≥ 0.75, mean_length > 70, n ≥ 10

Datasets failing the quality gate are excluded.

## Citation

If you use PhyloNaP in your research, please cite:

> [Citation to be added upon publication]

## License

[License to be specified]
