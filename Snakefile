"""
PhyloNaP Pipeline — Snakemake workflow
======================================

Automated construction of the PhyloNaP enzyme phylogeny database.

Usage:
    snakemake --cores 8 --use-conda
    snakemake --cores 1 -n          # dry-run
    snakemake --cores 8 --until annotate    # run up to annotation step

Steps:
    01  cluster           MMseqs2 clustering of input sequences
    02  filter_clusters   Filter clusters by size and BGC content
    03  align_trim        MAFFT alignment + trimAl trimming
    04  fast_trees        FastTree pre-trees for subclustering
    05  treecluster       TreeCluster subclustering
    06  subclusters       Split, re-align, trim, QC subclusters
    07  filter_seqs       Cascade filtration (Cr, Ca, retention, length)
    08  annotate          Comprehensive annotation merging
    09  deduplicate       Remove identical sequences
    10  filter_datasets   Final BGC-linkage filter
    11  build_trees       FastTree LG+gamma phylogenies
    12  organize          Organize into PT-numbered folders
    13  root_trees        Taxonomy + MAD tree rooting
    14  build_database    JSON + SQLite database generation
    15  reference_db      MMseqs2 reference database
"""

import yaml
from pathlib import Path

# ── Load configuration ────────────────────────────────────────────────────
configfile: "config/config.yaml"

OUTDIR = Path(config["output_dir"])
SCRIPTS = Path("workflow/scripts")


# ── Final target ──────────────────────────────────────────────────────────
rule all:
    input:
        OUTDIR / "14_database" / "db_structure.json",
        OUTDIR / "14_database" / "phylonap.db",
        OUTDIR / "15_reference_db" / "dereplicated.fasta",


# ── Step 01: Cluster ──────────────────────────────────────────────────────
rule cluster:
    output:
        directory(OUTDIR / "01_cluster_results"),
    log:
        OUTDIR / ".logs" / "01_cluster.log",
    shell:
        "bash {SCRIPTS}/01_cluster.sh 2>&1 | tee {log}"


# ── Step 02: Filter clusters ─────────────────────────────────────────────
rule filter_clusters:
    input:
        OUTDIR / "01_cluster_results",
    output:
        directory(OUTDIR / "02_filtered_clusters"),
    log:
        OUTDIR / ".logs" / "02_filter_clusters.log",
    shell:
        "python {SCRIPTS}/02_filter_clusters.py 2>&1 | tee {log}"


# ── Step 03: Align + trim ────────────────────────────────────────────────
rule align_trim:
    input:
        OUTDIR / "02_filtered_clusters",
    output:
        directory(OUTDIR / "03_trimmed_alignments"),
    log:
        OUTDIR / ".logs" / "03_align_trim.log",
    shell:
        "bash {SCRIPTS}/03_align_trim.sh 2>&1 | tee {log}"


# ── Step 04: Fast trees for subclustering ─────────────────────────────────
rule fast_trees:
    input:
        OUTDIR / "03_trimmed_alignments",
    output:
        directory(OUTDIR / "04_fast_trees"),
    log:
        OUTDIR / ".logs" / "04_fast_trees.log",
    shell:
        "bash {SCRIPTS}/04_fast_trees.sh 2>&1 | tee {log}"


# ── Step 05: TreeCluster subclustering ────────────────────────────────────
rule treecluster:
    input:
        OUTDIR / "04_fast_trees",
    output:
        directory(OUTDIR / "05_treecluster"),
    log:
        OUTDIR / ".logs" / "05_treecluster.log",
    shell:
        "bash {SCRIPTS}/05_treecluster.sh 2>&1 | tee {log}"


# ── Step 06: Create subclusters ──────────────────────────────────────────
rule subclusters:
    input:
        tc=OUTDIR / "05_treecluster",
        al=OUTDIR / "03_trimmed_alignments",
    output:
        directory(OUTDIR / "06_subclusters"),
    log:
        OUTDIR / ".logs" / "06_subclusters.log",
    shell:
        "python {SCRIPTS}/06_subclusters.py 2>&1 | tee {log}"


# ── Step 07: Cascade filtration ──────────────────────────────────────────
rule filter_sequences:
    input:
        OUTDIR / "06_subclusters",
    output:
        directory(OUTDIR / "07_filtered_datasets"),
    log:
        OUTDIR / ".logs" / "07_filter_sequences.log",
    shell:
        "python {SCRIPTS}/07_filter_sequences.py 2>&1 | tee {log}"


# ── Step 08: Annotate ────────────────────────────────────────────────────
rule annotate:
    input:
        OUTDIR / "07_filtered_datasets",
    output:
        directory(OUTDIR / "08_annotations"),
    log:
        OUTDIR / ".logs" / "08_annotate.log",
    shell:
        "python {SCRIPTS}/08_annotate.py 2>&1 | tee {log}"


# ── Step 09: Deduplicate ─────────────────────────────────────────────────
rule deduplicate:
    input:
        al=OUTDIR / "07_filtered_datasets",
        ann=OUTDIR / "08_annotations",
    output:
        al_out=directory(OUTDIR / "09_deduplicated_alignments"),
        ann_out=directory(OUTDIR / "09_deduplicated_annotations"),
    log:
        OUTDIR / ".logs" / "09_deduplicate.log",
    shell:
        "python {SCRIPTS}/09_deduplicate.py 2>&1 | tee {log}"


# ── Step 10: Filter datasets ─────────────────────────────────────────────
rule filter_datasets:
    input:
        al=OUTDIR / "09_deduplicated_alignments",
        ann=OUTDIR / "09_deduplicated_annotations",
    output:
        al_out=directory(OUTDIR / "10_final_alignments"),
        ann_out=directory(OUTDIR / "10_final_annotations"),
    log:
        OUTDIR / ".logs" / "10_filter_datasets.log",
    shell:
        "python {SCRIPTS}/10_filter_datasets.py 2>&1 | tee {log}"


# ── Step 11: Build phylogenetic trees ────────────────────────────────────
rule build_trees:
    input:
        OUTDIR / "10_final_alignments",
    output:
        directory(OUTDIR / "11_trees"),
    log:
        OUTDIR / ".logs" / "11_build_trees.log",
    shell:
        "python {SCRIPTS}/11_build_trees.py 2>&1 | tee {log}"


# ── Step 12: Organize datasets ───────────────────────────────────────────
rule organize:
    input:
        al=OUTDIR / "10_final_alignments",
        ann=OUTDIR / "10_final_annotations",
        trees=OUTDIR / "11_trees",
    output:
        directory(OUTDIR / "12_datasets"),
    log:
        OUTDIR / ".logs" / "12_organize.log",
    shell:
        "python {SCRIPTS}/12_organize_datasets.py 2>&1 | tee {log}"


# ── Step 13: Root trees ──────────────────────────────────────────────────
rule root_trees:
    input:
        OUTDIR / "12_datasets",
    output:
        OUTDIR / "13_rooting_summary.tsv",
    log:
        OUTDIR / ".logs" / "13_root_trees.log",
    shell:
        "python {SCRIPTS}/13_root_trees.py 2>&1 | tee {log}"


# ── Step 14: Build database ──────────────────────────────────────────────
rule build_database:
    input:
        ds=OUTDIR / "12_datasets",
        rooting=OUTDIR / "13_rooting_summary.tsv",
    output:
        json=OUTDIR / "14_database" / "db_structure.json",
        db=OUTDIR / "14_database" / "phylonap.db",
    log:
        OUTDIR / ".logs" / "14_build_database.log",
    shell:
        "python {SCRIPTS}/14_build_database.py 2>&1 | tee {log}"


# ── Step 15: Reference database ──────────────────────────────────────────
rule reference_db:
    input:
        OUTDIR / "12_datasets",
    output:
        OUTDIR / "15_reference_db" / "dereplicated.fasta",
    log:
        OUTDIR / ".logs" / "15_reference_db.log",
    shell:
        "bash {SCRIPTS}/15_build_reference_db.sh 2>&1 | tee {log}"
