#!/usr/bin/env bash
# Step 15: Build MMseqs2 reference database for sequence searching.
#
# The MMseqs2 database is built from ALL input sequences so that every
# original sequence ID is preserved and searchable.  Dereplication at
# 70 % identity / 80 % coverage is run separately to produce a compact
# FASTA (dereplicated.fasta) for reference, but it does NOT feed into
# the database – the full sequence set does.
#
# Output:
#   results/15_reference_db/
#     dereplicated.fasta          Compact dereplicated sequences (reference only)
#     mmseqs_db/                  MMseqs2 indexed database (all sequences)
#
# Usage: bash 15_build_reference_db.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG="$SCRIPT_DIR/../../config/config.yaml"

# Read config
read_cfg() { python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print($1)"; }

OUTPUT_DIR="$(read_cfg "c['output_dir']")"
MMSEQS="$(read_cfg "c.get('tools',{}).get('mmseqs') or 'mmseqs'")"
MMSEQS="${MMSEQS:-mmseqs}"
DEDUP_ID="$(read_cfg "c['params']['dereplication_identity']")"
DEDUP_COV="$(read_cfg "c['params']['dereplication_coverage']")"
CURATED_FASTA="$(read_cfg "c.get('curated_fasta') or ''")" || CURATED_FASTA=""

INPUT_FASTA="$OUTPUT_DIR/12_datasets/all_sequences_merged.fasta"
REFDB_DIR="$OUTPUT_DIR/15_reference_db"
DEDUP_FASTA="$REFDB_DIR/dereplicated.fasta"
MMSEQS_DB_DIR="$REFDB_DIR/mmseqs_db"

# The FASTA that goes into the MMseqs2 DB always uses the full sequence set
# (plus curated sequences if provided).  DB_INPUT_FASTA is set below.
DB_INPUT_FASTA="$INPUT_FASTA"

echo "========================================================================"
echo "STEP 15: BUILD REFERENCE DATABASE"
echo "========================================================================"
echo "  Input:         $INPUT_FASTA"
echo "  Output:        $REFDB_DIR"
echo "  Dereplication: ${DEDUP_ID} id / ${DEDUP_COV} cov (compact FASTA only, not used for DB)"
[ -n "$CURATED_FASTA" ] && echo "  Curated FASTA: $CURATED_FASTA" || echo "  Curated FASTA: (not provided)"
echo ""

mkdir -p "$REFDB_DIR" "$MMSEQS_DB_DIR"

# ── Dereplication ─────────────────────────────────────────────────────────
if [ ! -f "$DEDUP_FASTA" ]; then
    echo "--- Dereplicating sequences ---"
    TMP_DIR="$REFDB_DIR/tmp_mmseqs"
    mkdir -p "$TMP_DIR"

    "$MMSEQS" easy-cluster \
        "$INPUT_FASTA" \
        "$REFDB_DIR/dedup" \
        "$TMP_DIR" \
        --min-seq-id "$DEDUP_ID" \
        -c "$DEDUP_COV" \
        --cov-mode 1 \
        --threads "$(nproc 2>/dev/null || echo 4)"

    # easy-cluster produces dedup_rep_seq.fasta, dedup_all_seqs.fasta, dedup_cluster.tsv
    mv "$REFDB_DIR/dedup_rep_seq.fasta"  "$DEDUP_FASTA"
    mv "$REFDB_DIR/dedup_all_seqs.fasta" "$REFDB_DIR/all_seqs.fasta"
    # Keep cluster membership so every original ID can be mapped back to its representative.
    # Place it inside mmseqs_db/ so it is co-located with the database files and is
    # included whenever the DB directory is copied to a deployment location.
    mv "$REFDB_DIR/dedup_cluster.tsv"    "$MMSEQS_DB_DIR/cluster_membership.tsv"
    rm -rf "$TMP_DIR"

    N_IN=$(grep -c '^>' "$INPUT_FASTA" || true)
    N_OUT=$(grep -c '^>' "$DEDUP_FASTA" || true)
    echo "  Input sequences:  $N_IN"
    echo "  After dedup:      $N_OUT"
else
    echo "Dereplicated FASTA already exists, skipping."
fi

# combine curated and automated fasta

# ── Append curated sequences to the full input for the DB ────────────────
# Curated sequences are appended to the FULL input (not the dereplicated
# subset) so all original IDs remain in the database.
if [ -n "$CURATED_FASTA" ]; then
    if [ ! -f "$CURATED_FASTA" ]; then
        echo "ERROR: curated_fasta specified but file not found: $CURATED_FASTA" >&2
        exit 1
    fi
    DB_INPUT_FASTA="$REFDB_DIR/db_input.fasta"
    if [ ! -f "$DB_INPUT_FASTA" ]; then
        echo ""
        echo "--- Appending curated sequences to full input ---"
        cat "$INPUT_FASTA" "$CURATED_FASTA" > "$DB_INPUT_FASTA"
        N_CURATED=$(grep -c '^>' "$CURATED_FASTA" || true)
        N_DB=$(grep -c '^>' "$DB_INPUT_FASTA" || true)
        echo "  Curated sequences added: $N_CURATED"
        echo "  DB input total:          $N_DB"
    else
        echo "DB input FASTA already exists, skipping."
    fi
fi

# ── Build MMseqs2 database ───────────────────────────────────────────────
DB_PREFIX="$MMSEQS_DB_DIR/referenceDB"
if [ ! -f "${DB_PREFIX}.dbtype" ]; then
    echo ""
    echo "--- Building MMseqs2 database ---"
    "$MMSEQS" createdb "$DB_INPUT_FASTA" "$DB_PREFIX"
    echo "--- Creating index ---"
    "$MMSEQS" createindex "$DB_PREFIX" "$MMSEQS_DB_DIR/tmp" \
        --threads "$(nproc 2>/dev/null || echo 4)"
    rm -rf "$MMSEQS_DB_DIR/tmp"
    echo "  Database: $DB_PREFIX"
else
    echo "MMseqs2 database already exists, skipping."
fi

echo ""
echo "========================================================================"
echo "REFERENCE DATABASE BUILD COMPLETE"
echo "========================================================================"
echo "  Dereplicated FASTA:   $DEDUP_FASTA  (compact reference, not the DB source)"
[ -n "$CURATED_FASTA" ] && echo "  DB input FASTA:       $DB_INPUT_FASTA  (full set + curated)"
echo "  MMseqs2 database:     $DB_PREFIX  (all original sequence IDs preserved)"
