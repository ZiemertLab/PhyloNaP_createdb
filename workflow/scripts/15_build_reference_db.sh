#!/usr/bin/env bash
# Step 15: Build MMseqs2 reference database for sequence searching.
#
# Dereplicates all merged sequences at 70% identity / 80% coverage,
# then creates a searchable MMseqs2 database.
#
# Output:
#   results/15_reference_db/
#     dereplicated.fasta          Dereplicated sequences
#     mmseqs_db/                  MMseqs2 indexed database
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

echo "========================================================================"
echo "STEP 15: BUILD REFERENCE DATABASE"
echo "========================================================================"
echo "  Input:         $INPUT_FASTA"
echo "  Output:        $REFDB_DIR"
echo "  Dereplication: ${DEDUP_ID} id / ${DEDUP_COV} cov"
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

    # easy-cluster produces dedup_rep_seq.fasta
    mv "$REFDB_DIR/dedup_rep_seq.fasta" "$DEDUP_FASTA"
    rm -f "$REFDB_DIR"/dedup_all_seqs.fasta "$REFDB_DIR"/dedup_cluster.tsv
    rm -rf "$TMP_DIR"

    N_IN=$(grep -c '^>' "$INPUT_FASTA" || true)
    N_OUT=$(grep -c '^>' "$DEDUP_FASTA" || true)
    echo "  Input sequences:  $N_IN"
    echo "  After dedup:      $N_OUT"
else
    echo "Dereplicated FASTA already exists, skipping."
fi

# ── Append curated sequences (already dereplicated) ──────────────────────
COMBINED_FASTA="$DEDUP_FASTA"
if [ -n "$CURATED_FASTA" ]; then
    if [ ! -f "$CURATED_FASTA" ]; then
        echo "ERROR: curated_fasta specified but file not found: $CURATED_FASTA" >&2
        exit 1
    fi
    COMBINED_FASTA="$REFDB_DIR/combined.fasta"
    if [ ! -f "$COMBINED_FASTA" ]; then
        echo ""
        echo "--- Appending curated sequences ---"
        cat "$DEDUP_FASTA" "$CURATED_FASTA" > "$COMBINED_FASTA"
        N_CURATED=$(grep -c '^>' "$CURATED_FASTA" || true)
        N_COMBINED=$(grep -c '^>' "$COMBINED_FASTA" || true)
        echo "  Curated sequences added: $N_CURATED"
        echo "  Combined total:          $N_COMBINED"
    else
        echo "Combined FASTA already exists, skipping."
    fi
fi

# ── Build MMseqs2 database ───────────────────────────────────────────────
DB_PREFIX="$MMSEQS_DB_DIR/phylonap_refdb"
if [ ! -f "${DB_PREFIX}.dbtype" ]; then
    echo ""
    echo "--- Building MMseqs2 database ---"
    "$MMSEQS" createdb "$COMBINED_FASTA" "$DB_PREFIX"
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
echo "  Dereplicated FASTA: $DEDUP_FASTA"
[ -n "$CURATED_FASTA" ] && echo "  Combined FASTA:     $COMBINED_FASTA"
echo "  MMseqs2 database:   $DB_PREFIX"
