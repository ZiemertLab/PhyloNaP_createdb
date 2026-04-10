#!/bin/bash
# ===========================================================================
# Step 01: Merge input FASTA files and cluster with MMseqs2
# ===========================================================================
# Usage: bash 01_cluster.sh <config.yaml>
# Expects: MITE FASTA dir, AntiSMASH/SwissProt/MIBiG FASTA files in resources/
# Outputs: results/01_cluster_results/
# ===========================================================================
set -euo pipefail

CONFIG="${1:-config/config.yaml}"
OUTDIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output_dir'])")
MITE_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['input_data']['mite_fasta_dir'])")
AS_FASTA=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['input_data']['antismash_fasta'])")
SP_FASTA=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['input_data']['swissprot_fasta'])")
MB_FASTA=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['input_data']['mibig_fasta'])")

SENS=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['params']['mmseqs_sensitivity'])")
MIN_ALN=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['params']['mmseqs_min_aln_len'])")
THREADS=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['params']['mmseqs_threads'])")

MMSEQS=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); t=c.get('tools',{}).get('mmseqs',''); print(t if t else 'mmseqs')")

RESULT_DIR="$OUTDIR/01_cluster_results"
mkdir -p "$RESULT_DIR" "$OUTDIR/tmp"

echo "============================================================"
echo "STEP 01: MERGE INPUT DATA & MMSEQS2 CLUSTERING"
echo "============================================================"

# --- Merge MITE FASTAs into one file ---
MITE_MERGED="$RESULT_DIR/mite_combined.fasta"
echo "Merging MITE FASTA files from $MITE_DIR ..."
cat "$MITE_DIR"/*.fasta > "$MITE_MERGED" 2>/dev/null || cat "$MITE_DIR"/*.fa > "$MITE_MERGED" 2>/dev/null || true

# Clean MITE headers: replace spaces with underscores
sed -i.bak 's/ /_/g' "$MITE_MERGED" && rm -f "${MITE_MERGED}.bak"

# --- Concatenate all sources ---
COMBINED="$RESULT_DIR/combined_initial_data.fasta"
echo "Concatenating all sources..."
cat "$AS_FASTA" "$SP_FASTA" "$MITE_MERGED" "$MB_FASTA" > "$COMBINED"
echo "  Combined FASTA: $(grep -c '^>' "$COMBINED") sequences"

# --- Run MMseqs2 clustering ---
echo ""
echo "Running MMseqs2 easy-cluster..."
echo "  Sensitivity: $SENS"
echo "  Min alignment length: $MIN_ALN"
echo "  Threads: $THREADS"

$MMSEQS easy-cluster \
    "$COMBINED" \
    "$RESULT_DIR/cluster_results" \
    "$OUTDIR/tmp" \
    -s "$SENS" \
    --min-aln-len "$MIN_ALN" \
    --threads "$THREADS"

echo ""
echo "Clustering complete."
echo "  Cluster TSV:  $RESULT_DIR/cluster_results_cluster.tsv"
echo "  Rep seqs:     $RESULT_DIR/cluster_results_rep_seq.fasta"
echo "  All seqs:     $RESULT_DIR/cluster_results_all_seqs.fasta"
echo ""
echo "Done."
