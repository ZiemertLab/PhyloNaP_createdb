#!/bin/bash
# ===========================================================================
# Step 03: Align and trim each cluster
# ===========================================================================
# For each FASTA in 02_filtered_clusters/:
#   1. Clean non-standard amino acids (J,B,Z,O,U -> X)
#   2. Align with MAFFT (--auto --maxiterate N)
#   3. Trim with trimAl (-automated1)
# Output: results/03_trimmed_alignments/
# ===========================================================================
set -euo pipefail

CONFIG="${1:-config/config.yaml}"
OUTDIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output_dir'])")
MAX_ITER=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['params']['mafft_maxiterate'])")
TRIMAL_MODE=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['params']['trimal_mode'])")

MAFFT=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); t=c.get('tools',{}).get('mafft',''); print(t if t else 'mafft')")
TRIMAL=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); t=c.get('tools',{}).get('trimal',''); print(t if t else 'trimal')")

INPUT_DIR="$OUTDIR/02_filtered_clusters"
ALIGNED_DIR="$OUTDIR/03_aligned"
TRIMMED_DIR="$OUTDIR/03_trimmed_alignments"
mkdir -p "$ALIGNED_DIR" "$TRIMMED_DIR"

echo "============================================================"
echo "STEP 03: ALIGN, CLEAN, AND TRIM"
echo "============================================================"
echo "  MAFFT: $MAFFT (--maxiterate $MAX_ITER)"
echo "  trimAl: $TRIMAL ($TRIMAL_MODE)"

count=0
total=$(ls "$INPUT_DIR"/*.fasta 2>/dev/null | wc -l | tr -d ' ')
echo "  Input clusters: $total"

for FASTA in "$INPUT_DIR"/*.fasta; do
    BASE=$(basename "$FASTA" .fasta)
    CLEANED="$ALIGNED_DIR/${BASE}_cleaned.fasta"
    ALIGNED="$ALIGNED_DIR/${BASE}_aligned.fasta"
    TRIMMED="$TRIMMED_DIR/${BASE}_trimmed.fasta"

    count=$((count + 1))

    # Skip if already done
    if [ -f "$TRIMMED" ]; then
        continue
    fi

    # 1. Replace non-standard amino acids
    sed '/^>/!s/[JjBbZzOoUu]/X/g' "$FASTA" > "$CLEANED"

    # 2. Align with MAFFT
    $MAFFT --auto --maxiterate "$MAX_ITER" --quiet "$CLEANED" > "$ALIGNED" 2>/dev/null

    # 3. Trim with trimAl
    $TRIMAL -in "$ALIGNED" $TRIMAL_MODE -out "$TRIMMED" 2>/dev/null || {
        # If trimAl fails, try without automated1
        $TRIMAL -in "$ALIGNED" -gt 0.5 -out "$TRIMMED" 2>/dev/null || cp "$ALIGNED" "$TRIMMED"
    }

    # Clean up intermediate
    rm -f "$CLEANED"

    if [ $((count % 500)) -eq 0 ] || [ "$count" -le 5 ]; then
        echo "  [$count/$total] $BASE done"
    fi
done

echo ""
echo "Done. Trimmed alignments in: $TRIMMED_DIR"
echo "  Total processed: $count"
