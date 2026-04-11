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
TEST_MODE=$(python3 -c "import os,yaml; c=yaml.safe_load(open('$CONFIG')); print('1' if os.environ.get('PHYLONAP_TEST_MODE','')=='1' or c.get('test_mode') else '0')")
TEST_N=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c.get('test_n_datasets', 10))")

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
total=$(find "$INPUT_DIR" -maxdepth 1 -name '*.fasta' | wc -l | tr -d ' ')
echo "  Input clusters: $total"

# In test mode, limit to N datasets
FIND_CMD="find '$INPUT_DIR' -maxdepth 1 -name '*.fasta' -print0 | sort -z"
if [ "$TEST_MODE" = "1" ]; then
    echo "  ** TEST MODE: limiting to $TEST_N datasets **"
    total="$TEST_N"
fi

while IFS= read -r -d '' FASTA; do
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
    # In test mode, stop after N datasets
    if [ "$TEST_MODE" = "1" ] && [ "$count" -ge "$TEST_N" ]; then
        echo "  Test mode: stopping after $TEST_N datasets"
        break
    fi
done < <(find "$INPUT_DIR" -maxdepth 1 -name '*.fasta' -print0 | sort -z)

echo ""
echo "Done. Trimmed alignments in: $TRIMMED_DIR"
echo "  Total processed: $count"
