#!/bin/bash
# ===========================================================================
# Step 04: Build fast trees for tree-based subclustering
# ===========================================================================
# Runs FastTree (no gamma, for speed) on trimmed alignments.
# Output: results/04_fast_trees/
# ===========================================================================
set -euo pipefail

CONFIG="${1:-config/config.yaml}"
OUTDIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output_dir'])")
FASTTREE=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); t=c.get('tools',{}).get('fasttree',''); print(t if t else 'fasttree')")

INPUT_DIR="$OUTDIR/03_trimmed_alignments"
OUTPUT_DIR="$OUTDIR/04_fast_trees"
mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "STEP 04: BUILD FAST TREES FOR SUBCLUSTERING"
echo "============================================================"

count=0
total=$(ls "$INPUT_DIR"/*_trimmed.fasta 2>/dev/null | wc -l | tr -d ' ')
echo "  Input alignments: $total"

for TRIMMED in "$INPUT_DIR"/*_trimmed.fasta; do
    BASE=$(basename "$TRIMMED" _trimmed.fasta)
    TREE="$OUTPUT_DIR/${BASE}.tree"
    count=$((count + 1))

    [ -f "$TREE" ] && continue

    $FASTTREE -quiet "$TRIMMED" > "$TREE" 2>/dev/null || echo "  WARN: $BASE failed"

    if [ $((count % 500)) -eq 0 ] || [ "$count" -le 5 ]; then
        echo "  [$count/$total] $BASE"
    fi
done

echo "Done. Fast trees in: $OUTPUT_DIR ($count processed)"
