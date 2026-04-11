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
TEST_MODE=$(python3 -c "import os,yaml; c=yaml.safe_load(open('$CONFIG')); print('1' if os.environ.get('PHYLONAP_TEST_MODE','')=='1' or c.get('test_mode') else '0')")
TEST_N=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c.get('test_n_datasets', 10))")

INPUT_DIR="$OUTDIR/03_trimmed_alignments"
OUTPUT_DIR="$OUTDIR/04_fast_trees"
mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "STEP 04: BUILD FAST TREES FOR SUBCLUSTERING"
echo "============================================================"

count=0
total=$(find "$INPUT_DIR" -maxdepth 1 -name '*_trimmed.fasta' | wc -l | tr -d ' ')
echo "  Input alignments: $total"

if [ "$TEST_MODE" = "1" ]; then
    echo "  ** TEST MODE: limiting to $TEST_N datasets **"
    total="$TEST_N"
fi

while IFS= read -r -d '' TRIMMED; do
    BASE=$(basename "$TRIMMED" _trimmed.fasta)
    TREE="$OUTPUT_DIR/${BASE}.tree"
    count=$((count + 1))

    [ -f "$TREE" ] && continue

    $FASTTREE -quiet "$TRIMMED" > "$TREE" 2>/dev/null || echo "  WARN: $BASE failed"

    if [ $((count % 500)) -eq 0 ] || [ "$count" -le 5 ]; then
        echo "  [$count/$total] $BASE"
    fi

    if [ "$TEST_MODE" = "1" ] && [ "$count" -ge "$TEST_N" ]; then
        echo "  Test mode: stopping after $TEST_N datasets"
        break
    fi
done < <(find "$INPUT_DIR" -maxdepth 1 -name '*_trimmed.fasta' -print0 | sort -z)

echo "Done. Fast trees in: $OUTPUT_DIR ($count processed)"
