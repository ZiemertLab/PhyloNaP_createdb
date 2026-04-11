#!/bin/bash
# ===========================================================================
# Step 05: Tree-based subclustering with TreeCluster
# ===========================================================================
# Runs TreeCluster on each fast tree to detect sub-clusters.
# Output: results/05_treecluster/
# ===========================================================================
set -euo pipefail

CONFIG="${1:-config/config.yaml}"
OUTDIR=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output_dir'])")
THRESHOLD=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['params']['treecluster_threshold'])")
TREECLUSTER=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); t=c.get('tools',{}).get('treecluster',''); print(t if t else 'TreeCluster.py')")
TEST_MODE=$(python3 -c "import os,yaml; c=yaml.safe_load(open('$CONFIG')); print('1' if os.environ.get('PHYLONAP_TEST_MODE','')=='1' or c.get('test_mode') else '0')")
TEST_N=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c.get('test_n_datasets', 10))")
INPUT_DIR="$OUTDIR/04_fast_trees"
OUTPUT_DIR="$OUTDIR/05_treecluster"
mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "STEP 05: TREE-BASED SUBCLUSTERING"
echo "============================================================"
echo "  Threshold: $THRESHOLD"

count=0
total=$(find "$INPUT_DIR" -maxdepth 1 -name '*.tree' | wc -l | tr -d ' ')

if [ "$TEST_MODE" = "1" ]; then
    echo "  ** TEST MODE: limiting to $TEST_N datasets **"
    total="$TEST_N"
fi

while IFS= read -r -d '' TREE; do
    BASE=$(basename "$TREE" .tree)
    OUT="$OUTPUT_DIR/${BASE}_clusters.tsv"
    count=$((count + 1))

    [ -f "$OUT" ] && continue

    $TREECLUSTER -i "$TREE" -t "$THRESHOLD" -o "$OUT" 2>/dev/null || {
        echo "  WARN: $BASE failed"
        continue
    }

    if [ $((count % 500)) -eq 0 ] || [ "$count" -le 5 ]; then
        echo "  [$count/$total] $BASE"
    fi

    if [ "$TEST_MODE" = "1" ] && [ "$count" -ge "$TEST_N" ]; then
        echo "  Test mode: stopping after $TEST_N datasets"
        break
    fi
done < <(find "$INPUT_DIR" -maxdepth 1 -name '*.tree' -print0 | sort -z)

echo "Done. TreeCluster results in: $OUTPUT_DIR ($count processed)"
