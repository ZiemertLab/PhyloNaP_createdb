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

INPUT_DIR="$OUTDIR/04_fast_trees"
OUTPUT_DIR="$OUTDIR/05_treecluster"
mkdir -p "$OUTPUT_DIR"

echo "============================================================"
echo "STEP 05: TREE-BASED SUBCLUSTERING"
echo "============================================================"
echo "  Threshold: $THRESHOLD"

count=0
total=$(ls "$INPUT_DIR"/*.tree 2>/dev/null | wc -l | tr -d ' ')

for TREE in "$INPUT_DIR"/*.tree; do
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
done

echo "Done. TreeCluster results in: $OUTPUT_DIR ($count processed)"
