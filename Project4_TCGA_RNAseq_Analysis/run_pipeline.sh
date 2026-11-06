
#!/usr/bin/env bash

set -euo pipefail

# Project 4: TCGA-LUAD NKX2-1 analysis

echo "=================================="
echo " TCGA-LUAD NKX2-1 Analysis"
echo "=================================="

# Step 1: Extract expression
echo "[1/4] Extracting NKX2-1 expression..."

python3 scripts/extract_expression.py \
    --sheet data/gdc_sample_sheet.tsv \
    --datadir data/tcga_data \
    --gene NKX2-1 \
    --out results/NKX2-1_expression.tsv

echo "[2/4] Paired dataset must be prepared."
echo "Run the paired dataset preparation script."

echo "[3/4] Running paired statistical analysis..."

Rscript scripts/analyze_paired.R

echo "[4/4] Generating paired expression plots..."

Rscript scripts/plot_paired.R

echo "=================================="
echo " Pipeline completed successfully"
echo "=================================="
