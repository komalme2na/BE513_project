#!/usr/bin/env bash
# run_pipeline.sh — one-command demo of the whole Project 3 workflow.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p data results

echo "==> 1. generating demo fragments (skip if you downloaded the real data)"
[ -f data/fragments.tsv.gz ] || python3 scripts/make_demo_data.py --out data/fragments.tsv.gz

echo "==> 2. building the count matrices (stdin -> stdout, as the brief asks)"
zcat data/fragments.tsv.gz | python3 scripts/vplot.py --format wide   > results/vplot_wide.tsv
zcat data/fragments.tsv.gz | python3 scripts/vplot.py --format long   > results/vplot_long.tsv
zcat data/fragments.tsv.gz | python3 scripts/vplot.py --format matrix > results/vplot_matrix.tsv

echo "==> 3. plotting"
if command -v gnuplot >/dev/null 2>&1; then
    gnuplot gnuplot/vplot_long.gnuplot
    gnuplot gnuplot/vplot_matrix.gnuplot
    echo "    gnuplot wrote results/vplot_long.png and results/vplot_matrix.png"
else
    echo "    gnuplot not found -> using the matplotlib fallback"
    python3 scripts/plot_vplot.py < results/vplot_wide.tsv --out results/vplot.png
fi

echo "==> done. See results/"
