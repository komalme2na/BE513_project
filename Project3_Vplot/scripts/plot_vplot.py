#!/usr/bin/env python3
"""
plot_vplot.py — matplotlib alternative to gnuplot for the V-plot.

Reads the WIDE-form matrix produced by vplot.py (1st column = fragment
length, remaining columns = counts at offsets xmin..xmax) and renders a
heatmap PNG.  Density is shown on a log1p scale so the faint nucleosome arms
are visible next to the intense central footprint.

Usage:
    zcat data/fragments.tsv.gz | python3 scripts/vplot.py --format wide \
        | python3 scripts/plot_vplot.py --out results/vplot.png
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/vplot.png")
    ap.add_argument("--xmin", type=int, default=-500)
    ap.add_argument("--xmax", type=int, default=500)
    args = ap.parse_args()

    lengths, grid = [], []
    for line in sys.stdin:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        lengths.append(int(f[0]))
        grid.append([int(v) for v in f[1:]])
    Z = np.asarray(grid, dtype=float)
    ys = np.asarray(lengths)
    xs = np.arange(args.xmin, args.xmax + 1)

    fig, ax = plt.subplots(figsize=(10, 6))

    im = ax.imshow(
        np.log1p(Z),
        origin="lower",
        aspect="auto",
        extent=[
            args.xmin - 0.5,
            args.xmax + 0.5,
            ys[0] - 0.5,
            ys[-1] + 0.5
        ],
        cmap="magma",
        interpolation="nearest"
    )

    ax.set_xlabel("Fragment midpoint distance from site centre (bp)")
    ax.set_ylabel("Fragment length (bp)")
    ax.set_title("V-plot of cfDNA fragments around CTCF sites")

    fig.colorbar(im, ax=ax, label="log1p(count)")
    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
