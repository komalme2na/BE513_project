#!/usr/bin/env python3
"""
plot_distributions.py — verification plot for Project 2.

Overlays the fragment-length distribution of
  (a) the original query BED,
  (b) the rescaled (sub-sampled) BED,
  (c) the reference histogram,
and reports the total-variation distance before vs after rescaling.

Usage:
    python3 scripts/plot_distributions.py \
        --query data/query.bed.gz \
        --rescaled results/rescaled.bed.gz \
        --reference data/reference.hist \
        --out results/distribution_check.png
"""
import argparse
import gzip
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def length_hist(path: str) -> Counter:
    c: Counter = Counter()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            c[int(f[2]) - int(f[1])] += 1
    return c


def read_reference(path: str):
    xs, ys = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            a, b = line.split()[:2]
            xs.append(int(a)); ys.append(float(b))
    return np.array(xs), np.array(ys)


def tv_distance(p: np.ndarray, q: np.ndarray) -> float:
    """Total-variation distance between two probability vectors on the
    same support: 0 = identical, 1 = completely disjoint."""
    return 0.5 * np.abs(p - q).sum()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--rescaled", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ref_x, ref_y = read_reference(args.reference)
    ref_y = ref_y / ref_y.sum()

    qh = length_hist(args.query)
    rh = length_hist(args.rescaled)

    grid = ref_x
    q_y = np.array([qh.get(x, 0) for x in grid], dtype=float); q_y /= q_y.sum()
    r_y = np.array([rh.get(x, 0) for x in grid], dtype=float); r_y /= r_y.sum()

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(grid, q_y, lw=1.6, label=f"query (before)   TV={tv_distance(q_y, ref_y):.3f}")
    ax.plot(grid, r_y, lw=1.6, label=f"rescaled (after)  TV={tv_distance(r_y, ref_y):.3f}")
    ax.plot(grid, ref_y, lw=2.0, color="green", label="reference")
    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Normalised frequency")
    ax.set_title("Rescale-to-reference: fragment length distributions")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")
    print(f"TV distance before rescaling: {tv_distance(q_y, ref_y):.4f}")
    print(f"TV distance after  rescaling: {tv_distance(r_y, ref_y):.4f}")


if __name__ == "__main__":
    main()
