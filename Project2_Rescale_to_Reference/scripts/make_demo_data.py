#!/usr/bin/env python3
"""
make_demo_data.py — generate a small, realistic demo dataset for Project 2.

The real data lives on figshare (see README). This script fabricates a
miniature version with the SAME file formats so the pipeline can be tested
end-to-end without the real files:

  data/query.bed.gz   - gzipped BED3 file of cfDNA fragments whose fragment-
                        length distribution is WRONG on purpose (peaks ~250 bp).
  data/reference.hist - two-column (length, normalized frequency) target
                        distribution with the classic cfDNA peak at ~167 bp.

Usage:
    python3 scripts/make_demo_data.py            # uses defaults
    python3 scripts/make_demo_data.py --n 500000 --outdir data --seed 7
"""
import argparse
import gzip
import random

import numpy as np


def reference_pmf(lengths: np.ndarray) -> np.ndarray:
    """Toy 'reference' cfDNA distribution: sharp nucleosome peak at 167 bp
    plus a small shoulder for di-nucleosomes (~330 bp)."""
    p = np.exp(-0.5 * ((lengths - 167.0) / 25.0) ** 2)          # mono-nucleosome
    p += 0.25 * np.exp(-0.5 * ((lengths - 330.0) / 35.0) ** 2)  # di-nucleosome
    p += 0.05 * np.exp(-0.5 * ((lengths - 80.0) / 20.0) ** 2)   # short TF-protected
    return p / p.sum()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=200_000, help="number of fragments")
    ap.add_argument("--outdir", default="data", help="output directory")
    ap.add_argument("--seed", type=int, default=42, help="random seed")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    lengths = np.arange(30, 601)

    # ---- reference.hist -------------------------------------------------
    ref = reference_pmf(lengths)
    with open(f"{args.outdir}/reference.hist", "w") as fh:
        fh.write("#length\tnormalized_frequency\n")
        for L, p in zip(lengths, ref):
            fh.write(f"{L}\t{p:.8f}\n")

    # ---- query.bed.gz ----------------------------------------------------
    # Query fragments follow a DIFFERENT (flatter, longer) distribution,
    # e.g. a degraded or over-digested sample.
    q = np.exp(-0.5 * ((lengths - 260.0) / 80.0) ** 2)
    q += 0.3 * np.exp(-0.5 * ((lengths - 167.0) / 40.0) ** 2)
    q /= q.sum()
    frag_lengths = rng.choice(lengths, size=args.n, p=q)

    chroms = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    chrom_weights = np.array([len(c) for c in chroms], dtype=float)  # arbitrary
    chrom_weights /= chrom_weights.sum()

    with gzip.open(f"{args.outdir}/query.bed.gz", "wt") as fh:
        for L in frag_lengths:
            chrom = rng.choice(chroms, p=chrom_weights)
            start = rng.integers(1, 100_000_000)
            fh.write(f"{chrom}\t{start}\t{start + int(L)}\n")

    print(f"wrote {args.outdir}/query.bed.gz ({args.n} fragments)")
    print(f"wrote {args.outdir}/reference.hist ({len(lengths)} length bins)")


if __name__ == "__main__":
    main()
