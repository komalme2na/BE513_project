#!/usr/bin/env python3
"""
rescale.py — Rescale a query BED file's fragment-length distribution to a
reference distribution by weighted random sub-sampling (acceptance-rejection).

WHY THIS ALGORITHM
------------------
For every observed fragment length L we know

    q(L)   = fraction of the QUERY sample with length L      (what we have)
    ref(L) = fraction in the REFERENCE with length L         (what we want)

If we keep each length-L fragment with probability

    p(L) = ( ref(L) / q(L) ) / max_L' ( ref(L') / q(L') )

then the post-sampling distribution is exactly ref(L) wherever fragments of
that length exist, and every p(L) is in [0, 1] by construction.  Dividing by
the maximum ratio makes the sampling as *efficient* as possible: at least one
length class is kept in full, so we throw away the minimum amount of data.

This is the classic acceptance-rejection (a.k.a. rejection) sampler where the
query itself is the proposal distribution.  Lengths that are absent from the
query can never appear in the output — those are the points that "miss out"
(e.g. the 300-400 bp gap in the assignment figure).

TWO PASSES OVER THE FILE
------------------------
Pass 1: count fragments per length -> q(L).  Only an integer histogram is
        stored, so RAM stays O(#length bins) no matter how big the BED is.
Pass 2: stream the file again; keep a line if random() < p(length).

Usage:
    python3 scripts/rescale.py \
        --query data/query.bed.gz \
        --reference data/reference.hist \
        --output results/rescaled.bed.gz \
        --seed 1
"""
import argparse
import gzip
import math
import random
import sys
from collections import Counter


def open_maybe_gzip(path: str):
    """Open plain text or .gz transparently (magic-number sniffing would also
    work; here the extension convention is enough for course data)."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_reference(path: str) -> dict[int, float]:
    """Read a two-column histogram: length <TAB/space> normalized frequency.
    Lines starting with '#' are treated as comments/headers and skipped."""
    ref: dict[int, float] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            ref[int(parts[0])] = float(parts[1])
    return ref


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", required=True, help="query BED file (.bed or .bed.gz)")
    ap.add_argument("--reference", required=True, help="reference histogram file")
    ap.add_argument("--output", required=True, help="output BED file (use .gz to gzip)")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed for reproducibility")
    ap.add_argument("--min-len", type=int, default=None,
                    help="optional: drop fragments shorter than this before sampling")
    ap.add_argument("--max-len", type=int, default=None,
                    help="optional: drop fragments longer than this before sampling")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    ref = read_reference(args.reference)

    # Normalise the reference (defensive: real files may not sum to exactly 1)
    total_ref = sum(ref.values())
    if not math.isclose(total_ref, 1.0, rel_tol=1e-3):
        print(f"[info] reference frequencies sum to {total_ref:.6f}; renormalising",
              file=sys.stderr)
        ref = {L: p / total_ref for L, p in ref.items()}

    # ---------------- PASS 1: empirical query distribution q(L) --------------
    counts: Counter[int] = Counter()
    with open_maybe_gzip(args.query) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            length = int(fields[2]) - int(fields[1])   # BED: end - start
            if args.min_len is not None and length < args.min_len:
                continue
            if args.max_len is not None and length > args.max_len:
                continue
            counts[length] += 1

    n_total = sum(counts.values())
    if n_total == 0:
        sys.exit("error: query file contained no usable fragments")

    q = {L: c / n_total for L, c in counts.items()}

    # ---------------- acceptance probabilities -------------------------------
    ratios = {L: ref[L] / q[L] for L in ref if L in q and q[L] > 0}
    if not ratios:
        sys.exit("error: no overlap between reference lengths and query lengths")
    max_ratio = max(ratios.values())
    accept_p = {L: r / max_ratio for L, r in ratios.items()}
    # Lengths not present in `ref` get p = 0 (they are silently dropped).

    expected_keep = sum(q.get(L, 0.0) * p for L, p in accept_p.items())
    uncovered = sum(p for L, p in ref.items() if L not in counts)
    print(f"[info] query fragments : {n_total}", file=sys.stderr)
    print(f"[info] expected keep   : {expected_keep:.1%} of query", file=sys.stderr)
    print(f"[info] ref. mass uncovered by query: {uncovered:.2%} "
          "(these lengths can only 'miss out')", file=sys.stderr)

    # ---------------- PASS 2: stream + accept/reject --------------------------
    kept = 0
    with open_maybe_gzip(args.query) as fin, \
         (gzip.open(args.output, "wt") if args.output.endswith(".gz")
          else open(args.output, "w")) as fout:
        for line in fin:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            length = int(fields[2]) - int(fields[1])
            p = accept_p.get(length, 0.0)
            if p > 0.0 and rng.random() < p:
                fout.write(line)
                kept += 1

    print(f"[info] kept {kept} fragments ({kept / n_total:.1%})", file=sys.stderr)


if __name__ == "__main__":
    main()
