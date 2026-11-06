
# Sorting large genomic files — strategies compared

## Why not just `sort file.bed`?
Lexicographic order gives `chr1, chr10, chr11, …, chr2, chrM, chrX, chrY`.
Natural ("version") order via `sort -V` fixes digit runs but hard-codes the
chromosome order. Real reference orders (GRCh37 vs GRCh38, decoys, patches,
viral contigs) are arbitrary — so the robust general solution is exactly what
Project 1 does: **map every contig to a rank from an order file, sort by rank**.

## Memory behaviour
- Python `sorted()` / list.sort(): O(N) RAM. 10M short lines ≈ 1–2 GB. Fine
  until it isn't — the failure mode is the OOM killer with no warning.
- GNU `sort`: external merge sort — spills to `$TMPDIR` automatically. Slower
  but bounded. Tune with `-S` and `-T`.
- Rule of thumb: file size < ~20% of free RAM → Python is fine and simpler.
  Otherwise stream through awk+sort.

## `sort -m` (merge)
Only merge-presorted runs: `sort -m a.sorted b.sorted`. If inputs are unsorted
shuffles, `cat a b | sort` is identical in cost and simpler.

## Reproducibility checklist for pipelines
```bash
export LC_ALL=C          # byte-wise, locale-independent order
set -euo pipefail        # catch failures in any pipe stage
```
