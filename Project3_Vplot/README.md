# Project 3 — V-plot

**Task** write a Python program that takes cfDNA
fragments mapped around binding sites (here: CTCF sites) from **standard
input** and emits a **matrix file** that Gnuplot can draw as a heatmap — the
**V-plot**. Two matrix layouts are accepted: *wide* (1st column = fragment
length, columns 2…1001 = counts at offsets −500…+500) and *long*
(`x  y  z` triplets: offset, length, count).
---


##  The algorithm — a 2-D histogram, nothing more

For each input line (`offset`, `length`):

1. validate & clip to the window (`xmin…xmax`, `ymin…ymax`);
2. `counts[(x, y)] += 1` — a sparse dictionary keyed by the bin;
3. after EOF, walk the **full grid** (every offset × every length, zeros
   included) and print in the requested layout.

Key ideas to internalise:

- **Streaming**: input is read line-by-line from `sys.stdin` — memory is
  O(matrix size), independent of input size. You could pipe a billion
  fragments through it.
- **Sparse accumulation, dense output**: a `dict`/`Counter` avoids
  pre-allocating and never crashes on unexpected coordinates, but the *output*
  must be dense (every bin printed) because gnuplot's grid formats require it.
- **Unix philosophy**: the program is a *filter* (`stdin → stdout`), so it
  composes: `zcat data.gz | python3 vplot.py | gnuplot …`. This is exactly
  what the assignment's "take the data as standard input" is testing.

### Where do `offset` and `length` come from?

Upstream (not part of the deliverable), each cfDNA fragment BED interval is
intersected with CTCF-site intervals (e.g. `bedtools intersect` / `closest`):

```
site centre   = (site_start + site_end) / 2
fragment mid  = (frag_start + frag_end) / 2
offset        = fragment_mid − site_centre        (strand-corrected if needed)
length        = frag_end − frag_start
```

The figshare file (10.6084/m9.figshare.27073276, "DNA molecules mapped to
CTCF binding sites") provides these per-fragment numbers directly; the script
defaults to columns 1 and 2 but `--offset-col/--length-col` handle any layout.

## The three output layout

| layout | produced by | shape | use with |
|---|---|---|---|
| **wide** | `--format wide` | rows = lengths, col 1 = length, cols 2..N+1 = counts at −500…+500 | the format the brief shows; easy to eyeball/`awk`; needs coordinate bookkeeping in gnuplot |
| **long** | `--format long` | `x y z` triplets, x fastest-varying, **blank line between scan-lines** | `plot '…' u 1:2:3 with image` in gnuplot; tidy-data style; also imports straight into pandas/`ggplot` |
| **matrix** | `--format matrix` | gnuplot-native: first row = x coords, first col = y coords | `plot '…' matrix with image` — zero coordinate hacks |

Two details that bite students:

1. **Zeros must be printed.** If you skip empty bins, gnuplot's `with image`
   sees a ragged grid and either refuses or distorts the plot. Long form is
   *not* the same as sparse form.
2. **Blank lines matter.** In 3-column ASCII grid mode gnuplot uses a blank
   line as the *scan-line separator* (one per y value, x fastest-varying).
   `vplot.py --format long` emits exactly this.


## 6. File-format 

- **TSV vs CSV**: genomics defaults to tabs (`\t`) because commas can appear
  inside fields. Always set `-F'\t'` in awk / `delimiter="\t"` in Python —
  and always `split("\t")`, not `split()`, when columns might be empty…
  though for loose whitespace files `split()` is the tolerant choice. The
  script tries tab first, falls back to any-whitespace.
- **`.gz` on stdin**: you don't need to decompress in Python — `zcat` into the
  pipe; the filter never knows the difference.
- **Gnuplot heatmap vocabulary**:
  - `plot 'f' matrix with image` — file is a grid; first row/column can carry
    real coordinates;
  - `plot 'f' u 1:2:3 with image` — 3-column uniform grid, blank-line
    separated scan-lines;
  - `splot` + `set pm3d map` — the older 3-D-flavoured way; works but `with
    image` is faster and cleaner for 2-D density;
  - `set palette defined (…)` + `log10($3+1)` colour scaling — essential here
    because the centre dot is ~100× brighter than the arms (a linear scale
    shows only the dot).
- **pngcairo terminal** (`set terminal pngcairo`) is the modern gnuplot PNG
  driver; the old `png` driver lacks antialiasing and fonts.
- **stdout vs files**: writing results to stdout lets the *caller* decide the
  filename and enables piping; diagnostics go to **stderr**
  (`print(…, file=sys.stderr)`) so they never corrupt the data stream. This
  separation is a core Unix habit.
