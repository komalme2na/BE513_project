
# Project 1 — Quick BED Sort

**Task:** merge two (or more) gzipped BED files of one sample, keep only the
chromosomes listed in a selection file, and sort the result **in the order the
chromosomes appear in the selection file** (not alphabetical order!). Output goes
to `sorted_bed_file_per_sample/` with a name that encodes both sample and selection.

```
shuf.a.bed.gz ─┐
               ├─► filter by selection ─► sort by (selection-rank, start, end) ─► X.standard_selection.sorted.bed
shuf.b.bed.gz ─┘
```
## 2. logic


1. **Selection = a ranking, not a set.** If we only *filter*, order is whatever the
   input happened to have.convert the selection file
   into a lookup `chrom → rank` (line 1 of the file = rank 1, etc.).
2. **Sort by a key, not by the line.** We decorate each kept line with
   `(rank, start, end)` and sort numerically by those three fields. Sorting by the
   raw chromosome string would give `chr1, chr10, chr11, … chr2` (lexicographic),
   which is wrong here. This is the classic **Schwartzian transform / decorate-sort-undecorate**
   pattern.
3. **Merge = concatenate-then-sort**, not sort-each-then-merge. Since GNU `sort` /
   Python's `sorted` handle the whole stream, we can simply stream both files
   through one sort. (If the files were each *already* sorted we'd use `sort -m`,
   see alternatives.)



## File formats

### BED (Browser Extensible Data)
- **Plain text, TAB-delimited** (not CSV — tabs matter; Excel will silently mangle it).
- **Zero-based, half-open coordinates**: `start` is 0-based, `end` is exclusive.
  A fragment covering the first 100 bp is `chr1  0  100`. Length = `end − start`.
- **Column 1–3 are mandatory**: `chrom, start, end`. Optional: `name, score, strand, …` up to BED12.
- **Constraints**: `start ≥ 0`, `end > start`, coordinates are integers.
  Files violating this are silently mishandled by most tools — worth a `awk '$2<0 || $3<=$2{bad++} END{print bad+0}'` sanity check.
- No header line (convention), comment lines start with `#`.
- Coordinates **do not have to be sorted** — that's why this project exists.

### Selection file (`.tsv`)
- One chromosome name per line; **order is the specification**.
- Must not contain duplicates (our rank map silently keeps the first occurrence — documented behaviour).
- Chromosomes absent from it (e.g. `chrM`, `chrUn_*`, `chrEBV`) are dropped from output.

### gzip (`.gz`)
- `zcat`/`gzip -dc` decompresses to stdout — files never need unpacking on disk.
- Block-gzip (`bgzip`) files (common in genomics) are also readable by `zcat`, and additionally support random access via `.tbi` indices — irrelevant here but you'll meet them in samtools/tabix workflows.

---

## 5. bucketlist

- **`NR==FNR` idiom**
- **GNU `sort` tuning**
- **Locale trap**
- **Scaling test**
