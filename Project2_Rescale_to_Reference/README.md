# Project 2 — Rescale to Reference

**Task:** given `query.bed.gz` (cfDNA fragments) and
`reference.hist` (a target fragment-length distribution), sub-sample lines of
the query BED so that the *post-sampling fragment-length distribution matches
the reference*.


---


## The algo — acceptance-rejection (rejection) sampling

For each fragment length `L`:

```
q(L)   = fraction of the QUERY with length L        (measured from the file)
ref(L) = fraction in the REFERENCE with length L    (read from reference.hist)
```

Keep each fragment of length `L` with probability:

```
              ref(L) / q(L)
p(L)  =  --------------------------
          max over L' of ref(L') / q(L')
```

**Why this works (the 3-line proof):** after sampling, the probability of a
surviving fragment having length L is proportional to `q(L) · p(L)` ∝ `ref(L)`.
Renormalising gives exactly `ref(L)`. Dividing by the maximum ratio guarantees
every `p(L) ≤ 1` (a probability!), and because the *largest* ratio maps to
`p = 1`, at least one length class is kept in full — i.e. we discard the
**minimum possible amount of data**.

**Why some lengths "miss out":** if the query contains *no* fragments of some
length present in the reference (e.g. the 300–400 bp region in the assignment
figure), sampling cannot create them — acceptance-rejection can only *remove*
mass, never add. The script reports how much reference mass is uncovered.

### Why two passes over the file?

`p(L)` needs the *global* `q(L)`, which is unknown until the whole file has
been seen once. So:

| Pass | Action | Memory |
|------|--------|--------|
| 1 | count fragments per length → `q(L)` | one small histogram |
| 2 | stream again, keep line if `rand() < p(L)` | one line at a time |

Both passes are **streaming** (`gzip.open` in text mode), so a 50 GB BED file
uses the same RAM as a 5 MB one. Never read a genomics file fully into memory
unless you must.


## Code walkthrough — `scripts/rescale.py`

1. **`read_reference()`** — parses the 2-column histogram, skipping `#`
   comments, and *renormalises* defensively (never trust that input sums to 1).
2. **Pass 1** — for each line: `length = int(col3) - int(col2)`, increment a
   `collections.Counter`. Lines starting with `#`, `track`, `browser` are
   skipped (legal BED decorations).
3. **Ratios** — `ref(L)/q(L)` is computed only where `q(L) > 0`; the max ratio
   normalises everything into `[0,1]`. Lengths absent from the reference get
   `p = 0` and vanish.
4. **Pass 2** — stream the file again; `random.random() < p(length)` decides
   keep/drop; survivors are written out, re-gzipped if the name ends in `.gz`.
5. **Reporting** — expected keep rate vs actual, and the uncovered reference
   mass (a diagnostic your examiner will love).

The random generator is **seeded** (`--seed`), so results are exactly
reproducible — a habit that matters in research code.


**Rule of thumb:** if the reshaping depends on a per-bin weight, you want
rejection sampling (or its stratified equivalent); if you just need "fewer
reads", use uniform sampling — but that is *not* this question.

## File formats
### The `.hist` file
- Just two whitespace-separated columns: `length  normalized_frequency`.
- Constraints you should *check, not assume*: frequencies sum to 1, lengths are
  integers, no negative values, monotonic? (no — it's a distribution, not sorted
  probability). Our script renormalises defensively.

## Extensions / new concepts worth adding to your toolbox

1. **argparse + logging + `if __name__ == "__main__"`** — already used here;
   copy this skeleton for every future script.
2. **Reservoir sampling** — when you need *k* random lines from a stream of
   unknown length in ONE pass (Vitter's algorithm R).
3. **Reproducibility**: seed your RNG; record the seed in the output header or
   filename (`rescaled.seed1.bed.gz`).
4. **Parallelise by chromosome**: split BED by chrom → rescale chunks on
   multiple cores (`xargs -P`) → `cat` together. Acceptance probabilities are
   global, so compute them once first.
5. **WGS reality check**: real pipelines keep `mapq` filters, remove
   duplicates (Picard), and use `.bam`+`samtools` instead of BED.
6. Turn the script into a **proper package** (`pyproject.toml`, `pip install -e .`)
   or a **Snakemake/Nextflow rule** — the standard way real labs run pipelines.

