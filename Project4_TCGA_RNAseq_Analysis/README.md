# Project 4 — TCGA RNA-seq Analysis (NKX2-1 in Lung Adenocarcinoma)

**Task:** using the GDC data bundle (`gdc_sample_sheet.tsv`,
`gdc_manifest.txt`, `tcga_data.tar.gz`), quantify whether the biomarker
**NKX2-1** is differentially expressed between **Primary Tumor (PT)** and
**Solid Tissue Normal (STN)** samples of TCGA-LUAD, and draw a boxplot of
**log2(TPM + 1)** for the two groups.

This repo = working pipeline + study guide: §1 biology, §2 data formats,
§3 pipeline logic, §4 run instructions, §5 the awk one-liner alternative,
§6 statistics choices, §7 pitfalls, §8 extensions, §9 viva questions.


## files — what they are and their constraints

### `gdc_sample_sheet.tsv` — the metadata map
Tab-separated, one row per downloaded **file**. Key columns (resolve them
**by name**, not position — GDC has changed column order before):

| column | meaning |
|---|---|
| `File ID` | UUID = the **directory name** inside the untarred bundle |
| `File Name` | the counts file inside that directory |
| `Case ID` | the patient (TCGA barcode, e.g. `TCGA-05-4244`) |
| `Sample Type` | `Primary Tumor`, `Solid Tissue Normal`, `Metastatic`, `Recurrent Tumor`… — we keep only PT & STN |

Constraint that bites everyone: the sample sheet is the **only** link between
file UUIDs and biology. Lose it and the expression data is anonymous.

### `gdc_manifest.txt` — the download receipt
`id  filename  md5  size  state` — used by `gdc-client` to download/verify
files. The **md5** column lets you verify integrity:
`md5sum -c` style checks. For this project it is informational, but in real
work always verify checksums before analysis.

### `tcga_data.tar.gz` → `*.rna_seq.augmented_star_gene_counts.tsv`
One directory **per File ID UUID**, each holding one counts file. Layout of
each TSV (GDC "STAR - Counts" workflow, GENCODE v36, GRCh38):

```
# gene-model: GENCODE v36            ← '#' comment lines (STAR metadata)
# workflow-type: STAR - Counts
gene_id  gene_name  gene_type  unstranded  stranded_first  stranded_second  tpm_unstranded  fpkm_unstranded  fpkm_uq_unstranded
N_unmapped  N_unmapped  N_unmapped  …        ← 4 QC rows (not genes!)
ENSG00000136352.15  NKX2-1  protein_coding  81234  40610  40624  53.1718  …
```

Properties/constraints to remember:
- **`tpm_unstranded` is the 7th column** — the assignment's hint
  (`53.1718` for the example STN file). Robust code finds it **by header
  name** anyway (columns can move; names don't).
- **`gene_id` carries an Ensembl version suffix** (`.15`): strip it before
  joining to other databases (`ENSG00000136352.15 → ENSG00000136352`).
- The 4 `N_*` rows are QC counters, not genes — skip or you'll "find" weird
  genes. Matching by exact `gene_name == "NKX2-1"` avoids them naturally.
- **TPM** (Transcripts Per Million) = length-normalised, library-size-scaled
  abundance; **within-sample comparable**. FPKM is the older cousin;
  `unstranded` = raw read counts (needed for DESeq2-style methods).
- The tar is **gzipped tar** (`tar -zxvf`): archive of many files + DEFLATE
  compression. `z` = gzip, `x` = extract, `v` = verbose, `f` = file.

## The pipeline — logic walkthrough

```
sample sheet ──► keep PT & STN rows ──► {File ID → (Case, Sample Type)}
                                              │
tcga_data.tar.gz ──► tar -zxvf ──► <uuid>/<uuid>.rna_seq...tsv
                                              │
                       for each kept File ID: find counts file
                       find row gene_name == NKX2-1
                       take tpm_unstranded (by header name)
                                              │
                                              ▼
                          expression table (file, case, type, TPM, log2(TPM+1))
                                              │
                                              ▼
                     boxplot PT vs STN  +  Mann–Whitney U p-value
```

`scripts/extract_expression.py` does the extraction. Design choices that make
it *correct by construction* (contrast with the awk quickie in §5):

1. **`csv.DictReader`** — columns resolved by header name.
2. **Comment skipping** — `#` lines before the real header.
3. **Version-suffix stripping** — `gene_id.split(".")[0]`.
4. **Match by name *or* Ensembl ID** — `--gene NKX2-1` and
   `--gene ENSG00000136352` both work.
5. **Graceful degradation** — a missing directory or missing gene warns on
   stderr and is counted, never crashes the run (with 500+ files, one *will*
   be weird).
6. **`math.log2(tpm + 1)`** — the `+1` **pseudocount** avoids `log2(0) = −∞`
   and compresses the skewed TPM scale so the boxplot is readable.


## The awk limits

`scripts/extract_awk.sh` reproduces the whole extraction in ~15 lines of
Bash+awk, and its output is **byte-identical** to the Python version on the
demo bundle. Use it for a first glance; ship the Python version. Why:

| aspect | Bash+awk | Python |
|---|---|---|
| speed of writing | ✅ one-liner-ish | slower |
| column by header name | ❌ fixed `$7`, `$8` | ✅ `DictReader` |
| malformed lines / missing files | ❌ silent weirdness | ✅ try/except + counters |
| log2 transform | `log(x+1)/log(2)` trick | `math.log2` |
| maintainability over a semester | ❌ | ✅ |

The general lesson: **awk is a scalpel for exploration; Python is the
workshop for pipelines.** If an awk script grows a second feature, rewrite it
in Python before it grows a third.

