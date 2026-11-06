
# BED format deep dive

## The coordinate system (the #1 source of bugs in genomics)
BED is **0-based, half-open**: the interval `[start, end)` contains bases
`start, start+1, …, end−1`. Length = `end − start`.
GFF/GTF and 1-based formats use `[start, end]` — length = `end − start + 1`.
Converting: `bed_start = gff_start − 1`, `bed_end = gff_end`.
Whenever two files disagree about a coordinate convention by one base, an
entire analysis can shift — always check before `intersect`ing.

## BED3 … BED12
| cols | name | notes |
|---|---|---|
| 1–3 | chrom, start, end | required |
| 4 | name | any string; in cfDNA fragment files often encodes length/mapping info |
| 5 | score | 0–1000 convention |
| 6 | strand | `+` / `-` / `.` |
| 7–12 | thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts | gene-model display |

## Common companion formats
- **bedGraph** (4 cols: fixedStep-style signal) — same coordinate convention.
- **narrowPeak/broadPeak** (ENCODE peaks) — BED6 + 4 extra statistical columns.
- **bigBed / bigWig** — binary indexed versions; tools: `bedToBigBed`, `wigToBigWig`.

## Validation one-liners
```bash
awk -F'\t' '$2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2<0 || $3<=$2 {bad++} END{print bad+0" bad lines"}' file.bed
zcat file.bed.gz | awk -F'\t' 'NF<3{bad++} END{print bad+0" short lines"}'
```
