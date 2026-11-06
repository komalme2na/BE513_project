
#!/usr/bin/env python3
"""
quick_bed_sort.py
-----------------
Merge one or more BED files belonging to ONE sample, keep only the
chromosomes listed in a selection file, and sort the result in the
*order of the selection file* (NOT alphabetical / natural order).

Usage:
    python3 quick_bed_sort.py --sample X \
        --selection standard_selection.tsv \
        shuf.a.bed.gz shuf.b.bed.gz

Output:
    sorted_bed_file_per_sample/X.standard_selection.sorted.bed
"""
import argparse
import gzip
import os
import sys


def read_selection(path):
    """Return {chromosome: rank} where rank follows FILE order."""
    rank_of = {}
    with open(path) as f:
        for line in f:
            chrom = line.strip().split("\t")[0]
            if chrom and chrom not in rank_of:
                rank_of[chrom] = len(rank_of) + 1
    return rank_of


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True, help="sample name, e.g. X")
    ap.add_argument("--selection", required=True, help="tsv with one chromosome per line")
    ap.add_argument("bed", nargs="+", help="input bed(.gz) files of this sample")
    ap.add_argument("--outdir", default="sorted_bed_file_per_sample")
    args = ap.parse_args()

    rank_of = read_selection(args.selection)
    sel_name = os.path.splitext(os.path.basename(args.selection))[0]
    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, f"{args.sample}.{sel_name}.sorted.bed")

    kept = []                      # (chrom_rank, start, end, original_line)
    for path in args.bed:
        with open_maybe_gz(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                chrom = fields[0]
                if chrom in rank_of:                     # <-- the filter
                    kept.append((rank_of[chrom],          # <-- the ORDER
                                 int(fields[1]), int(fields[2]), line))

    kept.sort(key=lambda t: (t[0], t[1], t[2]))          # rank, start, end

    with open(out_path, "w") as out:
        for *_, original in kept:
            out.write(original + "\n")
    print(f"[quick_bed_sort] kept {len(kept)} lines -> {out_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
