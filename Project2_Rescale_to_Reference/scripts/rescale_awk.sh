#!/usr/bin/env bash
# rescale_awk.sh — pure Bash+awk implementation of the same algorithm as
# rescale.py.  Useful when Python is unavailable and as a study example of
# two-pass streaming with awk associative arrays.
#
# Usage:
#   bash scripts/rescale_awk.sh data/query.bed.gz data/reference.hist results/rescaled_awk.bed
#
# NOTE: awk's rand() is a simple LCG — fine for teaching, but prefer Python's
# random/numpy RNGs (better period, seedable across platforms) in real work.

set -euo pipefail

QUERY=${1:?usage: rescale_awk.sh QUERY.bed[.gz] REFERENCE.hist OUT.bed}
REF=${2:?usage: rescale_awk.sh QUERY.bed[.gz] REFERENCE.hist OUT.bed}
OUT=${3:?usage: rescale_awk.sh QUERY.bed[.gz] REFERENCE.hist OUT.bed}

# zcat for .gz, cat otherwise
CAT="cat"; [[ "$QUERY" == *.gz ]] && CAT="zcat"

# -------- pass 1: build the query length histogram --------
QHIST=$(mktemp)
trap 'rm -f "$QHIST"' EXIT
$CAT "$QUERY" | awk -F'\t' '!/^#/ && NF>=3 { c[$3-$2]++ ; n++ }
                 END { for (L in c) printf "%d\t%.10f\n", L, c[L]/n }' > "$QHIST"

# -------- pass 2: compute p(L) on the fly and accept/reject --------
$CAT "$QUERY" | awk -F'\t' -v OFS='\t' -v qhist="$QHIST" -v ref="$REF" '
BEGIN { srand(1)
        while ((getline line < qhist) > 0) { split(line,a,"\t"); Q[a[1]]=a[2] }
        close(qhist)
        RSUM=0
        while ((getline line < ref) > 0) {
            if (line ~ /^#/ || line == "") continue
            split(line,a,/[ \t]+/); REF[a[1]]=a[2]; RSUM+=a[2]
        }
        close(ref)
        for (L in REF) REF[L]/=RSUM
        MAXR=0
        for (L in REF) if (L in Q && Q[L]>0) { r=REF[L]/Q[L]; if (r>MAXR) MAXR=r }
        for (L in REF) if (L in Q && Q[L]>0) P[L]=(REF[L]/Q[L])/MAXR
      }
!/^#/ && NF>=3 { len=$3-$2
                 if (len in P && rand() < P[len]) print }' > "$OUT"

echo "wrote $OUT ($(wc -l < "$OUT") lines)" >&2
