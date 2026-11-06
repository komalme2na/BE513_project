
#!/usr/bin/env bash

# Usage:
# ./merge_sort_bed.sh <sample> <selection.tsv> <outdir> <bed1.gz> [bed2.gz ...]

set -euo pipefail

# 1. Check arguments
if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <sample> <selection.tsv> <outdir> <bed1.gz> [bed2.gz ...]" >&2
    exit 1
fi

SAMPLE=$1
SELECTION=$2
OUTDIR=$3
shift 3

# 2. Check input files
if [[ ! -f "$SELECTION" || ! -r "$SELECTION" ]]; then
    echo "Error: Selection file not found or unreadable." >&2
    exit 1
fi

for file in "$@"; do
    if [[ ! -f "$file" || ! -r "$file" ]]; then
        echo "Error: BED file not found or unreadable: $file" >&2
        exit 1
    fi
    gzip -t -- "$file" || {
        echo "Error: Invalid gzip file: $file" >&2
        exit 1
    }
done

# 3. Create output directory
mkdir -p "$OUTDIR"

# 4. Generate output filename
SEL_NAME=$(basename "$SELECTION" .tsv)
OUT="$OUTDIR/${SAMPLE}.${SEL_NAME}.sorted.bed"

# 5. Create temporary output
TMP=$(mktemp "$OUTDIR/.${SAMPLE}.${SEL_NAME}.XXXXXX")
trap 'rm -f "$TMP"' EXIT

# 6. Filter, rank and sort BED records
zcat -- "$@" |
awk -F '\t' -v OFS='\t' -v selection_file="$SELECTION" '
    FILENAME == selection_file {
        sub(/\r$/, "", $1)

        if ($1 != "" && !($1 in rank)) {
            rank[$1] = ++count
        }
        next
    }

    NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2 {
        print "Error: Invalid BED record at input line " FNR > "/dev/stderr"
        exit 1
    }

    $1 in rank {
        print rank[$1], $2, $3, $0
    }

    END {
        if (count == 0) {
            print "Error: Selection file is empty." > "/dev/stderr"
            exit 1
        }
    }
' "$SELECTION" - |
LC_ALL=C sort -t $'\t' -k1,1n -k2,2n -k3,3n |
cut -f4- > "$TMP"

# 7. Save final output
mv "$TMP" "$OUT"

echo "Output saved: $OUT"

