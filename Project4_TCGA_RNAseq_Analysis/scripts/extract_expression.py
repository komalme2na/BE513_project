
#!/usr/bin/env python3
"""
Extract NKX2-1 TPM expression from GDC STAR-counts files.

Input:
    GDC sample sheet
    Extracted TCGA STAR-counts directories

Output:
    Tidy TSV containing sample ID, patient ID,
    sample type, TPM and log2(TPM + 1).
"""

import argparse
import csv
import math
import sys
from pathlib import Path

KEEP_TYPES = {"Primary Tumor", "Solid Tissue Normal"}


def read_sample_sheet(path):
    """Map File ID to Case ID and Sample Type."""

    mapping = {}

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {"File ID", "Case ID", "Sample Type"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            sys.exit(f"Missing columns: {missing}")

        for row in reader:
            sample_type = row["Sample Type"].strip()

            if sample_type not in KEEP_TYPES:
                continue

            filename = row.get("File Name", "")

            # Ignore files that are not STAR gene counts
            if not filename.endswith(
                ".rna_seq.augmented_star_gene_counts.tsv"
            ):
                continue

            file_id = row["File ID"].strip()
            case_id = row["Case ID"].strip()

            mapping[file_id] = (case_id, sample_type)

    return mapping


def find_gene_tpm(counts_file, gene):
    """Extract tpm_unstranded for the requested gene."""

    with open(counts_file) as fh:

        header = None

        for line in fh:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if header is None:
                header = fields

                try:
                    tpm_idx = header.index("tpm_unstranded")
                except ValueError:
                    raise ValueError(
                        f"tpm_unstranded column not found: {counts_file}"
                    )

                continue

            gene_id = fields[0].split(".")[0]
            gene_name = fields[1] if len(fields) > 1 else ""

            if gene_name == gene or gene_id == gene.split(".")[0]:
                return float(fields[tpm_idx])

    raise KeyError(f"{gene} not found in {counts_file}")


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--sheet", required=True)
    parser.add_argument("--datadir", required=True)
    parser.add_argument("--gene", default="NKX2-1")
    parser.add_argument("--out", required=True)

    args = parser.parse_args()

    mapping = read_sample_sheet(args.sheet)

    print(
        f"[info] PT/STN STAR-counts files: {len(mapping)}",
        file=sys.stderr
    )

    datadir = Path(args.datadir)
    output = Path(args.out)

    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    missing = 0

    for file_id, (case_id, sample_type) in sorted(mapping.items()):

        matches = list(
            datadir.glob(
                f"{file_id}/*.rna_seq.augmented_star_gene_counts.tsv"
            )
        )

        if not matches:
            missing += 1
            continue

        try:
            tpm = find_gene_tpm(matches[0], args.gene)

        except (KeyError, ValueError) as error:
            print(f"[warn] {error}", file=sys.stderr)
            missing += 1
            continue

        rows.append([
            file_id,
            case_id,
            sample_type,
            tpm,
            math.log2(tpm + 1)
        ])

    with open(output, "w", newline="") as fh:

        writer = csv.writer(fh, delimiter="\t")

        writer.writerow([
            "file_id",
            "case_id",
            "sample_type",
            "tpm_unstranded",
            "log2_tpm_plus1"
        ])

        writer.writerows(rows)

    tumour = sum(r[2] == "Primary Tumor" for r in rows)
    normal = sum(r[2] == "Solid Tissue Normal" for r in rows)

    print(f"[info] Output: {output}")
    print(f"[info] Total samples: {len(rows)}")
    print(f"[info] Primary Tumor: {tumour}")
    print(f"[info] Solid Tissue Normal: {normal}")
    print(f"[info] Missing files/genes: {missing}")


if __name__ == "__main__":
    main()
