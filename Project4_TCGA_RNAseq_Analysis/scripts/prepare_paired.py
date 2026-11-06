
#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

input_file = "results/NKX2-1_expression.tsv"
output_file = "results/NKX2-1_paired_expression.tsv"

patients = defaultdict(list)

with open(input_file) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        patients[row["case_id"]].append(row)

paired = []

for case_id, rows in patients.items():

    tumour = [
        r for r in rows
        if r["sample_type"] == "Primary Tumor"
    ]

    normal = [
        r for r in rows
        if r["sample_type"] == "Solid Tissue Normal"
    ]

    if len(tumour) == 1 and len(normal) == 1:

        paired.append({
            "case_id": case_id,
            "tumour_tpm": tumour[0]["tpm_unstranded"],
            "normal_tpm": normal[0]["tpm_unstranded"],
            "tumour_log2": tumour[0]["log2_tpm_plus1"],
            "normal_log2": normal[0]["log2_tpm_plus1"]
        })

Path(output_file).parent.mkdir(parents=True, exist_ok=True)

with open(output_file, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "case_id",
            "tumour_tpm",
            "normal_tpm",
            "tumour_log2",
            "normal_log2"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(paired)

print("Matched patients:", len(paired))
print("Output:", output_file)
