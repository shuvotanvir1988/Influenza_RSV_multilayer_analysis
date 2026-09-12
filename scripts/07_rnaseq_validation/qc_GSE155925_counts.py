import pandas as pd
import numpy as np
from pathlib import Path

COUNTFILE = Path(
    "data/rnaseq_validation/GSE155925/raw_counts/"
    "GSE155925_Raw_counts_matrix.txt.gz"
)

ELIGFILE = Path(
    "results/rnaseq_validation/tables/"
    "GSE155925_sample_eligibility_FROZEN_v1.0.tsv"
)

OUTDIR = Path("results/rnaseq_validation/qc/GSE155925")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Load counts
# ----------------------------
counts = pd.read_csv(COUNTFILE, sep="\t")

counts = counts.rename(columns={counts.columns[0]: "gene_id"})

print("=== FULL COUNT MATRIX ===")
print("Genes:", counts.shape[0])
print("Samples:", counts.shape[1] - 1)

# Validate uniqueness
print("Unique gene IDs:", counts["gene_id"].nunique())
print("Duplicated gene IDs:", counts["gene_id"].duplicated().sum())

# ----------------------------
# Load frozen eligibility
# ----------------------------
meta = pd.read_csv(ELIGFILE, sep="\t")

eligible = meta.loc[
    meta["primary_eligible"] == True
].copy()

eligible_labels = eligible["count_matrix_label"].tolist()

missing = sorted(set(eligible_labels) - set(counts.columns))

if missing:
    raise RuntimeError(
        f"Eligible samples missing from count matrix: {missing}"
    )

sub = counts[["gene_id"] + eligible_labels].copy()

# ----------------------------
# Library-level QC
# ----------------------------
mat = sub.set_index("gene_id")

qc = pd.DataFrame(index=eligible_labels)

qc["library_size"] = mat.sum(axis=0)
qc["detected_genes_count_gt0"] = (mat > 0).sum(axis=0)
qc["genes_count_ge10"] = (mat >= 10).sum(axis=0)
qc["zero_fraction"] = (mat == 0).mean(axis=0)

qc.index.name = "count_matrix_label"
qc = qc.reset_index()

qc = qc.merge(
    eligible[
        [
            "count_matrix_label",
            "geo_accession",
            "primary_group",
            "age_months",
            "sex",
            "hospital_batch",
            "enrollment_year_batch",
        ]
    ],
    on="count_matrix_label",
    how="left",
    validate="one_to_one"
)

qc.to_csv(
    OUTDIR / "GSE155925_sample_level_QC.tsv",
    sep="\t",
    index=False
)

# ----------------------------
# Gene-level count filtering audit
# ----------------------------

n = mat.shape[1]

filter_summary = []

for min_count, min_samples in [
    (1, 1),
    (10, 1),
    (10, 2),
    (10, 5),
    (10, 10),
]:
    keep = (mat >= min_count).sum(axis=1) >= min_samples

    filter_summary.append({
        "minimum_count": min_count,
        "minimum_samples": min_samples,
        "genes_retained": int(keep.sum()),
        "genes_removed": int((~keep).sum()),
        "fraction_retained": float(keep.mean()),
    })

filter_df = pd.DataFrame(filter_summary)

filter_df.to_csv(
    OUTDIR / "GSE155925_count_filtering_audit.tsv",
    sep="\t",
    index=False
)

# ----------------------------
# Summary
# ----------------------------

print("\n=== ELIGIBLE MATRIX ===")
print("Genes:", mat.shape[0])
print("Samples:", mat.shape[1])

print("\n=== LIBRARY SIZE ===")
print(qc["library_size"].describe())

print("\n=== DETECTED GENES (>0) ===")
print(qc["detected_genes_count_gt0"].describe())

print("\n=== GENES WITH COUNT >=10 ===")
print(qc["genes_count_ge10"].describe())

print("\n=== ZERO FRACTION ===")
print(qc["zero_fraction"].describe())

print("\n=== FILTERING AUDIT ===")
print(filter_df.to_string(index=False))

print("\n=== SAMPLE QC TABLE ===")
print(
    qc.sort_values("library_size")[
        [
            "count_matrix_label",
            "geo_accession",
            "primary_group",
            "library_size",
            "detected_genes_count_gt0",
            "genes_count_ge10",
            "zero_fraction",
        ]
    ].to_string(index=False)
)
