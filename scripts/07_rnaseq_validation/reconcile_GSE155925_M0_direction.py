import pandas as pd
import numpy as np
from pathlib import Path

FROZEN = Path(
    "results/DS001_GSE38900/integrated_architecture/tables/"
    "Figure7A_high_confidence_gene_architecture.tsv"
)

M0 = Path(
    "results/rnaseq_validation/GSE155925_DE/model_sensitivity/"
    "GSE155925_M0_unadjusted_RSV_vs_negative.tsv"
)

NORM = Path(
    "results/rnaseq_validation/GSE155925_DE/"
    "GSE155925_normalized_counts.tsv"
)

META = Path(
    "results/rnaseq_validation/tables/"
    "GSE155925_sample_eligibility_FROZEN_v1.0.tsv"
)

OUT = Path(
    "results/rnaseq_validation/validation/"
    "GSE155925_M0_direction_reconciliation.tsv"
)

frozen = pd.read_csv(FROZEN, sep="\t")
m0 = pd.read_csv(M0, sep="\t")
norm = pd.read_csv(NORM, sep="\t")
meta = pd.read_csv(META, sep="\t")

meta = meta[meta["primary_eligible"] == True]

rsv = meta.loc[
    meta["primary_group"] == "RSV_ONLY",
    "count_matrix_label"
].tolist()

neg = meta.loc[
    meta["primary_group"] == "VIRUS_NEGATIVE",
    "count_matrix_label"
].tolist()

# ----------------------------------------------------------
# Feature parsing
# ----------------------------------------------------------

for d in [m0, norm]:
    d["gene_symbol"] = (
        d["gene_id"]
        .astype(str)
        .str.split(":", n=1)
        .str[0]
    )

# Use EXACT M0 feature chosen by smallest padj/pvalue
m0["_rank"] = m0["padj"].fillna(m0["pvalue"])

m0_unique = (
    m0.sort_values(["gene_symbol", "_rank"])
      .drop_duplicates("gene_symbol", keep="first")
      .copy()
)

# Join normalized counts using exact gene_id chosen in M0
z = frozen.merge(
    m0_unique[
        [
            "gene_symbol",
            "gene_id",
            "log2FoldChange",
            "padj"
        ]
    ],
    on="gene_symbol",
    how="left"
)

norm_exact = norm.drop(columns=["gene_symbol"])

z = z.merge(
    norm_exact,
    on="gene_id",
    how="left"
)

records = []

for _, row in z.iterrows():

    if pd.isna(row.get("log2FoldChange")):
        continue

    rsv_vals = pd.to_numeric(row[rsv], errors="coerce")
    neg_vals = pd.to_numeric(row[neg], errors="coerce")

    rsv_mean = rsv_vals.mean()
    neg_mean = neg_vals.mean()

    # log2 ratio of mean normalized counts
    raw_logratio = np.log2(
        (rsv_mean + 0.5) /
        (neg_mean + 0.5)
    )

    micro_sign = np.sign(row["rsv6884_logFC"])
    m0_sign = np.sign(row["log2FoldChange"])
    mean_sign = np.sign(raw_logratio)

    records.append({
        "gene_symbol": row["gene_symbol"],
        "Figure7A_class": row["Figure7A_class"],
        "gene_id": row["gene_id"],

        "microarray_RSV_logFC":
            row["rsv6884_logFC"],

        "M0_DESeq2_log2FC":
            row["log2FoldChange"],

        "normalized_mean_log2ratio":
            raw_logratio,

        "M0_matches_microarray":
            m0_sign == micro_sign,

        "mean_ratio_matches_microarray":
            mean_sign == micro_sign,

        "M0_matches_mean_ratio":
            m0_sign == mean_sign,

        "M0_padj":
            row["padj"],
    })

out = pd.DataFrame(records)
out.to_csv(OUT, sep="\t", index=False)

print("=== M0 DIRECTION RECONCILIATION ===")
print("Evaluable:", len(out))

print("\nM0 coefficient vs normalized-count mean ratio:")
print(out["M0_matches_mean_ratio"].value_counts())

print("\nM0 vs frozen microarray:")
print(out["M0_matches_microarray"].value_counts())

print("\nNormalized mean ratio vs frozen microarray:")
print(out["mean_ratio_matches_microarray"].value_counts())

print("\n=== BY CLASS ===")
print(
    out.groupby("Figure7A_class")[
        [
            "M0_matches_microarray",
            "mean_ratio_matches_microarray",
            "M0_matches_mean_ratio"
        ]
    ].mean()
)

print("\n=== DISAGREEMENT BETWEEN M0 AND MEAN RATIO ===")
bad = out[out["M0_matches_mean_ratio"] == False]

print(
    bad[
        [
            "gene_symbol",
            "Figure7A_class",
            "M0_DESeq2_log2FC",
            "normalized_mean_log2ratio",
            "M0_padj"
        ]
    ].to_string(index=False)
)

print("\nWritten:", OUT)
