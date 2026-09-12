import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from pathlib import Path

V = Path(
    "results/rnaseq_validation/validation/"
    "GSE155925_frozen170_validation_master.tsv"
)

E = Path(
    "results/rnaseq_validation/GSE155925_DE/"
    "GSE155925_VST_expression.tsv"
)

M = Path(
    "results/rnaseq_validation/tables/"
    "GSE155925_sample_eligibility_FROZEN_v1.0.tsv"
)

OUT = Path(
    "results/rnaseq_validation/validation/"
    "GSE155925_frozen170_raw_vs_adjusted_audit.tsv"
)

v = pd.read_csv(V, sep="\t")
e = pd.read_csv(E, sep="\t")
m = pd.read_csv(M, sep="\t")

m = m[m["primary_eligible"] == True].copy()

rsv = m.loc[
    m["primary_group"] == "RSV_ONLY",
    "count_matrix_label"
].tolist()

neg = m.loc[
    m["primary_group"] == "VIRUS_NEGATIVE",
    "count_matrix_label"
].tolist()

# Gene symbols from SYMBOL:ENSG
e["gene_symbol"] = (
    e["gene_id"]
    .astype(str)
    .str.split(":", n=1)
    .str[0]
)

records = []

for _, row in v.iterrows():

    symbol = row["gene_symbol"]

    z = e[e["gene_symbol"] == symbol]

    if z.empty:
        continue

    # If duplicate symbol, use first for this diagnostic
    z = z.iloc[0]

    rsv_values = pd.to_numeric(z[rsv])
    neg_values = pd.to_numeric(z[neg])

    raw_difference = (
        rsv_values.mean() - neg_values.mean()
    )

    median_difference = (
        rsv_values.median() - neg_values.median()
    )

    try:
        _, mw_p = mannwhitneyu(
            rsv_values,
            neg_values,
            alternative="two-sided"
        )
    except Exception:
        mw_p = np.nan

    records.append({
        "gene_symbol": symbol,
        "Figure7A_class": row["Figure7A_class"],

        "microarray_RSV_GPL6884_logFC":
            row["rsv6884_logFC"],

        "microarray_RSV_GPL10558_logFC":
            row["rsv10558_logFC"],

        "RNAseq_adjusted_log2FC":
            row["rnaseq_RSV_log2FC"],

        "RNAseq_unadjusted_VST_mean_difference":
            raw_difference,

        "RNAseq_unadjusted_VST_median_difference":
            median_difference,

        "MannWhitney_p":
            mw_p,

        "raw_direction_matches_microarray":
            np.sign(raw_difference) ==
            np.sign(row["rsv6884_logFC"]),

        "adjusted_direction_matches_microarray":
            np.sign(row["rnaseq_RSV_log2FC"]) ==
            np.sign(row["rsv6884_logFC"]),
    })

d = pd.DataFrame(records)

d.to_csv(OUT, sep="\t", index=False)

print("=== RAW VS ADJUSTED DIRECTION AUDIT ===")

print("\nGenes:", len(d))

print("\nRaw VST direction concordance:")
print(
    d["raw_direction_matches_microarray"]
    .value_counts()
)

print("\nAdjusted DESeq2 direction concordance:")
print(
    d["adjusted_direction_matches_microarray"]
    .value_counts()
)

print("\n=== BY FROZEN CLASS ===")

print(
    d.groupby("Figure7A_class")[
        [
            "raw_direction_matches_microarray",
            "adjusted_direction_matches_microarray",
        ]
    ].mean()
)

targets = [
    "IFI27","IFI44","IFI44L",
    "IFIT1","IFIT2","IFIT3",
    "ISG15","MX1","MX2",
    "RSAD2","IRF7","OAS3",
    "USP18","HERC5"
]

print("\n=== INTERFERON TARGET AUDIT ===")

print(
    d[d["gene_symbol"].isin(targets)][
        [
            "gene_symbol",
            "microarray_RSV_GPL6884_logFC",
            "RNAseq_unadjusted_VST_mean_difference",
            "RNAseq_adjusted_log2FC",
            "raw_direction_matches_microarray",
            "adjusted_direction_matches_microarray",
            "MannWhitney_p",
        ]
    ].to_string(index=False)
)

print("\nWritten:", OUT)
