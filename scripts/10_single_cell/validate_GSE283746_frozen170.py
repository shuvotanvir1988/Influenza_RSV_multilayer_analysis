import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

FROZEN = Path(
    "results/DS001_GSE38900/integrated_architecture/tables/"
    "Figure7A_high_confidence_gene_architecture.tsv"
)

DE_DIR = Path(
    "results/scrnaseq_validation/GSE283746/DE"
)

OUTDIR = Path(
    "results/scrnaseq_validation/GSE283746/validation"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Frozen architecture
# ------------------------------------------------------------

f = pd.read_csv(FROZEN, sep="\t")

print("=== FROZEN ARCHITECTURE ===")
print("Rows:", len(f))
print(f["Figure7A_class"].value_counts().to_string())

if len(f) != 170:
    raise RuntimeError(
        f"Expected 170 frozen genes, found {len(f)}"
    )

# ------------------------------------------------------------
# Cell types / filenames
# ------------------------------------------------------------

celltypes = {
    "Naive CD4+ T":
        "GSE283746_Naive_CD4plus_T_RSV_vs_Healthy_edgeR_primary.tsv",

    "Naive CD8+ T":
        "GSE283746_Naive_CD8plus_T_RSV_vs_Healthy_edgeR_primary.tsv",

    "TrB":
        "GSE283746_TrB_RSV_vs_Healthy_edgeR_primary.tsv",

    "B":
        "GSE283746_B_RSV_vs_Healthy_edgeR_primary.tsv",

    "NK":
        "GSE283746_NK_RSV_vs_Healthy_edgeR_primary.tsv",

    "Memory CD4+ T":
        "GSE283746_Memory_CD4plus_T_RSV_vs_Healthy_edgeR_primary.tsv",

    "CD14+ Monocyte":
        "GSE283746_CD14plus_Monocyte_RSV_vs_Healthy_edgeR_primary.tsv",

    "Cytotoxic T":
        "GSE283746_Cytotoxic_T_RSV_vs_Healthy_edgeR_primary.tsv",

    "Treg":
        "GSE283746_Treg_RSV_vs_Healthy_edgeR_primary.tsv",

    "Memory B":
        "GSE283746_Memory_B_RSV_vs_Healthy_edgeR_primary.tsv",

    "CD16+ Monocyte":
        "GSE283746_CD16plus_Monocyte_RSV_vs_Healthy_edgeR_primary.tsv",
}

summary_rows = []
master_rows = []

for ct, fn in celltypes.items():

    d = pd.read_csv(
        DE_DIR / fn,
        sep="\t"
    )

    # Defensive collapse if duplicate symbols ever occur
    d["_rank"] = d["FDR"].fillna(d["PValue"])

    d = (
        d.sort_values(
            ["gene_symbol", "_rank"]
        )
        .drop_duplicates(
            "gene_symbol",
            keep="first"
        )
    )

    z = f.merge(
        d[
            [
                "gene_symbol",
                "logFC",
                "FDR",
                "PValue",
                "logCPM"
            ]
        ],
        on="gene_symbol",
        how="left"
    )

    z["cell_type"] = ct

    z["direction_matches_GPL6884"] = (
        np.sign(z["logFC"]) ==
        np.sign(z["rsv6884_logFC"])
    )

    z["direction_matches_GPL10558"] = (
        np.sign(z["logFC"]) ==
        np.sign(z["rsv10558_logFC"])
    )

    z["direction_matches_both"] = (
        z["direction_matches_GPL6884"] &
        z["direction_matches_GPL10558"]
    )

    master_rows.append(z)

    groups = (
        [("ALL_170", z)] +
        list(
            z.groupby(
                "Figure7A_class",
                sort=False
            )
        )
    )

    for label, sub in groups:

        e = sub.dropna(
            subset=["logFC"]
        ).copy()

        n = len(e)

        match6884 = (
            np.sign(e["logFC"]) ==
            np.sign(e["rsv6884_logFC"])
        )

        match10558 = (
            np.sign(e["logFC"]) ==
            np.sign(e["rsv10558_logFC"])
        )

        both = (
            match6884 & match10558
        )

        if n >= 3:

            p6884 = pearsonr(
                e["rsv6884_logFC"],
                e["logFC"]
            )

            s6884 = spearmanr(
                e["rsv6884_logFC"],
                e["logFC"]
            )

            p10558 = pearsonr(
                e["rsv10558_logFC"],
                e["logFC"]
            )

            s10558 = spearmanr(
                e["rsv10558_logFC"],
                e["logFC"]
            )

        else:

            p6884 = (np.nan, np.nan)
            s6884 = (np.nan, np.nan)
            p10558 = (np.nan, np.nan)
            s10558 = (np.nan, np.nan)

        summary_rows.append({
            "cell_type": ct,
            "class": label,

            "evaluable": n,

            "direction_concordant_GPL6884":
                int(match6884.sum()),

            "direction_concordance_fraction_GPL6884":
                float(match6884.mean())
                if n else np.nan,

            "direction_concordant_GPL10558":
                int(match10558.sum()),

            "direction_concordance_fraction_GPL10558":
                float(match10558.mean())
                if n else np.nan,

            "direction_matches_both":
                int(both.sum()),

            "direction_matches_both_fraction":
                float(both.mean())
                if n else np.nan,

            "pearson_r_GPL6884":
                p6884.statistic
                if hasattr(p6884, "statistic")
                else p6884[0],

            "pearson_p_GPL6884":
                p6884.pvalue
                if hasattr(p6884, "pvalue")
                else p6884[1],

            "spearman_rho_GPL6884":
                s6884.statistic
                if hasattr(s6884, "statistic")
                else s6884[0],

            "spearman_p_GPL6884":
                s6884.pvalue
                if hasattr(s6884, "pvalue")
                else s6884[1],

            "pearson_r_GPL10558":
                p10558.statistic
                if hasattr(p10558, "statistic")
                else p10558[0],

            "pearson_p_GPL10558":
                p10558.pvalue
                if hasattr(p10558, "pvalue")
                else p10558[1],

            "spearman_rho_GPL10558":
                s10558.statistic
                if hasattr(s10558, "statistic")
                else s10558[0],

            "spearman_p_GPL10558":
                s10558.pvalue
                if hasattr(s10558, "pvalue")
                else s10558[1],

            "FDR05":
                int(
                    (e["FDR"] < 0.05)
                    .sum()
                ),

            "FDR05_and_GPL6884_direction":
                int(
                    (
                        (e["FDR"] < 0.05) &
                        match6884
                    ).sum()
                ),

            "FDR05_and_both_microarray_direction":
                int(
                    (
                        (e["FDR"] < 0.05) &
                        both
                    ).sum()
                ),

            "median_scRNA_logFC":
                float(
                    e["logFC"].median()
                )
                if n else np.nan
        })

summary = pd.DataFrame(summary_rows)

master = pd.concat(
    master_rows,
    ignore_index=True
)

summary.to_csv(
    OUTDIR /
    "GSE283746_frozen170_celltype_validation_summary.tsv",
    sep="\t",
    index=False
)

master.to_csv(
    OUTDIR /
    "GSE283746_frozen170_celltype_validation_master.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)

# ------------------------------------------------------------
# Main readable summary
# ------------------------------------------------------------

all170 = summary[
    summary["class"] == "ALL_170"
].copy()

all170 = all170.sort_values(
    "direction_concordance_fraction_GPL6884",
    ascending=False
)

print("\n=== ALL 170 GENES ===")

print(
    all170[
        [
            "cell_type",
            "evaluable",
            "direction_concordant_GPL6884",
            "direction_concordance_fraction_GPL6884",
            "pearson_r_GPL6884",
            "pearson_p_GPL6884",
            "spearman_rho_GPL6884",
            "spearman_p_GPL6884",
            "FDR05",
            "FDR05_and_GPL6884_direction"
        ]
    ].to_string(index=False)
)

print("\n=== SHARED CORE ===")

core = summary[
    summary["class"] == "Shared_core"
].sort_values(
    "direction_concordance_fraction_GPL6884",
    ascending=False
)

print(
    core[
        [
            "cell_type",
            "evaluable",
            "direction_concordance_fraction_GPL6884",
            "pearson_r_GPL6884",
            "spearman_rho_GPL6884",
            "FDR05",
            "FDR05_and_GPL6884_direction"
        ]
    ].to_string(index=False)
)

print("\n=== INFLUENZA-AMPLIFIED SHARED ===")

fluamp = summary[
    summary["class"] ==
    "Influenza_amplified_shared"
].sort_values(
    "direction_concordance_fraction_GPL6884",
    ascending=False
)

print(
    fluamp[
        [
            "cell_type",
            "evaluable",
            "direction_concordance_fraction_GPL6884",
            "pearson_r_GPL6884",
            "spearman_rho_GPL6884",
            "FDR05"
        ]
    ].to_string(index=False)
)

print("\n=== RSV-AMPLIFIED SHARED ===")

rsvamp = summary[
    summary["class"] ==
    "RSV_amplified_shared"
]

print(
    rsvamp[
        [
            "cell_type",
            "evaluable",
            "direction_concordance_fraction_GPL6884",
            "FDR05"
        ]
    ].to_string(index=False)
)

print("\nWritten:")
print(
    OUTDIR /
    "GSE283746_frozen170_celltype_validation_summary.tsv"
)

print(
    OUTDIR /
    "GSE283746_frozen170_celltype_validation_master.tsv.gz"
)
