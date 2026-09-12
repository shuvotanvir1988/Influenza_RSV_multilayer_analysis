import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

FROZEN = Path(
    "results/DS001_GSE38900/integrated_architecture/tables/"
    "Figure7A_high_confidence_gene_architecture.tsv"
)

DE_DIR = Path(
    "results/scrnaseq_validation/GSE283746/isg_state/DE"
)

OUTDIR = Path(
    "results/scrnaseq_validation/GSE283746/isg_state/validation"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

states = {
    "ISG_Naive_CD4_T":
        DE_DIR /
        "GSE283746_ISG_Naive_CD4_T_RSV_vs_Healthy_edgeR.tsv",

    "Reference_Naive_CD4_T":
        DE_DIR /
        "GSE283746_Reference_Naive_CD4_T_RSV_vs_Healthy_edgeR.tsv"
}

f = pd.read_csv(FROZEN, sep="\t")

print("=== FROZEN ARCHITECTURE ===")
print("Rows:", len(f))
print(f["Figure7A_class"].value_counts().to_string())

if len(f) != 170:
    raise RuntimeError(
        f"Expected 170 frozen genes, found {len(f)}"
    )

summary_rows = []
master_rows = []

for state_id, de_file in states.items():

    print("\n========================================")
    print(state_id)
    print("========================================")

    if not de_file.exists():
        raise FileNotFoundError(de_file)

    d = pd.read_csv(
        de_file,
        sep="\t"
    )

    d = (
        d.sort_values(
            ["gene_symbol", "FDR", "PValue"]
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

    z["state_id"] = state_id

    z["direction_matches_microarray"] = (
        np.sign(z["logFC"]) ==
        np.sign(z["rsv6884_logFC"])
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

        if n:
            match = (
                np.sign(e["logFC"]) ==
                np.sign(e["rsv6884_logFC"])
            )
        else:
            match = pd.Series(
                dtype=bool
            )

        if n >= 3:

            pr = pearsonr(
                e["rsv6884_logFC"],
                e["logFC"]
            )

            sr = spearmanr(
                e["rsv6884_logFC"],
                e["logFC"]
            )

            pearson_r = pr.statistic
            pearson_p = pr.pvalue
            spearman_rho = sr.statistic
            spearman_p = sr.pvalue

        else:

            pearson_r = np.nan
            pearson_p = np.nan
            spearman_rho = np.nan
            spearman_p = np.nan

        summary_rows.append({
            "state_id": state_id,
            "class": label,
            "evaluable": n,

            "direction_concordant":
                int(match.sum())
                if n else 0,

            "direction_concordance_fraction":
                float(match.mean())
                if n else np.nan,

            "pearson_r":
                pearson_r,

            "pearson_p":
                pearson_p,

            "spearman_rho":
                spearman_rho,

            "spearman_p":
                spearman_p,

            "FDR05":
                int(
                    (e["FDR"] < 0.05)
                    .sum()
                ),

            "FDR05_direction_concordant":
                int(
                    (
                        (e["FDR"] < 0.05) &
                        match
                    ).sum()
                ),

            "median_logFC":
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
    "GSE283746_NaiveCD4_ISG_frozen170_validation_summary.tsv",
    sep="\t",
    index=False
)

master.to_csv(
    OUTDIR /
    "GSE283746_NaiveCD4_ISG_frozen170_validation_master.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)

print("\n=== FINAL COMPARISON ===")

show_classes = [
    "ALL_170",
    "Shared_core",
    "Influenza_amplified_shared",
    "RSV_amplified_shared"
]

for cls in show_classes:

    x = summary[
        summary["class"] == cls
    ].copy()

    print("\n----------------------------------------")
    print(cls)
    print("----------------------------------------")

    print(
        x[
            [
                "state_id",
                "evaluable",
                "direction_concordant",
                "direction_concordance_fraction",
                "pearson_r",
                "spearman_rho",
                "FDR05",
                "FDR05_direction_concordant",
                "median_logFC"
            ]
        ].to_string(index=False)
    )

# Direct ISG-minus-reference summary
wide = (
    summary[
        summary["class"].isin(
            [
                "ALL_170",
                "Shared_core",
                "Influenza_amplified_shared"
            ]
        )
    ]
    .pivot(
        index="class",
        columns="state_id",
        values=[
            "direction_concordance_fraction",
            "pearson_r",
            "spearman_rho",
            "FDR05_direction_concordant"
        ]
    )
)

print("\n=== SIDE-BY-SIDE METRICS ===")
print(wide.to_string())

print("\nWritten:")
print(
    OUTDIR /
    "GSE283746_NaiveCD4_ISG_frozen170_validation_summary.tsv"
)

print(
    OUTDIR /
    "GSE283746_NaiveCD4_ISG_frozen170_validation_master.tsv.gz"
)
