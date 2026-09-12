#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path("results/DS001_GSE38900/signature_analysis")
TABLES = ROOT / "tables"
SENS = ROOT / "sensitivity"

TABLES.mkdir(parents=True, exist_ok=True)
SENS.mkdir(parents=True, exist_ok=True)

MASTER = TABLES / "DS001_GSE38900_gene_evidence_master.tsv"

x = pd.read_csv(MASTER, sep="\t")


def corr_stats(df, xcol, ycol, label):
    z = df[[xcol, ycol]].dropna()

    if len(z) < 3:
        return {
            "comparison": label,
            "n": len(z),
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
        }

    pr = pearsonr(z[xcol], z[ycol])
    sr = spearmanr(z[xcol], z[ycol])

    return {
        "comparison": label,
        "n": len(z),
        "pearson_r": pr.statistic,
        "pearson_p": pr.pvalue,
        "spearman_rho": sr.statistic,
        "spearman_p": sr.pvalue,
    }


rows = []

# ---------------------------------------------------------
# 1. Genome-wide influenza vs RSV effect-size concordance
# ---------------------------------------------------------

rows.append(
    corr_stats(
        x,
        "flu_logFC",
        "rsv6884_logFC",
        "Influenza_vs_control_vs_RSV_GPL6884_vs_control_all_genes"
    )
)

# Shared significant genes on GPL6884
rows.append(
    corr_stats(
        x.loc[x["shared_GPL6884"]],
        "flu_logFC",
        "rsv6884_logFC",
        "Influenza_vs_RSV_shared_GPL6884"
    )
)

# High-confidence replicated shared genes
rows.append(
    corr_stats(
        x.loc[x["shared_replicated"]],
        "flu_logFC",
        "rsv6884_logFC",
        "Influenza_vs_RSV_shared_replicated"
    )
)

# ---------------------------------------------------------
# 2. RSV cross-platform effect-size concordance
# ---------------------------------------------------------

rows.append(
    corr_stats(
        x.loc[x["rsv_crossplatform_evaluable"]],
        "rsv10558_logFC",
        "rsv6884_logFC",
        "RSV_GPL10558_vs_GPL6884_all_evaluable"
    )
)

rows.append(
    corr_stats(
        x.loc[x["rsv_replicated"]],
        "rsv10558_logFC",
        "rsv6884_logFC",
        "RSV_GPL10558_vs_GPL6884_replicated"
    )
)

rows.append(
    corr_stats(
        x.loc[x["shared_replicated"]],
        "rsv10558_logFC",
        "rsv6884_logFC",
        "RSV_GPL10558_vs_GPL6884_shared_replicated"
    )
)

corr = pd.DataFrame(rows)

corr.to_csv(
    TABLES / "DS001_GSE38900_signature_effect_size_correlations.tsv",
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# 3. Directional concordance
# ---------------------------------------------------------

def direction_summary(df, label, col1, col2):
    z = df[[col1, col2]].dropna().copy()

    s1 = np.sign(z[col1])
    s2 = np.sign(z[col2])

    same = (s1 == s2)
    opposite = (s1 == -s2)

    return {
        "comparison": label,
        "n": len(z),
        "same_direction_n": int(same.sum()),
        "same_direction_fraction": float(same.mean()),
        "opposite_direction_n": int(opposite.sum()),
        "opposite_direction_fraction": float(opposite.mean()),
    }


direction_rows = []

direction_rows.append(
    direction_summary(
        x,
        "Influenza_vs_RSV_GPL6884_all_genes",
        "flu_logFC",
        "rsv6884_logFC"
    )
)

direction_rows.append(
    direction_summary(
        x.loc[x["rsv_crossplatform_evaluable"]],
        "RSV_crossplatform_all_evaluable",
        "rsv10558_logFC",
        "rsv6884_logFC"
    )
)

direction_rows.append(
    direction_summary(
        x.loc[x["rsv_replicated"]],
        "RSV_crossplatform_replicated",
        "rsv10558_logFC",
        "rsv6884_logFC"
    )
)

direction = pd.DataFrame(direction_rows)

direction.to_csv(
    TABLES / "DS001_GSE38900_signature_direction_concordance.tsv",
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# 4. Algebraic audit of direct contrast
#
# Expected:
# flu_vs_rsv_logFC =
# flu_logFC - rsv6884_logFC
# ---------------------------------------------------------

x["expected_direct_logFC"] = (
    x["flu_logFC"] - x["rsv6884_logFC"]
)

x["direct_logFC_difference"] = (
    x["flu_vs_rsv_logFC"] - x["expected_direct_logFC"]
)

x["direct_logFC_abs_difference"] = \
    x["direct_logFC_difference"].abs()

audit = pd.DataFrame({
    "metric": [
        "n_genes",
        "maximum_absolute_difference",
        "mean_absolute_difference",
        "median_absolute_difference",
        "n_abs_difference_gt_1e-10",
        "n_abs_difference_gt_1e-8",
        "n_abs_difference_gt_1e-6",
    ],
    "value": [
        len(x),
        x["direct_logFC_abs_difference"].max(),
        x["direct_logFC_abs_difference"].mean(),
        x["direct_logFC_abs_difference"].median(),
        int((x["direct_logFC_abs_difference"] > 1e-10).sum()),
        int((x["direct_logFC_abs_difference"] > 1e-8).sum()),
        int((x["direct_logFC_abs_difference"] > 1e-6).sum()),
    ]
})

audit.to_csv(
    SENS / "DS001_GSE38900_direct_contrast_algebra_audit.tsv",
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# 5. Shared-signature robustness counts
# ---------------------------------------------------------

shared = x.loc[x["shared_GPL6884"]].copy()

robustness_rows = [
    ("shared_GPL6884_total", len(shared)),
    (
        "shared_GPL6884_crossplatform_evaluable",
        int(shared["rsv_crossplatform_evaluable"].sum())
    ),
    (
        "shared_GPL6884_not_crossplatform_evaluable",
        int((~shared["rsv_crossplatform_evaluable"]).sum())
    ),
    (
        "shared_GPL6884_RSV_replicated",
        int(shared["rsv_replicated"].sum())
    ),
    (
        "shared_GPL6884_RSV_not_replicated",
        int(
            (
                shared["rsv_crossplatform_evaluable"]
                &
                ~shared["rsv_replicated"]
            ).sum()
        )
    ),
    (
        "shared_replicated_direct_not_significant",
        int(
            (
                shared["rsv_replicated"]
                &
                ~shared["direct_significant"]
            ).sum()
        )
    ),
    (
        "shared_replicated_direct_significant",
        int(
            (
                shared["rsv_replicated"]
                &
                shared["direct_significant"]
            ).sum()
        )
    ),
    (
        "shared_replicated_influenza_biased",
        int(
            (
                shared["rsv_replicated"]
                &
                shared["direct_significant"]
                &
                (shared["flu_vs_rsv_logFC"] > 0)
            ).sum()
        )
    ),
    (
        "shared_replicated_RSV_biased",
        int(
            (
                shared["rsv_replicated"]
                &
                shared["direct_significant"]
                &
                (shared["flu_vs_rsv_logFC"] < 0)
            ).sum()
        )
    ),
]

robustness = pd.DataFrame(
    robustness_rows,
    columns=["metric", "value"]
)

robustness.to_csv(
    TABLES / "DS001_GSE38900_shared_signature_robustness_summary.tsv",
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# 6. Biological/pathway-coherence counts
# ---------------------------------------------------------

groups = {
    "shared_GPL6884":
        x["shared_GPL6884"],

    "shared_replicated":
        x["shared_replicated"],

    "shared_replicated_no_pathogen_bias":
        (
            x["shared_replicated"]
            &
            ~x["direct_significant"]
        ),

    "shared_replicated_influenza_biased":
        (
            x["shared_replicated"]
            &
            x["direct_significant"]
            &
            (x["flu_vs_rsv_logFC"] > 0)
        ),

    "shared_replicated_RSV_biased":
        (
            x["shared_replicated"]
            &
            x["direct_significant"]
            &
            (x["flu_vs_rsv_logFC"] < 0)
        ),

    "influenza_biased_direct_all":
        x["influenza_biased_direct"],

    "RSV_biased_direct_all":
        x["rsv_biased_direct"],
}

coherence_rows = []

for label, mask in groups.items():
    z = x.loc[mask]

    n = len(z)
    n_le = int(z["session17_leading_edge"].sum())

    coherence_rows.append({
        "group": label,
        "n_genes": n,
        "leading_edge_genes": n_le,
        "leading_edge_fraction":
            n_le / n if n else np.nan
    })

coherence = pd.DataFrame(coherence_rows)

coherence.to_csv(
    TABLES / "DS001_GSE38900_signature_pathway_coherence.tsv",
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# 7. Export principal gene sets
# ---------------------------------------------------------

exports = {
    "DS001_GSE38900_shared_GPL6884_genes.tsv":
        x.loc[x["shared_GPL6884"]],

    "DS001_GSE38900_shared_replicated_genes.tsv":
        x.loc[x["shared_replicated"]],

    "DS001_GSE38900_shared_replicated_core.tsv":
        x.loc[
            x["shared_replicated"]
            &
            ~x["direct_significant"]
        ],

    "DS001_GSE38900_shared_replicated_influenza_biased.tsv":
        x.loc[
            x["shared_replicated"]
            &
            x["direct_significant"]
            &
            (x["flu_vs_rsv_logFC"] > 0)
        ],

    "DS001_GSE38900_shared_replicated_RSV_biased.tsv":
        x.loc[
            x["shared_replicated"]
            &
            x["direct_significant"]
            &
            (x["flu_vs_rsv_logFC"] < 0)
        ],

    "DS001_GSE38900_influenza_biased_direct.tsv":
        x.loc[x["influenza_biased_direct"]],

    "DS001_GSE38900_RSV_biased_direct.tsv":
        x.loc[x["rsv_biased_direct"]],
}

for filename, df in exports.items():
    df.to_csv(TABLES / filename, sep="\t", index=False)

print("\nEFFECT-SIZE CORRELATIONS")
print(corr.to_string(index=False))

print("\nDIRECTION CONCORDANCE")
print(direction.to_string(index=False))

print("\nDIRECT-CONTRAST ALGEBRA AUDIT")
print(audit.to_string(index=False))

print("\nSHARED SIGNATURE ROBUSTNESS")
print(robustness.to_string(index=False))

print("\nPATHWAY COHERENCE")
print(coherence.to_string(index=False))

print("\nExported principal signature tables:")
for filename, df in exports.items():
    print(f"{filename}: {len(df)} genes")
