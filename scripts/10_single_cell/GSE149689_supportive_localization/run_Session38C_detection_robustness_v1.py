from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest

ROOT = Path.home() / "Influenza_RSV_Project"

TAB = (
    ROOT /
    "results/session38_additional_analysis/Session38C/tables"
)

OUT = TAB

long = pd.read_csv(
    TAB / "Session38C_PSEUDOBULK_RAW_COUNTS_LONG_v1.0.tsv",
    sep="\t"
)

effects = pd.read_csv(
    TAB / "Session38C_GENE_LEVEL_CELLTYPE_EFFECTS_v1.0.tsv",
    sep="\t"
)

architecture_primary = pd.read_csv(
    TAB / "Session38C_ARCHITECTURE_LOCALIZATION_v1.0.tsv",
    sep="\t"
)

celltypes = architecture_primary["cell_type"].astype(str).tolist()

assert len(celltypes) == 5

flu_donors = {
    "Flu 1",
    "Flu 2",
    "Flu 3",
    "Flu 4",
    "Flu 5",
}

normal_donors = {
    "Normal 1",
    "Normal 2",
    "Normal 3",
    "Normal 4",
}

# ------------------------------------------------------------
# BH helper
# ------------------------------------------------------------

def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)

    order = np.argsort(p)
    ranked = p[order]

    q_ranked = ranked * n / np.arange(1, n + 1)

    q_ranked = np.minimum.accumulate(
        q_ranked[::-1]
    )[::-1]

    q_ranked = np.minimum(q_ranked, 1.0)

    q = np.empty(n, dtype=float)
    q[order] = q_ranked

    return q


# ------------------------------------------------------------
# Detection eligibility
# ------------------------------------------------------------

detection_rows = []

for ct in celltypes:

    z = long[
        long["cell_type"] == ct
    ].copy()

    genes = sorted(z["gene_symbol"].unique())

    for gene in genes:

        g = z[
            z["gene_symbol"] == gene
        ].copy()

        flu = g[
            g["donor_id"].isin(flu_donors)
        ]

        normal = g[
            g["donor_id"].isin(normal_donors)
        ]

        assert len(flu) == 5
        assert len(normal) == 4

        flu_detected = int(
            (flu["raw_count"] > 0).sum()
        )

        normal_detected = int(
            (normal["raw_count"] > 0).sum()
        )

        eligible = (
            flu_detected >= 3
            and normal_detected >= 3
        )

        detection_rows.append({
            "cell_type": ct,
            "gene_symbol": gene,
            "influenza_donors_detected_n":
                flu_detected,
            "normal_donors_detected_n":
                normal_detected,
            "detection_eligible":
                eligible,
        })

detection = pd.DataFrame(detection_rows)

assert len(detection) == 5 * 161

# ------------------------------------------------------------
# Recalculate architecture metrics
# ------------------------------------------------------------

rows = []

for ct in celltypes:

    elig = detection[
        (detection["cell_type"] == ct)
        & detection["detection_eligible"]
    ][
        ["gene_symbol"]
    ]

    z = effects[
        effects["cell_type"] == ct
    ].merge(
        elig,
        on="gene_symbol",
        how="inner"
    )

    n = len(z)

    if n < 20:
        raise RuntimeError(
            f"{ct}: too few genes after detection filter ({n})"
        )

    rho, p_rho = spearmanr(
        z["sc_effect_logCPM"],
        z["bulk_flu_logFC"]
    )

    same = int(z["same_direction"].sum())

    bt = binomtest(
        same,
        n,
        p=0.5,
        alternative="two-sided"
    )

    primary = architecture_primary[
        architecture_primary["cell_type"] == ct
    ].iloc[0]

    rows.append({
        "cell_type": ct,
        "primary_n_genes":
            int(primary["n_genes"]),
        "robust_n_genes":
            n,
        "genes_retained_fraction":
            n / int(primary["n_genes"]),
        "primary_spearman_rho":
            float(primary["spearman_rho"]),
        "robust_spearman_rho":
            rho,
        "robust_spearman_p":
            p_rho,
        "primary_direction_concordance":
            float(primary["direction_concordance"]),
        "robust_direction_concordant_n":
            same,
        "robust_direction_total_n":
            n,
        "robust_direction_concordance":
            same / n,
        "robust_direction_binomial_p":
            bt.pvalue,
    })

robust = pd.DataFrame(rows)

robust["robust_spearman_FDR"] = bh_adjust(
    robust["robust_spearman_p"].to_numpy()
)

robust["robust_direction_FDR"] = bh_adjust(
    robust["robust_direction_binomial_p"].to_numpy()
)

# ------------------------------------------------------------
# Rank stability
# ------------------------------------------------------------

robust["primary_rho_rank"] = (
    robust["primary_spearman_rho"]
    .rank(
        method="min",
        ascending=False
    )
    .astype(int)
)

robust["robust_rho_rank"] = (
    robust["robust_spearman_rho"]
    .rank(
        method="min",
        ascending=False
    )
    .astype(int)
)

robust["rho_rank_change"] = (
    robust["robust_rho_rank"]
    - robust["primary_rho_rank"]
)

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------

detection.to_csv(
    OUT /
    "Session38C_DETECTION_ELIGIBILITY_BY_GENE_CELLTYPE_v1.0.tsv",
    sep="\t",
    index=False
)

robust.to_csv(
    OUT /
    "Session38C_ARCHITECTURE_DETECTION_ROBUSTNESS_v1.0.tsv",
    sep="\t",
    index=False
)

print("=== SESSION 38C DETECTION ROBUSTNESS ===")
print()

print(
    robust.sort_values(
        "robust_spearman_rho",
        ascending=False
    ).to_string(index=False)
)

print()
print("=== ELIGIBLE GENE COUNTS ===")

for ct in celltypes:

    z = detection[
        (detection["cell_type"] == ct)
        & detection["detection_eligible"]
    ]

    print(
        f"{ct}: "
        f"{len(z)}/161 "
        f"({len(z)/161:.1%})"
    )

print()
print("=== PRIMARY VS ROBUST RANK ===")

print(
    robust[
        [
            "cell_type",
            "primary_rho_rank",
            "robust_rho_rank",
            "rho_rank_change",
            "primary_spearman_rho",
            "robust_spearman_rho",
        ]
    ]
    .sort_values("robust_rho_rank")
    .to_string(index=False)
)

print()
print("Written:")
print(
    OUT /
    "Session38C_DETECTION_ELIGIBILITY_BY_GENE_CELLTYPE_v1.0.tsv"
)
print(
    OUT /
    "Session38C_ARCHITECTURE_DETECTION_ROBUSTNESS_v1.0.tsv"
)
