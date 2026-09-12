import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr, binomtest

FROZEN = Path(
    "results/DS001_GSE38900/integrated_architecture/tables/"
    "Figure7A_high_confidence_gene_architecture.tsv"
)

DE_DIR = Path(
    "results/rnaseq_validation/influenza_SIG/DE"
)

OUTDIR = Path(
    "results/rnaseq_validation/influenza_SIG/validation"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

SEASONS = {
    "2018": "GSE158592_influenza_vs_control_limma_unadjusted.tsv",
    "2019": "GSE155635_influenza_vs_control_limma_unadjusted.tsv",
    "2020": "GSE196350_influenza_vs_control_limma_unadjusted.tsv",
    "2022": "GSE213168_influenza_vs_control_limma_unadjusted.tsv",
}

frozen = pd.read_csv(FROZEN, sep="\t")

print("=== FROZEN TABLE ===")
print("Rows:", len(frozen))
print("Classes:")
print(frozen["Figure7A_class"].value_counts())

# ------------------------------------------------------------
# Identify authoritative frozen influenza effect-size column
# ------------------------------------------------------------

candidates = [
    "flu_logFC",
    "influenza6884_logFC",
    "flu6884_logFC",
    "influenza_logFC",
    "Influenza_logFC",
    "influenza_vs_control_logFC",
    "inf6884_logFC",
]

influenza_col = None

for c in candidates:
    if c in frozen.columns:
        influenza_col = c
        break

if influenza_col is None:
    # Conservative fallback: identify columns containing both
    # influenza/flu and logFC.
    matches = [
        c for c in frozen.columns
        if (
            ("influenza" in c.lower() or "flu" in c.lower())
            and "logfc" in c.lower()
        )
    ]

    if len(matches) == 1:
        influenza_col = matches[0]
    else:
        print("\nAvailable columns:")
        for c in frozen.columns:
            print(c)
        raise RuntimeError(
            "Could not uniquely identify frozen influenza logFC column."
        )

print("\nFrozen influenza reference column:")
print(influenza_col)

assert len(frozen) == 170
assert frozen["gene_symbol"].is_unique

# ------------------------------------------------------------
# Load seasonal limma results
# ------------------------------------------------------------

merged = frozen.copy()

season_summaries = []

for season, filename in SEASONS.items():

    d = pd.read_csv(DE_DIR / filename, sep="\t")

    if "gene_symbol" not in d.columns:
        raise RuntimeError(f"{season}: gene_symbol missing")

    # Should already be unique because seasonal expression matrices
    # were collapsed before limma.
    if d["gene_symbol"].duplicated().any():
        d = (
            d.sort_values("P.Value")
             .drop_duplicates("gene_symbol", keep="first")
        )

    rename = {
        "logFC": f"logFC_{season}",
        "AveExpr": f"AveExpr_{season}",
        "t": f"t_{season}",
        "P.Value": f"P_{season}",
        "adj.P.Val": f"FDR_{season}",
        "B": f"B_{season}",
    }

    d = d.rename(columns=rename)

    keep = [
        "gene_symbol",
        f"logFC_{season}",
        f"AveExpr_{season}",
        f"t_{season}",
        f"P_{season}",
        f"FDR_{season}",
    ]

    merged = merged.merge(
        d[keep],
        on="gene_symbol",
        how="left",
        validate="one_to_one"
    )

    merged[f"evaluable_{season}"] = (
        merged[f"logFC_{season}"].notna()
    )

    merged[f"direction_match_{season}"] = np.where(
        merged[f"evaluable_{season}"],
        np.sign(merged[f"logFC_{season}"]) ==
        np.sign(merged[influenza_col]),
        np.nan
    )

    merged[f"FDR05_{season}"] = (
        merged[f"FDR_{season}"].notna() &
        (merged[f"FDR_{season}"] < 0.05)
    )

    merged[f"FDR05_direction_match_{season}"] = (
        merged[f"FDR05_{season}"] &
        (merged[f"direction_match_{season}"] == True)
    )

# ------------------------------------------------------------
# Cross-season consistency
# ------------------------------------------------------------

direction_cols = [
    f"direction_match_{s}"
    for s in SEASONS
]

merged["n_seasons_evaluable"] = (
    merged[
        [f"evaluable_{s}" for s in SEASONS]
    ].sum(axis=1)
)

merged["n_seasons_direction_match"] = (
    merged[direction_cols]
    .fillna(False)
    .sum(axis=1)
)

merged["all_evaluable_seasons_match"] = (
    merged["n_seasons_direction_match"] ==
    merged["n_seasons_evaluable"]
)

merged["match_3plus_seasons"] = (
    merged["n_seasons_direction_match"] >= 3
)

merged["match_all_4_seasons"] = (
    (merged["n_seasons_evaluable"] == 4) &
    (merged["n_seasons_direction_match"] == 4)
)

# Mean RNA-seq effect across seasons
merged["mean_RNAseq_logFC"] = merged[
    [f"logFC_{s}" for s in SEASONS]
].mean(axis=1, skipna=True)

merged["median_RNAseq_logFC"] = merged[
    [f"logFC_{s}" for s in SEASONS]
].median(axis=1, skipna=True)

# ------------------------------------------------------------
# Summary function
# ------------------------------------------------------------

def summarize(sub, label):

    rows = []

    for season in SEASONS:

        e = sub[
            sub[f"evaluable_{season}"]
        ].copy()

        n = len(e)

        if n == 0:
            continue

        conc = int(
            (e[f"direction_match_{season}"] == True).sum()
        )

        pr = (np.nan, np.nan)
        sr = (np.nan, np.nan)

        if n >= 3:
            pr = pearsonr(
                e[influenza_col],
                e[f"logFC_{season}"]
            )

            sr = spearmanr(
                e[influenza_col],
                e[f"logFC_{season}"]
            )

        rows.append({
            "class": label,
            "season": season,
            "frozen_genes": len(sub),
            "evaluable": n,
            "direction_concordant": conc,
            "direction_concordance_fraction": conc / n,
            "direction_binomial_p":
                binomtest(
                    conc,
                    n,
                    p=0.5,
                    alternative="greater"
                ).pvalue,
            "RNAseq_FDR05":
                int(e[f"FDR05_{season}"].sum()),
            "RNAseq_FDR05_direction_match":
                int(
                    e[
                        f"FDR05_direction_match_{season}"
                    ].sum()
                ),
            "pearson_r": pr[0],
            "pearson_p": pr[1],
            "spearman_rho": sr[0],
            "spearman_p": sr[1],
            "median_RNAseq_logFC":
                e[f"logFC_{season}"].median(),
            "median_microarray_influenza_logFC":
                e[influenza_col].median(),
        })

    return rows


summary_rows = []

summary_rows += summarize(
    merged,
    "ALL_170"
)

for cls, sub in merged.groupby(
    "Figure7A_class",
    sort=False
):
    summary_rows += summarize(sub, cls)

summary = pd.DataFrame(summary_rows)

# ------------------------------------------------------------
# Cross-season summary by class
# ------------------------------------------------------------

cross_rows = []

groups = [("ALL_170", merged)] + list(
    merged.groupby("Figure7A_class", sort=False)
)

for label, sub in groups:

    complete = sub[
        sub["n_seasons_evaluable"] == 4
    ].copy()

    cross_rows.append({
        "class": label,
        "frozen_genes": len(sub),
        "evaluable_all_4_seasons": len(complete),

        "direction_match_all_4":
            int(complete["match_all_4_seasons"].sum()),

        "direction_match_all_4_fraction":
            (
                complete["match_all_4_seasons"].mean()
                if len(complete) else np.nan
            ),

        "direction_match_3plus":
            int(
                (
                    complete["n_seasons_direction_match"] >= 3
                ).sum()
            ),

        "direction_match_3plus_fraction":
            (
                (
                    complete["n_seasons_direction_match"] >= 3
                ).mean()
                if len(complete) else np.nan
            ),

        "median_number_matching_seasons":
            (
                complete[
                    "n_seasons_direction_match"
                ].median()
                if len(complete) else np.nan
            ),
    })

cross = pd.DataFrame(cross_rows)

# ------------------------------------------------------------
# Mean effect correlation across seasons
# ------------------------------------------------------------

effect_rows = []

for label, sub in groups:

    e = sub.dropna(
        subset=[
            influenza_col,
            "mean_RNAseq_logFC"
        ]
    )

    if len(e) >= 3:

        pr = pearsonr(
            e[influenza_col],
            e["mean_RNAseq_logFC"]
        )

        sr = spearmanr(
            e[influenza_col],
            e["mean_RNAseq_logFC"]
        )

    else:
        pr = (np.nan, np.nan)
        sr = (np.nan, np.nan)

    effect_rows.append({
        "class": label,
        "evaluable": len(e),
        "pearson_r_microarray_vs_mean_RNAseq":
            pr[0],
        "pearson_p": pr[1],
        "spearman_rho_microarray_vs_mean_RNAseq":
            sr[0],
        "spearman_p": sr[1],
    })

effects = pd.DataFrame(effect_rows)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

master_file = (
    OUTDIR /
    "SIG_four_season_frozen170_validation_master.tsv"
)

summary_file = (
    OUTDIR /
    "SIG_four_season_frozen170_seasonal_summary.tsv"
)

cross_file = (
    OUTDIR /
    "SIG_four_season_frozen170_crossseason_summary.tsv"
)

effect_file = (
    OUTDIR /
    "SIG_four_season_frozen170_effect_correlation.tsv"
)

merged.to_csv(
    master_file,
    sep="\t",
    index=False
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)

cross.to_csv(
    cross_file,
    sep="\t",
    index=False
)

effects.to_csv(
    effect_file,
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

print("\n========================================")
print("SEASONAL VALIDATION")
print("========================================")

print(summary.to_string(index=False))

print("\n========================================")
print("CROSS-SEASON CONSISTENCY")
print("========================================")

print(cross.to_string(index=False))

print("\n========================================")
print("MEAN RNA-seq EFFECT CORRELATION")
print("========================================")

print(effects.to_string(index=False))

print("\n========================================")
print("TOP CONSISTENT GENES")
print("========================================")

cols = [
    "gene_symbol",
    "Figure7A_class",
    influenza_col,
    "logFC_2018",
    "logFC_2019",
    "logFC_2020",
    "logFC_2022",
    "n_seasons_direction_match",
]

print(
    merged.sort_values(
        [
            "n_seasons_direction_match",
            "mean_RNAseq_logFC"
        ],
        ascending=[False, False]
    )[cols]
    .head(40)
    .to_string(index=False)
)

print("\nWritten:")
for x in [
    master_file,
    summary_file,
    cross_file,
    effect_file
]:
    print(x)
