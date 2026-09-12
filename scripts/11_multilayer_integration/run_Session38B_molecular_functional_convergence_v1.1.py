from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import (
    pearsonr,
    spearmanr,
    binomtest,
)


def benjamini_hochberg(pvalues, alpha=0.05):
    """Benjamini-Hochberg FDR adjustment without external dependency."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)

    order = np.argsort(p)
    ranked = p[order]

    adjusted_ranked = ranked * n / np.arange(1, n + 1)

    # Enforce monotonicity from largest rank toward smallest.
    adjusted_ranked = np.minimum.accumulate(
        adjusted_ranked[::-1]
    )[::-1]

    adjusted_ranked = np.minimum(adjusted_ranked, 1.0)

    adjusted = np.empty(n, dtype=float)
    adjusted[order] = adjusted_ranked

    reject = adjusted <= alpha

    return reject, adjusted


ROOT = Path.home() / "Influenza_RSV_Project"

INFILE = ROOT / (
    "results/influenza_multilayer_integration/gene_level/"
    "SESSION27_INFLUENZA_OBJECTIVE_CANDIDATE_CLASSIFICATION_v1.0.tsv"
)

OUT = (
    ROOT /
    "results/session38_additional_analysis/Session38B/tables"
)
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INFILE, sep="\t")

assert len(df) == 170
assert df["gene_symbol"].nunique() == 170

# ============================================================
# DERIVED VARIABLES
# ============================================================

df["discovery_rna_protein_same_sign"] = np.where(
    df["proteomics_evaluable"],
    (
        df["discovery_iav_logFC"]
        * df["proteomics_median_logFC"]
        > 0
    ),
    np.nan
)

df["validation_rna_protein_same_sign"] = np.where(
    df["proteomics_evaluable"]
    & df["rnaseq_mean_logFC"].notna(),
    (
        df["rnaseq_mean_logFC"]
        * df["proteomics_median_logFC"]
        > 0
    ),
    np.nan
)

df["crispr_dependency_strength"] = np.where(
    df["crispr_iav_evaluable"],
    1.0 - df["crispr_iav_percentile"],
    np.nan
)

# ============================================================
# OBJECTIVE 1 — RNA/PROTEIN
# ============================================================

rp_discovery = df[
    df["proteomics_evaluable"].astype(bool)
    & df["discovery_iav_logFC"].notna()
    & df["proteomics_median_logFC"].notna()
].copy()

rp_validation = df[
    df["proteomics_evaluable"].astype(bool)
    & df["rnaseq_mean_logFC"].notna()
    & df["proteomics_median_logFC"].notna()
].copy()

results = []

def add_corr_results(label, x, y):
    pear_r, pear_p = pearsonr(x, y)
    spear_r, spear_p = spearmanr(x, y)

    results.extend([
        {
            "analysis": label,
            "test": "Pearson",
            "n": len(x),
            "effect": pear_r,
            "p_value": pear_p,
        },
        {
            "analysis": label,
            "test": "Spearman",
            "n": len(x),
            "effect": spear_r,
            "p_value": spear_p,
        },
    ])

add_corr_results(
    "Discovery_RNA_vs_protein",
    rp_discovery["discovery_iav_logFC"],
    rp_discovery["proteomics_median_logFC"],
)

add_corr_results(
    "Independent_RNA_vs_protein",
    rp_validation["rnaseq_mean_logFC"],
    rp_validation["proteomics_median_logFC"],
)

# Sign concordance exact tests.
for label, dat, col in [
    (
        "Discovery_RNA_vs_protein_sign",
        rp_discovery,
        "discovery_rna_protein_same_sign",
    ),
    (
        "Independent_RNA_vs_protein_sign",
        rp_validation,
        "validation_rna_protein_same_sign",
    ),
]:
    same = int(dat[col].astype(bool).sum())
    n = len(dat)

    bt = binomtest(
        same,
        n,
        p=0.5,
        alternative="two-sided",
    )

    results.append({
        "analysis": label,
        "test": "Exact_binomial_vs_0.5",
        "n": n,
        "effect": same / n,
        "p_value": bt.pvalue,
    })

# ============================================================
# OBJECTIVE 2 — MOLECULAR RESPONSE VS CRISPR
# ============================================================

crispr_tests = []

def crispr_spearman(label, data, molecular_col):
    z = data[
        data[molecular_col].notna()
        & data["crispr_dependency_strength"].notna()
    ].copy()

    rho, p = spearmanr(
        z[molecular_col].abs(),
        z["crispr_dependency_strength"],
    )

    crispr_tests.append({
        "analysis": label,
        "test": "Spearman",
        "n": len(z),
        "effect": rho,
        "p_value": p,
    })

crispr_eval = df[
    df["crispr_iav_evaluable"].astype(bool)
].copy()

crispr_spearman(
    "Abs_discovery_RNA_vs_CRISPR_dependency",
    crispr_eval,
    "discovery_iav_logFC",
)

crispr_spearman(
    "Abs_independent_RNA_vs_CRISPR_dependency",
    crispr_eval,
    "rnaseq_mean_logFC",
)

both = df[
    df["crispr_iav_evaluable"].astype(bool)
    & df["proteomics_evaluable"].astype(bool)
].copy()

crispr_spearman(
    "Abs_protein_vs_CRISPR_dependency",
    both,
    "proteomics_median_logFC",
)

pvals = [x["p_value"] for x in crispr_tests]

reject, qvals = benjamini_hochberg(
    pvals,
    alpha=0.05,
)

for row, q, rej in zip(
    crispr_tests,
    qvals,
    reject,
):
    row["FDR"] = q
    row["FDR_significant"] = bool(rej)

results.extend(crispr_tests)

# ============================================================
# OBJECTIVE 3 — CONVERGENCE CLASSES
# ============================================================

def classify(row):

    crispr = bool(
        row["crispr_iav_dependency_support"]
    )

    protein_eval = bool(
        row["proteomics_evaluable"]
    )

    if crispr and not protein_eval:
        return "RNA_CRISPR_NO_PROTEIN_COVERAGE"

    if protein_eval:
        concordant = bool(
            row["discovery_rna_protein_same_sign"]
        )

        if crispr and concordant:
            return "RNA_PROTEIN_CRISPR_CONVERGENT"

        if crispr and not concordant:
            return "RNA_CRISPR_PROTEIN_DISCORDANT"

        if (not crispr) and concordant:
            return (
                "RNA_PROTEIN_CONCORDANT_"
                "NO_CRISPR_SUPPORT"
            )

        if (not crispr) and not concordant:
            return (
                "RNA_PROTEIN_DISCORDANT_"
                "NO_CRISPR_SUPPORT"
            )

    return "OTHER"

df["session38B_convergence_class"] = (
    df.apply(classify, axis=1)
)

# ============================================================
# OUTPUT
# ============================================================

stats = pd.DataFrame(results)

master_cols = [
    "gene_symbol",
    "proteomics_architecture_class",
    "discovery_iav_logFC",
    "rnaseq_mean_logFC",
    "rnaseq_median_logFC",
    "proteomics_evaluable",
    "proteomics_median_logFC",
    "proteomics_median_standardized_effect",
    "proteomics_validation_category",
    "discovery_rna_protein_same_sign",
    "validation_rna_protein_same_sign",
    "crispr_iav_evaluable",
    "crispr_iav_percentile",
    "crispr_dependency_strength",
    "crispr_iav_dependency_support",
    "crispr_iav_dependency_high_confidence",
    "session38B_convergence_class",
]

master = df[master_cols].copy()

master.to_csv(
    OUT /
    "Session38B_MOLECULAR_FUNCTIONAL_CONVERGENCE_MASTER_v1.0.tsv",
    sep="\t",
    index=False,
)

stats.to_csv(
    OUT /
    "Session38B_STATISTICAL_TESTS_v1.0.tsv",
    sep="\t",
    index=False,
)

class_summary = (
    master["session38B_convergence_class"]
    .value_counts()
    .rename_axis("convergence_class")
    .reset_index(name="n_genes")
)

class_summary.to_csv(
    OUT /
    "Session38B_CONVERGENCE_CLASS_SUMMARY_v1.0.tsv",
    sep="\t",
    index=False,
)

three_layer = master[
    master["session38B_convergence_class"].isin([
        "RNA_PROTEIN_CRISPR_CONVERGENT",
        "RNA_CRISPR_PROTEIN_DISCORDANT",
        "RNA_CRISPR_NO_PROTEIN_COVERAGE",
    ])
].copy()

three_layer.to_csv(
    OUT /
    "Session38B_CRISPR_BRIDGE_MOLECULAR_CONTEXT_v1.0.tsv",
    sep="\t",
    index=False,
)

# ============================================================
# TERMINAL REPORT
# ============================================================

print("=== SESSION 38B COMPLETE ===")
print()

print("=== STATISTICAL TESTS ===")
print(stats.to_string(index=False))

print()
print("=== CONVERGENCE CLASS COUNTS ===")
print(class_summary.to_string(index=False))

print()
print("=== CRISPR BRIDGE MOLECULAR CONTEXT ===")

show = [
    "gene_symbol",
    "discovery_iav_logFC",
    "rnaseq_mean_logFC",
    "proteomics_evaluable",
    "proteomics_median_logFC",
    "proteomics_validation_category",
    "crispr_iav_percentile",
    "crispr_iav_dependency_high_confidence",
    "session38B_convergence_class",
]

print(
    three_layer[show]
    .sort_values(
        [
            "session38B_convergence_class",
            "crispr_iav_percentile",
        ]
    )
    .to_string(index=False)
)

print()
print("=== RNA/PROTEIN CONCORDANCE ===")

print(
    "Discovery:",
    int(
        rp_discovery[
            "discovery_rna_protein_same_sign"
        ].astype(bool).sum()
    ),
    "/",
    len(rp_discovery),
)

print(
    "Independent RNA:",
    int(
        rp_validation[
            "validation_rna_protein_same_sign"
        ].astype(bool).sum()
    ),
    "/",
    len(rp_validation),
)

print()
print("Written outputs:")
for p in [
    "Session38B_MOLECULAR_FUNCTIONAL_CONVERGENCE_MASTER_v1.0.tsv",
    "Session38B_STATISTICAL_TESTS_v1.0.tsv",
    "Session38B_CONVERGENCE_CLASS_SUMMARY_v1.0.tsv",
    "Session38B_CRISPR_BRIDGE_MOLECULAR_CONTEXT_v1.0.tsv",
]:
    print(OUT / p)
