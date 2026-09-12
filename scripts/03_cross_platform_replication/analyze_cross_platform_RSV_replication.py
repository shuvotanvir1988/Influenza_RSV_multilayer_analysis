#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.stats import hypergeom

ROOT = Path.home() / "Influenza_RSV_Project"

TABLE_DIR = (
    ROOT /
    "results/DS001_GSE38900/"
    "differential_expression/tables"
)

f10558 = (
    TABLE_DIR /
    "GPL10558_RSV_acute_vs_control_full_DE.tsv"
)

f6884 = (
    TABLE_DIR /
    "GPL6884_RSV_acute_vs_control_full_DE.tsv"
)

a = pd.read_csv(f10558, sep="\t")
b = pd.read_csv(f6884, sep="\t")

a = a[
    [
        "gene_symbol",
        "logFC",
        "P.Value",
        "adj.P.Val"
    ]
].copy()

b = b[
    [
        "gene_symbol",
        "logFC",
        "P.Value",
        "adj.P.Val"
    ]
].copy()

a.columns = [
    "gene_symbol",
    "logFC_GPL10558",
    "P_GPL10558",
    "FDR_GPL10558"
]

b.columns = [
    "gene_symbol",
    "logFC_GPL6884",
    "P_GPL6884",
    "FDR_GPL6884"
]

m = a.merge(
    b,
    on="gene_symbol",
    how="inner",
    validate="one_to_one"
)

print("=" * 80)
print("SESSION 16 — CROSS-PLATFORM RSV REPLICATION")
print("=" * 80)

print("\nShared genes:", len(m))

# ============================================================
# Global effect-size correlations
# ============================================================

pearson = pearsonr(
    m["logFC_GPL10558"],
    m["logFC_GPL6884"]
)

spearman = spearmanr(
    m["logFC_GPL10558"],
    m["logFC_GPL6884"]
)

print("\nGLOBAL logFC CORRELATION")
print("-" * 80)

print(
    "Pearson r:",
    pearson.statistic,
    "P:",
    pearson.pvalue
)

print(
    "Spearman rho:",
    spearman.statistic,
    "P:",
    spearman.pvalue
)

# ============================================================
# Direction concordance
# ============================================================

nonzero = m[
    (m["logFC_GPL10558"] != 0) &
    (m["logFC_GPL6884"] != 0)
].copy()

nonzero["same_direction"] = (
    np.sign(nonzero["logFC_GPL10558"]) ==
    np.sign(nonzero["logFC_GPL6884"])
)

direction_pct = (
    100 *
    nonzero["same_direction"].mean()
)

print("\nGLOBAL DIRECTION CONCORDANCE")
print("-" * 80)

print(
    "Same direction:",
    int(nonzero["same_direction"].sum()),
    "/",
    len(nonzero)
)

print(
    "Percent:",
    direction_pct
)

# ============================================================
# Significant sets
# ============================================================

sig10558 = set(
    m.loc[
        m["FDR_GPL10558"] < 0.05,
        "gene_symbol"
    ]
)

sig6884 = set(
    m.loc[
        m["FDR_GPL6884"] < 0.05,
        "gene_symbol"
    ]
)

overlap = sig10558 & sig6884

print("\nFDR < 0.05 REPLICATION")
print("-" * 80)

print("GPL10558 significant:", len(sig10558))
print("GPL6884 significant:", len(sig6884))
print("Significant in both:", len(overlap))

if len(sig10558) > 0:

    replication_rate = (
        100 * len(overlap) / len(sig10558)
    )

    print(
        "GPL10558 replication rate in GPL6884:",
        replication_rate,
        "%"
    )

# ============================================================
# Hypergeometric overlap
# ============================================================

M = len(m)
n = len(sig6884)
N = len(sig10558)
k = len(overlap)

if N > 0 and n > 0:

    overlap_p = hypergeom.sf(
        k - 1,
        M,
        n,
        N
    )

    print(
        "Hypergeometric overlap P:",
        overlap_p
    )

# ============================================================
# Direction among replicated significant genes
# ============================================================

both = m[
    m["gene_symbol"].isin(overlap)
].copy()

if len(both) > 0:

    both["same_direction"] = (
        np.sign(both["logFC_GPL10558"]) ==
        np.sign(both["logFC_GPL6884"])
    )

    print(
        "\nDirection concordance among genes "
        "significant on BOTH platforms:"
    )

    print(
        int(both["same_direction"].sum()),
        "/",
        len(both),
        "=",
        100 * both["same_direction"].mean(),
        "%"
    )

# ============================================================
# Stronger-effect replication
# ============================================================

strong10558 = set(
    m.loc[
        (m["FDR_GPL10558"] < 0.05) &
        (m["logFC_GPL10558"].abs() >= 0.5),
        "gene_symbol"
    ]
)

strong6884 = set(
    m.loc[
        (m["FDR_GPL6884"] < 0.05) &
        (m["logFC_GPL6884"].abs() >= 0.5),
        "gene_symbol"
    ]
)

strong_overlap = (
    strong10558 &
    strong6884
)

print("\nFDR < 0.05 + |logFC| >= 0.5")
print("-" * 80)

print("GPL10558:", len(strong10558))
print("GPL6884:", len(strong6884))
print("Overlap:", len(strong_overlap))

# ============================================================
# Save full cross-platform table
# ============================================================

m["significant_GPL10558"] = (
    m["FDR_GPL10558"] < 0.05
)

m["significant_GPL6884"] = (
    m["FDR_GPL6884"] < 0.05
)

m["significant_both"] = (
    m["significant_GPL10558"] &
    m["significant_GPL6884"]
)

m["same_direction"] = (
    np.sign(m["logFC_GPL10558"]) ==
    np.sign(m["logFC_GPL6884"])
)

m["abs_logFC_GPL10558"] = (
    m["logFC_GPL10558"].abs()
)

m["abs_logFC_GPL6884"] = (
    m["logFC_GPL6884"].abs()
)

m["minimum_abs_logFC"] = np.minimum(
    m["abs_logFC_GPL10558"],
    m["abs_logFC_GPL6884"]
)

m = m.sort_values(
    [
        "significant_both",
        "same_direction",
        "minimum_abs_logFC"
    ],
    ascending=[
        False,
        False,
        False
    ]
)

outfile = (
    TABLE_DIR /
    "DS001_GSE38900_RSV_cross_platform_replication.tsv"
)

m.to_csv(
    outfile,
    sep="\t",
    index=False
)

# ============================================================
# Save replicated genes
# ============================================================

replicated = m[
    m["significant_both"] &
    m["same_direction"]
].copy()

replicated_file = (
    TABLE_DIR /
    "DS001_GSE38900_RSV_replicated_DE_genes.tsv"
)

replicated.to_csv(
    replicated_file,
    sep="\t",
    index=False
)

print("\nTOP REPLICATED GENES")
print("-" * 80)

if len(replicated) > 0:

    print(
        replicated[
            [
                "gene_symbol",
                "logFC_GPL10558",
                "FDR_GPL10558",
                "logFC_GPL6884",
                "FDR_GPL6884"
            ]
        ]
        .head(30)
        .to_string(index=False)
    )

print("\nSaved:")
print(outfile)
print(replicated_file)

print(
    "\nSESSION 16 CROSS-PLATFORM "
    "RSV REPLICATION COMPLETE"
)
