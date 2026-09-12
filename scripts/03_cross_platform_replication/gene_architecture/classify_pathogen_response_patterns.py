#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("results/DS001_GSE38900/signature_analysis")
TABLES = ROOT / "tables"

MASTER = TABLES / "DS001_GSE38900_gene_evidence_master.tsv"

FDR = 0.05

x = pd.read_csv(MASTER, sep="\t")


# ---------------------------------------------------------
# Infection-vs-control state assignment
# ---------------------------------------------------------

def infection_state(logfc, fdr):

    if pd.isna(logfc) or pd.isna(fdr):
        return "not_evaluable"

    if fdr >= FDR:
        return "not_significant"

    if logfc > 0:
        return "up"

    if logfc < 0:
        return "down"

    return "zero"


x["flu_state"] = [
    infection_state(fc, fdr)
    for fc, fdr in zip(x["flu_logFC"], x["flu_FDR"])
]

x["rsv_state"] = [
    infection_state(fc, fdr)
    for fc, fdr in zip(x["rsv6884_logFC"], x["rsv6884_FDR"])
]


# ---------------------------------------------------------
# Detailed direct-differential response classification
# ---------------------------------------------------------

def classify_pattern(row):

    # Only interpret genes with significant direct contrast.
    if not row["direct_significant"]:
        return "no_significant_pathogen_difference"

    fs = row["flu_state"]
    rs = row["rsv_state"]
    direct = row["flu_vs_rsv_logFC"]

    # -----------------------------------------------------
    # Opposite-direction responses
    # -----------------------------------------------------

    if fs == "up" and rs == "down":
        return "opposite_influenza_up_RSV_down"

    if fs == "down" and rs == "up":
        return "opposite_influenza_down_RSV_up"

    # -----------------------------------------------------
    # Both induced
    # -----------------------------------------------------

    if fs == "up" and rs == "up":

        if direct > 0:
            return "shared_induction_influenza_stronger"

        if direct < 0:
            return "shared_induction_RSV_stronger"

    # -----------------------------------------------------
    # Both suppressed
    # -----------------------------------------------------

    if fs == "down" and rs == "down":

        # Positive Influenza-RSV means influenza is less
        # negative / RSV is more strongly suppressed.
        if direct > 0:
            return "shared_suppression_RSV_stronger"

        if direct < 0:
            return "shared_suppression_influenza_stronger"

    # -----------------------------------------------------
    # Influenza-only significant response
    # -----------------------------------------------------

    if fs == "up" and rs == "not_significant":
        return "influenza_specific_induction"

    if fs == "down" and rs == "not_significant":
        return "influenza_specific_suppression"

    # -----------------------------------------------------
    # RSV-only significant response
    # -----------------------------------------------------

    if fs == "not_significant" and rs == "up":
        return "RSV_specific_induction"

    if fs == "not_significant" and rs == "down":
        return "RSV_specific_suppression"

    # -----------------------------------------------------
    # Neither infection-vs-control contrast significant,
    # despite significant direct contrast.
    #
    # This can occur because the direct comparison has
    # greater precision than either control comparison.
    # -----------------------------------------------------

    if fs == "not_significant" and rs == "not_significant":

        if direct > 0:
            return "direct_only_influenza_higher"

        if direct < 0:
            return "direct_only_RSV_higher"

    return "other_direct_pattern"


x["pathogen_response_pattern"] = x.apply(
    classify_pattern,
    axis=1
)


# ---------------------------------------------------------
# Add broad interpretive class
# ---------------------------------------------------------

def broad_class(pattern):

    mapping = {

        "shared_induction_influenza_stronger":
            "shared_response_quantitatively_influenza_biased",

        "shared_induction_RSV_stronger":
            "shared_response_quantitatively_RSV_biased",

        "shared_suppression_influenza_stronger":
            "shared_response_quantitatively_influenza_biased",

        "shared_suppression_RSV_stronger":
            "shared_response_quantitatively_RSV_biased",

        "influenza_specific_induction":
            "influenza_specific",

        "influenza_specific_suppression":
            "influenza_specific",

        "RSV_specific_induction":
            "RSV_specific",

        "RSV_specific_suppression":
            "RSV_specific",

        "opposite_influenza_up_RSV_down":
            "opposite_direction",

        "opposite_influenza_down_RSV_up":
            "opposite_direction",

        "direct_only_influenza_higher":
            "direct_only",

        "direct_only_RSV_higher":
            "direct_only",

        "no_significant_pathogen_difference":
            "no_significant_pathogen_difference",
    }

    return mapping.get(pattern, "other")


x["pathogen_response_broad_class"] = \
    x["pathogen_response_pattern"].map(broad_class)


# ---------------------------------------------------------
# Export updated master table
# ---------------------------------------------------------

updated = (
    TABLES /
    "DS001_GSE38900_gene_evidence_master_with_response_patterns.tsv"
)

x.to_csv(updated, sep="\t", index=False)


# ---------------------------------------------------------
# Detailed pattern summary
# ---------------------------------------------------------

detail = (
    x.loc[x["direct_significant"]]
    .groupby("pathogen_response_pattern", dropna=False)
    .agg(
        n_genes=("gene_symbol", "size"),
        leading_edge_genes=("session17_leading_edge", "sum"),
        rsv_replicated_genes=("rsv_replicated", "sum"),
        shared_replicated_genes=("shared_replicated", "sum"),
        median_flu_logFC=("flu_logFC", "median"),
        median_rsv_logFC=("rsv6884_logFC", "median"),
        median_direct_logFC=("flu_vs_rsv_logFC", "median"),
    )
    .reset_index()
)

detail["leading_edge_fraction"] = (
    detail["leading_edge_genes"] /
    detail["n_genes"]
)

detail["rsv_replicated_fraction"] = (
    detail["rsv_replicated_genes"] /
    detail["n_genes"]
)

detail.to_csv(
    TABLES /
    "DS001_GSE38900_pathogen_response_pattern_summary.tsv",
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# Broad category summary
# ---------------------------------------------------------

broad = (
    x.loc[x["direct_significant"]]
    .groupby("pathogen_response_broad_class", dropna=False)
    .agg(
        n_genes=("gene_symbol", "size"),
        leading_edge_genes=("session17_leading_edge", "sum"),
        rsv_replicated_genes=("rsv_replicated", "sum"),
        shared_replicated_genes=("shared_replicated", "sum"),
    )
    .reset_index()
)

broad["leading_edge_fraction"] = (
    broad["leading_edge_genes"] /
    broad["n_genes"]
)

broad.to_csv(
    TABLES /
    "DS001_GSE38900_pathogen_response_broad_summary.tsv",
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# Direction matrix
# ---------------------------------------------------------

matrix = pd.crosstab(
    x.loc[x["direct_significant"], "flu_state"],
    x.loc[x["direct_significant"], "rsv_state"],
    margins=True
)

matrix.to_csv(
    TABLES /
    "DS001_GSE38900_direct_DE_infection_state_matrix.tsv",
    sep="\t"
)


# ---------------------------------------------------------
# Export each interpretable pattern
# ---------------------------------------------------------

for pattern in sorted(
    x.loc[x["direct_significant"],
          "pathogen_response_pattern"].unique()
):

    z = x.loc[
        x["pathogen_response_pattern"] == pattern
    ].copy()

    safe = pattern.replace("/", "_").replace(" ", "_")

    z.to_csv(
        TABLES / f"pattern_{safe}.tsv",
        sep="\t",
        index=False
    )


# ---------------------------------------------------------
# Print
# ---------------------------------------------------------

print("UPDATED MASTER TABLE")
print(updated)

print("\nDETAILED RESPONSE PATTERNS")
print(
    detail.sort_values(
        "n_genes",
        ascending=False
    ).to_string(index=False)
)

print("\nBROAD RESPONSE CLASSES")
print(
    broad.sort_values(
        "n_genes",
        ascending=False
    ).to_string(index=False)
)

print("\nINFECTION STATE MATRIX")
print(matrix.to_string())
