#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

ROOT = Path.home() / "Influenza_RSV_Project"

DE_DIR = (
    ROOT /
    "results/DS001_GSE38900/"
    "differential_expression/tables"
)

OUTDIR = (
    ROOT /
    "results/DS001_GSE38900/"
    "pathway_analysis/tables/ranked_lists"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

FILES = {
    "GPL10558_RSV_acute_vs_control":
        DE_DIR / "GPL10558_RSV_acute_vs_control_full_DE.tsv",

    "GPL6884_RSV_acute_vs_control":
        DE_DIR / "GPL6884_RSV_acute_vs_control_full_DE.tsv",

    "GPL6884_InfluenzaA_vs_control":
        DE_DIR / "GPL6884_InfluenzaA_vs_control_full_DE.tsv",

    "GPL6884_InfluenzaA_vs_RSVacute":
        DE_DIR / "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv",
}

summary_rows = []

print("=" * 90)
print("SESSION 17 — CREATE RANKED GENE LISTS")
print("PRIMARY RANKING METRIC: MODERATED LIMMA t STATISTIC")
print("=" * 90)

for name, infile in FILES.items():

    print("\n" + "=" * 90)
    print(name)
    print("=" * 90)

    df = pd.read_csv(infile, sep="\t")

    required = [
        "gene_symbol",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val"
    ]

    missing = [x for x in required if x not in df.columns]

    if missing:
        raise ValueError(
            f"{name}: missing required columns {missing}"
        )

    ranked = df[required].copy()

    ranked["gene_symbol"] = (
        ranked["gene_symbol"]
        .astype("string")
        .str.strip()
    )

    for col in ["logFC", "t", "P.Value", "adj.P.Val"]:
        ranked[col] = pd.to_numeric(
            ranked[col],
            errors="coerce"
        )

    # Remove unusable records only.
    ranked = ranked[
        ranked["gene_symbol"].notna() &
        (ranked["gene_symbol"] != "") &
        ranked["t"].notna()
    ].copy()

    # This should already be guaranteed by the gene-level DE tables.
    if ranked["gene_symbol"].duplicated().any():
        raise ValueError(
            f"{name}: duplicated gene symbols detected. "
            "Do not silently collapse them."
        )

    # Highest positive t -> strongest evidence for positive direction.
    # Lowest negative t -> strongest evidence for negative direction.
    ranked = ranked.sort_values(
        by=["t", "gene_symbol"],
        ascending=[False, True]
    ).reset_index(drop=True)

    ranked["rank"] = range(1, len(ranked) + 1)

    ranked = ranked[
        [
            "rank",
            "gene_symbol",
            "t",
            "logFC",
            "P.Value",
            "adj.P.Val"
        ]
    ]

    outfile = OUTDIR / f"{name}_ranked_by_limma_t.tsv"

    ranked.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    # Minimal two-column file suitable for many enrichment tools.
    rnkfile = OUTDIR / f"{name}_limma_t.rnk"

    ranked[
        ["gene_symbol", "t"]
    ].to_csv(
        rnkfile,
        sep="\t",
        index=False,
        header=False
    )

    positive = int((ranked["t"] > 0).sum())
    negative = int((ranked["t"] < 0).sum())
    zero = int((ranked["t"] == 0).sum())

    print("Input genes:", len(df))
    print("Ranked genes:", len(ranked))
    print("Positive t:", positive)
    print("Negative t:", negative)
    print("Zero t:", zero)
    print("Maximum t:", ranked["t"].max())
    print("Minimum t:", ranked["t"].min())

    print("\nTop 5 positive:")
    print(
        ranked.head(5)[
            ["gene_symbol", "t", "logFC", "adj.P.Val"]
        ].to_string(index=False)
    )

    print("\nTop 5 negative:")
    print(
        ranked.tail(5)[
            ["gene_symbol", "t", "logFC", "adj.P.Val"]
        ].to_string(index=False)
    )

    print("\nSaved:")
    print(outfile)
    print(rnkfile)

    summary_rows.append({
        "contrast": name,
        "input_genes": len(df),
        "ranked_genes": len(ranked),
        "positive_t": positive,
        "negative_t": negative,
        "zero_t": zero,
        "maximum_t": ranked["t"].max(),
        "minimum_t": ranked["t"].min(),
        "ranking_metric": "moderated_limma_t"
    })


summary = pd.DataFrame(summary_rows)

summary_file = (
    OUTDIR /
    "DS001_GSE38900_ranked_gene_list_summary.tsv"
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)

print("\n" + "=" * 90)
print("RANKED GENE LIST SUMMARY")
print("=" * 90)

print(summary.to_string(index=False))

print("\nSaved summary:")
print(summary_file)

print("\nSESSION 17 RANKED GENE LIST CREATION COMPLETE")
