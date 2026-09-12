#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path.home() / "Influenza_RSV_Project"

TABLE_DIR = (
    ROOT /
    "results/DS001_GSE38900/"
    "differential_expression/tables"
)

FILES = {
    "GPL10558_RSV_vs_control":
        TABLE_DIR /
        "GPL10558_RSV_acute_vs_control_full_DE.tsv",

    "GPL6884_RSV_vs_control":
        TABLE_DIR /
        "GPL6884_RSV_acute_vs_control_full_DE.tsv",

    "GPL6884_Influenza_vs_control":
        TABLE_DIR /
        "GPL6884_InfluenzaA_vs_control_full_DE.tsv",

    "GPL6884_Influenza_vs_RSV":
        TABLE_DIR /
        "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv",
}

print("=" * 90)
print("SESSION 16 — PRIMARY DE RESULT AUDIT")
print("=" * 90)

summary_rows = []

for contrast, path in FILES.items():

    print("\n" + "=" * 90)
    print(contrast)
    print("=" * 90)

    df = pd.read_csv(path, sep="\t")

    required = [
        "gene_symbol",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B",
    ]

    missing = [
        x for x in required
        if x not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{contrast}: missing columns {missing}"
        )

    print("Genes:", len(df))

    print("\nMissing values:")
    print(
        df[required]
        .isna()
        .sum()
        .to_string()
    )

    print("\nlogFC summary:")
    print(
        df["logFC"]
        .describe(
            percentiles=[
                0.01,
                0.05,
                0.25,
                0.50,
                0.75,
                0.95,
                0.99,
            ]
        )
        .to_string()
    )

    print("\nRaw P-value summary:")
    print(
        df["P.Value"]
        .describe(
            percentiles=[
                0.01,
                0.05,
                0.25,
                0.50,
                0.75,
                0.95,
                0.99,
            ]
        )
        .to_string()
    )

    print("\nAdjusted P-value summary:")
    print(
        df["adj.P.Val"]
        .describe(
            percentiles=[
                0.01,
                0.05,
                0.25,
                0.50,
                0.75,
                0.95,
                0.99,
            ]
        )
        .to_string()
    )

    # --------------------------------------------------------
    # Significance thresholds
    # --------------------------------------------------------

    thresholds = {}

    for fc in [0, 0.25, 0.5, 1.0]:

        mask = (
            (df["adj.P.Val"] < 0.05)
            & (df["logFC"].abs() >= fc)
        )

        thresholds[fc] = int(mask.sum())

    print("\nFDR < 0.05 counts:")

    for fc, n in thresholds.items():
        print(
            f"  |logFC| >= {fc:0.2f}: {n}"
        )

    # --------------------------------------------------------
    # Direction
    # --------------------------------------------------------

    sig = df[df["adj.P.Val"] < 0.05]

    print("\nDirection among FDR-significant genes:")
    print("  Up:", int((sig["logFC"] > 0).sum()))
    print("  Down:", int((sig["logFC"] < 0).sum()))

    # --------------------------------------------------------
    # Top genes by statistical significance
    # --------------------------------------------------------

    top_p = (
        df.sort_values("P.Value")
        .head(20)[
            [
                "gene_symbol",
                "logFC",
                "AveExpr",
                "t",
                "P.Value",
                "adj.P.Val",
            ]
        ]
    )

    print("\nTOP 20 GENES BY P VALUE")
    print("-" * 90)
    print(
        top_p.to_string(
            index=False
        )
    )

    # --------------------------------------------------------
    # Largest positive effects
    # --------------------------------------------------------

    top_up = (
        df.sort_values(
            "logFC",
            ascending=False
        )
        .head(15)[
            [
                "gene_symbol",
                "logFC",
                "AveExpr",
                "P.Value",
                "adj.P.Val",
            ]
        ]
    )

    print("\nTOP 15 POSITIVE logFC")
    print("-" * 90)
    print(
        top_up.to_string(
            index=False
        )
    )

    # --------------------------------------------------------
    # Largest negative effects
    # --------------------------------------------------------

    top_down = (
        df.sort_values(
            "logFC",
            ascending=True
        )
        .head(15)[
            [
                "gene_symbol",
                "logFC",
                "AveExpr",
                "P.Value",
                "adj.P.Val",
            ]
        ]
    )

    print("\nTOP 15 NEGATIVE logFC")
    print("-" * 90)
    print(
        top_down.to_string(
            index=False
        )
    )

    # --------------------------------------------------------
    # P-value histogram bins
    # --------------------------------------------------------

    bins = np.linspace(
        0,
        1,
        11
    )

    counts, _ = np.histogram(
        df["P.Value"].dropna(),
        bins=bins
    )

    print("\nRAW P-VALUE HISTOGRAM")
    print("-" * 90)

    for i, count in enumerate(counts):

        print(
            f"{bins[i]:.1f}–{bins[i+1]:.1f}: "
            f"{count}"
        )

    # --------------------------------------------------------
    # Compact summary
    # --------------------------------------------------------

    summary_rows.append({
        "contrast": contrast,
        "genes_tested": len(df),

        "median_logFC":
            df["logFC"].median(),

        "median_abs_logFC":
            df["logFC"].abs().median(),

        "max_abs_logFC":
            df["logFC"].abs().max(),

        "median_raw_p":
            df["P.Value"].median(),

        "median_FDR":
            df["adj.P.Val"].median(),

        "FDR_lt_0_05":
            int(
                (df["adj.P.Val"] < 0.05)
                .sum()
            ),

        "FDR_lt_0_05_absFC_0_5":
            thresholds[0.5],

        "FDR_lt_0_05_absFC_1":
            thresholds[1.0],
    })


summary = pd.DataFrame(
    summary_rows
)

summary_file = (
    TABLE_DIR /
    "DS001_GSE38900_primary_DE_diagnostic_summary.tsv"
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)

print("\n" + "=" * 90)
print("MASTER DIAGNOSTIC SUMMARY")
print("=" * 90)

print(
    summary.to_string(
        index=False
    )
)

print("\nSaved:")
print(summary_file)

print("\nSESSION 16 PRIMARY DE AUDIT COMPLETE")
