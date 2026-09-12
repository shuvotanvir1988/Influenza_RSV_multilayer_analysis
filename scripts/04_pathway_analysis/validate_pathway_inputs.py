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
    "pathway_analysis/tables"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

FILES = {

    "GPL10558_RSV_vs_control":
        DE_DIR /
        "GPL10558_RSV_acute_vs_control_full_DE.tsv",

    "GPL6884_RSV_vs_control":
        DE_DIR /
        "GPL6884_RSV_acute_vs_control_full_DE.tsv",

    "GPL6884_Influenza_vs_control":
        DE_DIR /
        "GPL6884_InfluenzaA_vs_control_full_DE.tsv",

    "GPL6884_Influenza_vs_RSV":
        DE_DIR /
        "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv",

    "RSV_cross_platform_replication":
        DE_DIR /
        "DS001_GSE38900_RSV_cross_platform_replication.tsv",
}

print("=" * 85)
print("SESSION 17 — PATHWAY ANALYSIS INPUT VALIDATION")
print("=" * 85)

rows = []

for name, path in FILES.items():

    print("\n" + "=" * 85)
    print(name)
    print("=" * 85)

    print("File:")
    print(path)

    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(
        path,
        sep="\t"
    )

    print("\nShape:", df.shape)

    print("\nColumns:")
    for col in df.columns:
        print(" ", col)

    if "gene_symbol" not in df.columns:
        raise ValueError(
            f"{name}: gene_symbol missing"
        )

    gene = (
        df["gene_symbol"]
        .astype("string")
        .str.strip()
    )

    missing_gene = int(
        gene.isna().sum()
        +
        (gene == "").sum()
    )

    duplicated_gene = int(
        gene.dropna().duplicated().sum()
    )

    print("\nGene symbols:")
    print("  Total rows:", len(df))
    print("  Missing:", missing_gene)
    print("  Duplicated:", duplicated_gene)
    print(
        "  Unique:",
        gene.dropna().nunique()
    )

    # --------------------------------------------
    # DE tables
    # --------------------------------------------

    if "logFC" in df.columns:

        required = [
            "logFC",
            "t",
            "P.Value",
            "adj.P.Val",
        ]

        missing_required = [
            c for c in required
            if c not in df.columns
        ]

        if missing_required:
            raise ValueError(
                f"{name}: missing "
                f"{missing_required}"
            )

        print("\nDE statistics:")

        for col in required:

            numeric = pd.to_numeric(
                df[col],
                errors="coerce"
            )

            print(
                f"  {col}: "
                f"nonmissing={numeric.notna().sum()}, "
                f"missing={numeric.isna().sum()}, "
                f"min={numeric.min()}, "
                f"max={numeric.max()}"
            )

        finite_t = pd.to_numeric(
            df["t"],
            errors="coerce"
        ).notna().sum()

        print(
            "  Genes with usable t statistic:",
            finite_t
        )

    rows.append({
        "dataset": name,
        "rows": len(df),
        "unique_gene_symbols":
            gene.dropna().nunique(),
        "missing_gene_symbols":
            missing_gene,
        "duplicate_gene_symbols":
            duplicated_gene,
        "has_logFC":
            "logFC" in df.columns,
        "has_t":
            "t" in df.columns,
        "has_p_value":
            "P.Value" in df.columns,
        "has_FDR":
            "adj.P.Val" in df.columns,
    })


summary = pd.DataFrame(rows)

outfile = (
    OUTDIR /
    "DS001_GSE38900_pathway_input_validation.tsv"
)

summary.to_csv(
    outfile,
    sep="\t",
    index=False
)

print("\n" + "=" * 85)
print("PATHWAY INPUT SUMMARY")
print("=" * 85)

print(
    summary.to_string(
        index=False
    )
)

print("\nSaved:")
print(outfile)

print(
    "\nSESSION 17 PATHWAY INPUT "
    "VALIDATION COMPLETE"
)
