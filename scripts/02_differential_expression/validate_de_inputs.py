#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from datetime import datetime

ROOT = Path.home() / "Influenza_RSV_Project"

EXPR_FILES = {
    "GPL10558": ROOT / "data/processed/DS001_GSE38900/analysis_ready/GPL10558_gene_expression.tsv.gz",
    "GPL6884": ROOT / "data/processed/DS001_GSE38900/analysis_ready/GPL6884_gene_expression.tsv.gz",
}

METADATA_FILE = (
    ROOT
    / "data/processed/DS001_GSE38900"
    / "DS001_GSE38900_master_sample_metadata_harmonized.tsv"
)

OUTPUT_DIR = (
    ROOT
    / "results/DS001_GSE38900/differential_expression/tables"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("SESSION 16 — DIFFERENTIAL EXPRESSION INPUT VALIDATION")
print("=" * 80)
print("Timestamp:", datetime.now().isoformat())
print()

# ------------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------------

metadata = pd.read_csv(METADATA_FILE, sep="\t")

print("METADATA")
print("-" * 80)
print("File:", METADATA_FILE)
print("Rows:", metadata.shape[0])
print("Columns:", metadata.shape[1])

print("\nMetadata columns:")
for col in metadata.columns:
    print(" ", repr(col))

print("\nFirst five metadata rows:")
print(metadata.head().to_string(index=False))

# ------------------------------------------------------------------
# Validate expression matrices
# ------------------------------------------------------------------

summary = []

for platform, path in EXPR_FILES.items():

    print("\n" + "=" * 80)
    print(platform)
    print("=" * 80)

    df = pd.read_csv(path, sep="\t", compression="gzip")

    annotation_cols = [
        "gene_symbol",
        "entrez_gene_id",
        "probe_id",
    ]

    missing_annotation_cols = [
        col for col in annotation_cols if col not in df.columns
    ]

    if missing_annotation_cols:
        raise ValueError(
            f"{platform}: missing expected annotation columns: "
            f"{missing_annotation_cols}"
        )

    sample_cols = [
        col for col in df.columns
        if col not in annotation_cols
    ]

    expression = df[sample_cols]

    missing_expression = int(expression.isna().sum().sum())
    missing_gene_symbol = int(df["gene_symbol"].isna().sum())
    missing_entrez = int(df["entrez_gene_id"].isna().sum())
    missing_probe = int(df["probe_id"].isna().sum())

    duplicate_gene_symbols = int(
        df["gene_symbol"].duplicated().sum()
    )

    duplicate_sample_ids = (
        len(sample_cols) - len(set(sample_cols))
    )

    print("Genes:", len(df))
    print("Samples:", len(sample_cols))
    print("Expression missing values:", missing_expression)
    print("Missing gene symbols:", missing_gene_symbol)
    print("Missing Entrez IDs:", missing_entrez)
    print("Missing probe IDs:", missing_probe)
    print("Duplicate gene symbols:", duplicate_gene_symbols)
    print("Duplicate sample IDs:", duplicate_sample_ids)

    print("\nSample IDs:")
    for sample in sample_cols:
        print(" ", sample)

    summary.append({
        "platform": platform,
        "gene_count": len(df),
        "sample_count": len(sample_cols),
        "missing_expression_values": missing_expression,
        "missing_gene_symbols": missing_gene_symbol,
        "missing_entrez_ids": missing_entrez,
        "missing_probe_ids": missing_probe,
        "duplicate_gene_symbols": duplicate_gene_symbols,
        "duplicate_sample_ids": duplicate_sample_ids,
    })

summary_df = pd.DataFrame(summary)

summary_file = (
    OUTPUT_DIR
    / "DS001_GSE38900_DE_input_validation_summary.tsv"
)

summary_df.to_csv(
    summary_file,
    sep="\t",
    index=False,
)

print("\n" + "=" * 80)
print("VALIDATION SUMMARY")
print("=" * 80)
print(summary_df.to_string(index=False))

print("\nSaved:")
print(summary_file)

print("\nSESSION 16 INPUT VALIDATION COMPLETE")
