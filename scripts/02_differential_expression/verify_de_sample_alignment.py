#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

ROOT = Path.home() / "Influenza_RSV_Project"

META_FILE = (
    ROOT /
    "data/processed/DS001_GSE38900/"
    "DS001_GSE38900_master_sample_metadata_harmonized.tsv"
)

EXPR_FILES = {
    "GPL10558":
        ROOT /
        "data/processed/DS001_GSE38900/analysis_ready/"
        "GPL10558_gene_expression.tsv.gz",

    "GPL6884":
        ROOT /
        "data/processed/DS001_GSE38900/analysis_ready/"
        "GPL6884_gene_expression.tsv.gz",
}

OUTDIR = (
    ROOT /
    "results/DS001_GSE38900/"
    "differential_expression/tables"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(META_FILE, sep="\t")

annotation_cols = [
    "gene_symbol",
    "entrez_gene_id",
    "probe_id",
]

all_mappings = []

print("=" * 80)
print("SESSION 16 — EXPRESSION / METADATA SAMPLE ALIGNMENT")
print("=" * 80)

overall_pass = True

for platform, expr_file in EXPR_FILES.items():

    print("\n" + "=" * 80)
    print(platform)
    print("=" * 80)

    expr = pd.read_csv(
        expr_file,
        sep="\t",
        compression="gzip"
    )

    sample_cols = [
        c for c in expr.columns
        if c not in annotation_cols
    ]

    m = meta[
        meta["platform_id"] == platform
    ].copy()

    # Expression sample identifiers correspond to
    # the GEO metadata description field in this dataset.
    metadata_ids = (
        m["description"]
        .astype(str)
        .str.strip()
        .tolist()
    )

    expression_ids = [
        str(x).strip()
        for x in sample_cols
    ]

    expr_set = set(expression_ids)
    meta_set = set(metadata_ids)

    missing_in_metadata = sorted(
        expr_set - meta_set
    )

    missing_in_expression = sorted(
        meta_set - expr_set
    )

    duplicate_metadata_ids = (
        m["description"]
        .astype(str)
        .str.strip()
        .duplicated()
        .sum()
    )

    duplicate_expression_ids = (
        len(expression_ids)
        - len(set(expression_ids))
    )

    print("Expression samples:", len(expression_ids))
    print("Metadata samples:", len(metadata_ids))

    print(
        "Duplicate expression IDs:",
        duplicate_expression_ids
    )

    print(
        "Duplicate metadata description IDs:",
        duplicate_metadata_ids
    )

    print(
        "Expression IDs missing from metadata:",
        len(missing_in_metadata)
    )

    if missing_in_metadata:
        for x in missing_in_metadata:
            print("  ", repr(x))

    print(
        "Metadata IDs missing from expression:",
        len(missing_in_expression)
    )

    if missing_in_expression:
        for x in missing_in_expression:
            print("  ", repr(x))

    platform_pass = (
        len(expression_ids) == len(metadata_ids)
        and duplicate_expression_ids == 0
        and duplicate_metadata_ids == 0
        and len(missing_in_metadata) == 0
        and len(missing_in_expression) == 0
    )

    print(
        "\nALIGNMENT STATUS:",
        "PASS" if platform_pass else "FAIL"
    )

    if not platform_pass:
        overall_pass = False

    # Build explicit mapping in expression-column order
    lookup = (
        m.assign(
            expression_sample_id=
            m["description"]
            .astype(str)
            .str.strip()
        )
        .set_index("expression_sample_id")
    )

    if platform_pass:

        for position, sample_id in enumerate(
            expression_ids,
            start=1
        ):

            row = lookup.loc[sample_id]

            all_mappings.append({
                "platform_id": platform,
                "expression_column_position": position,
                "expression_sample_id": sample_id,
                "canonical_sample_id":
                    row["canonical_sample_id"],
                "geo_accession":
                    row["geo_accession"],
                "harmonized_group":
                    row["harmonized_group"],
                "pathogen":
                    row["pathogen"],
                "infection_phase":
                    row["infection_phase"],
                "sex_standardized":
                    row["sex_standardized"],
                "age_months_harmonized":
                    row["age_months_harmonized"],
            })


mapping_df = pd.DataFrame(all_mappings)

mapping_file = (
    OUTDIR /
    "DS001_GSE38900_expression_metadata_sample_mapping.tsv"
)

mapping_df.to_csv(
    mapping_file,
    sep="\t",
    index=False
)

print("\n" + "=" * 80)
print("FINAL ALIGNMENT RESULT")
print("=" * 80)

print(
    "OVERALL STATUS:",
    "PASS" if overall_pass else "FAIL"
)

print("Mapped samples:", len(mapping_df))
print("Mapping file:")
print(mapping_file)

if not overall_pass:
    raise SystemExit(
        "\nERROR: Expression/metadata alignment failed. "
        "Differential expression analysis must not proceed."
    )

print(
    "\nSESSION 16 SAMPLE ALIGNMENT COMPLETE"
)
