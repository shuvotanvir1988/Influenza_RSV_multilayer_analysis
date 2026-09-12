#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from datetime import datetime

ROOT = Path.home() / "Influenza_RSV_Project"

META = (
    ROOT / "data/processed/DS001_GSE38900/"
    "DS001_GSE38900_master_sample_metadata_harmonized.tsv"
)

OUTDIR = (
    ROOT / "results/DS001_GSE38900/"
    "differential_expression/tables"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(META, sep="\t")

print("=" * 80)
print("SESSION 16 — DE SAMPLE INVENTORY")
print("=" * 80)
print("Timestamp:", datetime.now().isoformat())
print("Total metadata rows:", len(df))


def show_counts(column):
    print("\n" + "-" * 80)
    print(column)
    print("-" * 80)

    if column not in df.columns:
        print("COLUMN NOT FOUND")
        return

    print(
        df[column]
        .fillna("<MISSING>")
        .value_counts(dropna=False)
        .to_string()
    )


# ------------------------------------------------------------
# Core biological variables
# ------------------------------------------------------------

for column in [
    "platform_id",
    "harmonized_group",
    "pathogen",
    "infection_phase",
    "is_control",
    "sex_standardized",
    "tissue_standardized",
    "age_unit_status",
]:
    show_counts(column)


# ------------------------------------------------------------
# Platform × biological group
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("PLATFORM × HARMONIZED GROUP")
print("=" * 80)

platform_group = pd.crosstab(
    df["platform_id"],
    df["harmonized_group"],
    dropna=False
)

print(platform_group.to_string())

platform_group.to_csv(
    OUTDIR / "DS001_GSE38900_platform_by_harmonized_group.tsv",
    sep="\t"
)


# ------------------------------------------------------------
# Platform × pathogen
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("PLATFORM × PATHOGEN")
print("=" * 80)

platform_pathogen = pd.crosstab(
    df["platform_id"],
    df["pathogen"],
    dropna=False
)

print(platform_pathogen.to_string())

platform_pathogen.to_csv(
    OUTDIR / "DS001_GSE38900_platform_by_pathogen.tsv",
    sep="\t"
)


# ------------------------------------------------------------
# Platform × infection phase
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("PLATFORM × INFECTION PHASE")
print("=" * 80)

platform_phase = pd.crosstab(
    df["platform_id"],
    df["infection_phase"],
    dropna=False
)

print(platform_phase.to_string())

platform_phase.to_csv(
    OUTDIR / "DS001_GSE38900_platform_by_infection_phase.tsv",
    sep="\t"
)


# ------------------------------------------------------------
# Eligibility flags
# ------------------------------------------------------------

eligibility_columns = [
    "eligible_primary_infection_vs_control",
    "eligible_direct_influenza_vs_rsv",
    "eligible_rsv_acute_vs_control",
    "eligible_rsv_recovery_analysis",
    "eligible_age_adjusted_model",
]

print("\n" + "=" * 80)
print("ELIGIBILITY SUMMARY")
print("=" * 80)

eligibility_rows = []

for col in eligibility_columns:

    if col not in df.columns:
        continue

    counts = (
        df.groupby("platform_id")[col]
        .value_counts(dropna=False)
        .unstack(fill_value=0)
    )

    print("\n", col)
    print(counts.to_string())

    for platform in df["platform_id"].dropna().unique():

        sub = df[df["platform_id"] == platform]

        eligible = int(
            sub[col]
            .astype(str)
            .str.lower()
            .eq("true")
            .sum()
        )

        eligibility_rows.append({
            "platform_id": platform,
            "eligibility_flag": col,
            "eligible_samples": eligible,
            "total_platform_samples": len(sub),
        })


eligibility_df = pd.DataFrame(eligibility_rows)

eligibility_df.to_csv(
    OUTDIR / "DS001_GSE38900_DE_eligibility_summary.tsv",
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# QC flags
# ------------------------------------------------------------

qc_columns = [
    "qc_unclassified_group",
    "qc_unknown_sex",
    "qc_unknown_tissue",
    "qc_duplicate_geo_accession",
]

print("\n" + "=" * 80)
print("METADATA QC FLAGS")
print("=" * 80)

for col in qc_columns:
    show_counts(col)


# ------------------------------------------------------------
# Duplicate GEO accessions
# ------------------------------------------------------------

duplicate_geo = df[
    df["geo_accession"].duplicated(keep=False)
].sort_values("geo_accession")

print("\nDuplicate GEO accession rows:",
      len(duplicate_geo))

if len(duplicate_geo) > 0:
    print(
        duplicate_geo[
            [
                "canonical_sample_id",
                "geo_accession",
                "platform_id",
                "harmonized_group",
            ]
        ].to_string(index=False)
    )


# ------------------------------------------------------------
# Age audit
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("AGE AUDIT")
print("=" * 80)

for platform, sub in df.groupby("platform_id"):

    print("\nPlatform:", platform)
    print("Samples:", len(sub))

    if "age" in sub.columns:
        print(
            "Raw age non-missing:",
            int(sub["age"].notna().sum())
        )

    if "age_months_harmonized" in sub.columns:
        print(
            "Harmonized age non-missing:",
            int(sub["age_months_harmonized"].notna().sum())
        )

    if "eligible_age_adjusted_model" in sub.columns:
        print(
            "Age-adjusted eligible:",
            int(
                sub["eligible_age_adjusted_model"]
                .astype(str)
                .str.lower()
                .eq("true")
                .sum()
            )
        )


# ------------------------------------------------------------
# Save compact sample inventory
# ------------------------------------------------------------

inventory_columns = [
    "canonical_sample_id",
    "geo_accession",
    "platform_id",
    "title",
    "harmonized_group",
    "pathogen",
    "infection_phase",
    "is_control",
    "sex_standardized",
    "tissue_standardized",
    "age",
    "age_months_harmonized",
    "age_unit_status",
    "eligible_primary_infection_vs_control",
    "eligible_direct_influenza_vs_rsv",
    "eligible_rsv_acute_vs_control",
    "eligible_rsv_recovery_analysis",
    "eligible_age_adjusted_model",
]

inventory_columns = [
    x for x in inventory_columns
    if x in df.columns
]

inventory = df[inventory_columns].copy()

inventory.to_csv(
    OUTDIR / "DS001_GSE38900_DE_sample_inventory.tsv",
    sep="\t",
    index=False
)


print("\n" + "=" * 80)
print("FILES SAVED")
print("=" * 80)

for filename in [
    "DS001_GSE38900_DE_sample_inventory.tsv",
    "DS001_GSE38900_DE_eligibility_summary.tsv",
    "DS001_GSE38900_platform_by_harmonized_group.tsv",
    "DS001_GSE38900_platform_by_pathogen.tsv",
    "DS001_GSE38900_platform_by_infection_phase.tsv",
]:
    print(OUTDIR / filename)

print("\nSESSION 16 SAMPLE INVENTORY COMPLETE")
