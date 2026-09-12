#!/usr/bin/env python3

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path.home() / "Influenza_RSV_Project"
INPUT_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "DS001_GSE38900"
    / "DS001_GSE38900_master_sample_metadata_initial.tsv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "DS001_GSE38900"
)
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "logs"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_METADATA = OUTPUT_DIR / "DS001_GSE38900_master_sample_metadata_harmonized.tsv"
OUTPUT_QC = TABLE_DIR / "DS001_GSE38900_phenotype_harmonization_qc.tsv"
OUTPUT_GROUP_COUNTS = TABLE_DIR / "DS001_GSE38900_harmonized_group_counts.tsv"
OUTPUT_ELIGIBILITY = TABLE_DIR / "DS001_GSE38900_analysis_eligibility_counts.tsv"
OUTPUT_SUMMARY = LOG_DIR / "DS001_GSE38900_phenotype_harmonization_summary.json"


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def normalize_sex(value: object) -> str:
    text = clean_text(value).lower()
    mapping = {
        "male": "Male",
        "m": "Male",
        "female": "Female",
        "f": "Female",
    }
    return mapping.get(text, "Unknown" if text == "" else text.title())


def normalize_tissue(value: object) -> str:
    text = clean_text(value).lower()
    if text == "whole blood":
        return "Whole blood"
    if text == "":
        return "Unknown"
    return clean_text(value)


def classify_group(row: pd.Series) -> tuple[str, str, str, bool]:
    combined = " | ".join(
        clean_text(row.get(column, ""))
        for column in ["sample_group", "diagnosis", "title"]
    ).lower()

    if "healthy matched control" in combined or "healthy control" in combined:
        return "Healthy control", "Control", "Control", True

    if "1-2 month after" in combined and "rsv" in combined:
        return "RSV recovery", "RSV", "Recovery", False

    if "influenza a" in combined:
        return "Influenza A acute", "Influenza A", "Acute", False

    if "rhinovirus" in combined or "hrv" in combined:
        return "HRV acute", "HRV", "Acute", False

    if "respiratory syncytial virus" in combined or "rsv" in combined:
        return "RSV acute", "RSV", "Acute", False

    return "Unclassified", "Unknown", "Unknown", False


def harmonize_age(row: pd.Series) -> tuple[float, str, str]:
    platform = clean_text(row.get("platform_id", ""))

    if platform == "GPL6884":
        value = pd.to_numeric(row.get("age_months"), errors="coerce")
        if pd.notna(value):
            return float(value), "months", "age_months"
        return np.nan, "unknown", "age_months_missing"

    if platform == "GPL10558":
        value = pd.to_numeric(row.get("age"), errors="coerce")
        if pd.notna(value):
            # Unit intentionally not assumed. Preserve value in a separate column
            # and do not use it as age_months_harmonized until unit verification.
            return np.nan, "unverified", "age_unit_unverified"
        return np.nan, "unknown", "age_missing"

    return np.nan, "unknown", "platform_unknown"


def main() -> None:
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Missing input metadata: {INPUT_PATH}")

    df = pd.read_csv(INPUT_PATH, sep="\t", dtype=str)

    # Preserve all source fields and add canonical fields only.
    group_info = df.apply(classify_group, axis=1, result_type="expand")
    group_info.columns = [
        "harmonized_group",
        "pathogen",
        "infection_phase",
        "is_control",
    ]
    df = pd.concat([df, group_info], axis=1)

    df["sex_standardized"] = df["gender"].apply(normalize_sex)
    df["tissue_standardized"] = df["tissue"].apply(normalize_tissue)

    age_info = df.apply(harmonize_age, axis=1, result_type="expand")
    age_info.columns = [
        "age_months_harmonized",
        "age_unit_status",
        "age_harmonization_note",
    ]
    df = pd.concat([df, age_info], axis=1)

    # Additional preserved raw-age field for GPL10558.
    df["age_raw_gpl10558"] = np.where(
        df["platform_id"].eq("GPL10558"),
        pd.to_numeric(df["age"], errors="coerce"),
        np.nan,
    )

    # Canonical analysis flags.
    df["eligible_primary_infection_vs_control"] = (
        df["harmonized_group"].isin(
            ["Healthy control", "RSV acute", "Influenza A acute"]
        )
        & df["tissue_standardized"].eq("Whole blood")
    )

    df["eligible_direct_influenza_vs_rsv"] = (
        df["harmonized_group"].isin(["RSV acute", "Influenza A acute"])
        & df["platform_id"].eq("GPL6884")
        & df["tissue_standardized"].eq("Whole blood")
    )

    df["eligible_rsv_acute_vs_control"] = (
        df["harmonized_group"].isin(["RSV acute", "Healthy control"])
        & df["tissue_standardized"].eq("Whole blood")
    )

    df["eligible_rsv_recovery_analysis"] = (
        df["harmonized_group"].isin(["RSV acute", "RSV recovery"])
        & df["platform_id"].eq("GPL6884")
        & df["tissue_standardized"].eq("Whole blood")
    )

    # Conservative age-adjustment eligibility:
    # only GPL6884 age_months is treated as verified for now.
    df["eligible_age_adjusted_model"] = (
        df["age_months_harmonized"].notna()
        & df["harmonized_group"].ne("Unclassified")
    )

    # QC flags.
    df["qc_unclassified_group"] = df["harmonized_group"].eq("Unclassified")
    df["qc_unknown_sex"] = df["sex_standardized"].eq("Unknown")
    df["qc_unknown_tissue"] = df["tissue_standardized"].eq("Unknown")
    df["qc_duplicate_geo_accession"] = df.duplicated(
        subset=["geo_accession"],
        keep=False,
    )

    # Stable sample key.
    df.insert(
        0,
        "canonical_sample_id",
        df["dataset_id"].astype(str)
        + "_"
        + df["platform_id"].astype(str)
        + "_"
        + df["geo_accession"].astype(str),
    )

    # Save harmonized metadata.
    df.to_csv(OUTPUT_METADATA, sep="\t", index=False)

    # Group counts.
    group_counts = (
        df.groupby(
            ["platform_id", "harmonized_group", "pathogen", "infection_phase"],
            dropna=False,
        )
        .size()
        .reset_index(name="n_samples")
        .sort_values(["platform_id", "harmonized_group"])
    )
    group_counts.to_csv(OUTPUT_GROUP_COUNTS, sep="\t", index=False)

    # Eligibility counts.
    eligibility_columns = [
        "eligible_primary_infection_vs_control",
        "eligible_direct_influenza_vs_rsv",
        "eligible_rsv_acute_vs_control",
        "eligible_rsv_recovery_analysis",
        "eligible_age_adjusted_model",
    ]
    eligibility_rows = []
    for column in eligibility_columns:
        eligibility_rows.append(
            {
                "analysis_flag": column,
                "eligible_samples": int(df[column].sum()),
                "ineligible_samples": int((~df[column]).sum()),
                "total_samples": int(len(df)),
            }
        )
    pd.DataFrame(eligibility_rows).to_csv(
        OUTPUT_ELIGIBILITY,
        sep="\t",
        index=False,
    )

    # QC table.
    qc_rows = [
        {
            "qc_metric": "total_samples",
            "value": int(len(df)),
        },
        {
            "qc_metric": "unique_geo_accessions",
            "value": int(df["geo_accession"].nunique()),
        },
        {
            "qc_metric": "duplicate_geo_accession_rows",
            "value": int(df["qc_duplicate_geo_accession"].sum()),
        },
        {
            "qc_metric": "unclassified_group_samples",
            "value": int(df["qc_unclassified_group"].sum()),
        },
        {
            "qc_metric": "unknown_sex_samples",
            "value": int(df["qc_unknown_sex"].sum()),
        },
        {
            "qc_metric": "unknown_tissue_samples",
            "value": int(df["qc_unknown_tissue"].sum()),
        },
        {
            "qc_metric": "verified_age_months_samples",
            "value": int(df["age_months_harmonized"].notna().sum()),
        },
        {
            "qc_metric": "gpl10558_age_unit_unverified_samples",
            "value": int(
                df["age_harmonization_note"].eq("age_unit_unverified").sum()
            ),
        },
    ]
    qc_df = pd.DataFrame(qc_rows)
    qc_df.to_csv(OUTPUT_QC, sep="\t", index=False)

    summary = {
        "dataset_id": "DS001",
        "gse_accession": "GSE38900",
        "input_metadata": str(INPUT_PATH),
        "output_metadata": str(OUTPUT_METADATA),
        "total_samples": int(len(df)),
        "group_counts": (
            df["harmonized_group"].value_counts().sort_index().to_dict()
        ),
        "platform_counts": (
            df["platform_id"].value_counts().sort_index().to_dict()
        ),
        "sex_counts": (
            df["sex_standardized"].value_counts().sort_index().to_dict()
        ),
        "unclassified_samples": int(df["qc_unclassified_group"].sum()),
        "duplicate_geo_accession_rows": int(
            df["qc_duplicate_geo_accession"].sum()
        ),
        "age_policy": (
            "GPL6884 age_months accepted as verified months. "
            "GPL10558 age preserved but excluded from age_months_harmonized "
            "until unit verification."
        ),
        "outputs": {
            "harmonized_metadata": str(OUTPUT_METADATA),
            "group_counts": str(OUTPUT_GROUP_COUNTS),
            "eligibility_counts": str(OUTPUT_ELIGIBILITY),
            "qc_table": str(OUTPUT_QC),
        },
    }

    OUTPUT_SUMMARY.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print("Canonical phenotype harmonization complete.")
    print(f"Total samples: {len(df)}")
    print(f"Unclassified samples: {int(df['qc_unclassified_group'].sum())}")
    print(
        "Duplicate GEO accession rows: "
        f"{int(df['qc_duplicate_geo_accession'].sum())}"
    )
    print("\nHarmonized group counts:")
    print(df["harmonized_group"].value_counts().to_string())
    print("\nEligibility counts:")
    for column in eligibility_columns:
        print(f"  {column}: {int(df[column].sum())}")
    print(f"\nHarmonized metadata: {OUTPUT_METADATA}")
    print(f"QC table: {OUTPUT_QC}")
    print(f"Group counts: {OUTPUT_GROUP_COUNTS}")
    print(f"Eligibility table: {OUTPUT_ELIGIBILITY}")


if __name__ == "__main__":
    main()
