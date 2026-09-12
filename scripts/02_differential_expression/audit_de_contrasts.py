#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

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
print("SESSION 16 — DE CONTRAST AND CONFOUNDER AUDIT")
print("=" * 80)

# ------------------------------------------------------------
# Group × platform × sex
# ------------------------------------------------------------

print("\nGROUP × PLATFORM × SEX")
print("-" * 80)

sex_table = pd.crosstab(
    [df["platform_id"], df["harmonized_group"]],
    df["sex_standardized"],
    margins=True
)

print(sex_table.to_string())

sex_table.to_csv(
    OUTDIR / "DS001_GSE38900_group_platform_sex.tsv",
    sep="\t"
)

# ------------------------------------------------------------
# Age distribution for GPL6884
# ------------------------------------------------------------

g = df[df["platform_id"] == "GPL6884"].copy()

g["age_months_harmonized"] = pd.to_numeric(
    g["age_months_harmonized"],
    errors="coerce"
)

print("\n" + "=" * 80)
print("GPL6884 AGE DISTRIBUTION BY GROUP")
print("=" * 80)

age_summary = (
    g.groupby("harmonized_group")["age_months_harmonized"]
    .agg(
        n="count",
        mean="mean",
        std="std",
        median="median",
        minimum="min",
        maximum="max"
    )
    .sort_index()
)

print(age_summary.to_string())

age_summary.to_csv(
    OUTDIR / "GPL6884_age_distribution_by_group.tsv",
    sep="\t"
)

# ------------------------------------------------------------
# Age distribution by pathogen
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("GPL6884 AGE DISTRIBUTION BY PATHOGEN")
print("=" * 80)

age_pathogen = (
    g.groupby("pathogen")["age_months_harmonized"]
    .agg(
        n="count",
        mean="mean",
        std="std",
        median="median",
        minimum="min",
        maximum="max"
    )
)

print(age_pathogen.to_string())

# ------------------------------------------------------------
# Detailed candidate contrasts
# ------------------------------------------------------------

candidate_contrasts = {
    "GPL10558_RSV_acute_vs_control": (
        (df["platform_id"] == "GPL10558")
        & (df["harmonized_group"].isin(
            ["RSV acute", "Healthy control"]
        ))
    ),

    "GPL6884_RSV_acute_vs_control": (
        (df["platform_id"] == "GPL6884")
        & (df["harmonized_group"].isin(
            ["RSV acute", "Healthy control"]
        ))
    ),

    "GPL6884_InfluenzaA_vs_control": (
        (df["platform_id"] == "GPL6884")
        & (df["harmonized_group"].isin(
            ["Influenza A acute", "Healthy control"]
        ))
    ),

    "GPL6884_InfluenzaA_vs_RSVacute": (
        (df["platform_id"] == "GPL6884")
        & (df["harmonized_group"].isin(
            ["Influenza A acute", "RSV acute"]
        ))
    ),

    "GPL6884_RSVacute_vs_recovery": (
        (df["platform_id"] == "GPL6884")
        & (df["harmonized_group"].isin(
            ["RSV acute", "RSV recovery"]
        ))
    ),
}

rows = []

print("\n" + "=" * 80)
print("CANDIDATE CONTRASTS")
print("=" * 80)

for name, mask in candidate_contrasts.items():

    sub = df.loc[mask].copy()

    print("\n", name)
    print("-" * 80)

    print(
        sub["harmonized_group"]
        .value_counts()
        .to_string()
    )

    print("\nSex:")
    print(
        sub["sex_standardized"]
        .value_counts()
        .to_string()
    )

    if sub["age_months_harmonized"].notna().any():

        age = pd.to_numeric(
            sub["age_months_harmonized"],
            errors="coerce"
        )

        print("\nAge months:")
        print(age.describe().to_string())

    for group, group_df in sub.groupby("harmonized_group"):

        age = pd.to_numeric(
            group_df["age_months_harmonized"],
            errors="coerce"
        )

        rows.append({
            "contrast": name,
            "group": group,
            "n": len(group_df),
            "male_n": int(
                (group_df["sex_standardized"] == "Male").sum()
            ),
            "female_n": int(
                (group_df["sex_standardized"] == "Female").sum()
            ),
            "age_n": int(age.notna().sum()),
            "age_mean_months": age.mean(),
            "age_median_months": age.median(),
            "age_min_months": age.min(),
            "age_max_months": age.max(),
        })


contrast_summary = pd.DataFrame(rows)

contrast_summary.to_csv(
    OUTDIR / "DS001_GSE38900_candidate_contrast_summary.tsv",
    sep="\t",
    index=False
)

print("\n" + "=" * 80)
print("CONTRAST SUMMARY")
print("=" * 80)

print(contrast_summary.to_string(index=False))

print("\nSaved:")
print(
    OUTDIR /
    "DS001_GSE38900_candidate_contrast_summary.tsv"
)

print("\nSESSION 16 CONTRAST AUDIT COMPLETE")
