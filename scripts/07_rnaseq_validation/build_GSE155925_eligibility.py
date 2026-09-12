import pandas as pd
import re
from pathlib import Path

INFILE = Path(
    "data/rnaseq_validation/GSE155925/metadata/"
    "GSE155925_sample_metadata_raw.tsv"
)

OUTDIR = Path("results/rnaseq_validation/tables")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INFILE, sep="\t")

def extract(pattern, text):
    m = re.search(pattern, str(text), flags=re.I)
    return m.group(1).strip() if m else None

df["sex"] = df["characteristics"].apply(
    lambda x: extract(r"Sex:\s*([^|]+)", x)
)

df["age_months"] = pd.to_numeric(
    df["characteristics"].apply(
        lambda x: extract(r"age \(months\):\s*([^|]+)", x)
    ),
    errors="coerce"
)

df["pathogen"] = df["characteristics"].apply(
    lambda x: extract(r"pathogen:\s*([^|]+)", x)
)

df["hospital_batch"] = pd.to_numeric(
    df["characteristics"].apply(
        lambda x: extract(r"batch \(hospital\):\s*([^|]+)", x)
    ),
    errors="coerce"
)

df["enrollment_year_batch"] = pd.to_numeric(
    df["characteristics"].apply(
        lambda x: extract(r"batch \(year enrollment\):\s*([^|]+)", x)
    ),
    errors="coerce"
)

# ------------------------------------------------------------
# Frozen primary eligibility definition
#
# RSV_ONLY:
#   pathogen field exactly "RSV"
#
# VIRUS_NEGATIVE:
#   pathogen field exactly "negative"
#
# All coinfections and other pathogens are excluded from the
# primary RSV-vs-negative comparison.
# ------------------------------------------------------------

path_upper = df["pathogen"].fillna("").str.strip().str.upper()

df["primary_group"] = "EXCLUDED_OTHER"

df.loc[path_upper.eq("RSV"), "primary_group"] = "RSV_ONLY"
df.loc[path_upper.eq("NEGATIVE"), "primary_group"] = "VIRUS_NEGATIVE"

df["primary_eligible"] = df["primary_group"].isin(
    ["RSV_ONLY", "VIRUS_NEGATIVE"]
)

def exclusion_reason(row):
    if row["primary_eligible"]:
        return ""

    p = str(row["pathogen"])

    if "," in p and "RSV" in p.upper():
        return "RSV_COINFECTION"

    return "OTHER_PATHOGEN"

df["exclusion_reason"] = df.apply(exclusion_reason, axis=1)

# Preserve only explicit, auditable fields
cols = [
    "sample_index",
    "count_matrix_label",
    "geo_accession",
    "description",
    "sex",
    "age_months",
    "pathogen",
    "hospital_batch",
    "enrollment_year_batch",
    "primary_group",
    "primary_eligible",
    "exclusion_reason",
]

elig = df[cols].copy()

outfile = OUTDIR / "GSE155925_sample_eligibility_FROZEN_v1.0.tsv"
elig.to_csv(outfile, sep="\t", index=False)

print("Written:", outfile)

print("\n=== PATHOGEN COUNTS ===")
print(elig["pathogen"].value_counts(dropna=False).sort_index())

print("\n=== PRIMARY GROUP COUNTS ===")
print(elig["primary_group"].value_counts())

print("\n=== EXCLUSION REASONS ===")
print(
    elig.loc[~elig["primary_eligible"], "exclusion_reason"]
        .value_counts()
)

print("\n=== PRIMARY COHORT ===")
primary = elig[elig["primary_eligible"]].copy()

print(primary.groupby("primary_group").agg(
    n=("geo_accession", "size"),
    age_median=("age_months", "median"),
    age_min=("age_months", "min"),
    age_max=("age_months", "max"),
))

print("\n=== SEX BY PRIMARY GROUP ===")
print(pd.crosstab(primary["primary_group"], primary["sex"]))

print("\n=== HOSPITAL BATCH ===")
print(pd.crosstab(
    primary["primary_group"],
    primary["hospital_batch"]
))

print("\n=== ENROLLMENT-YEAR BATCH ===")
print(pd.crosstab(
    primary["primary_group"],
    primary["enrollment_year_batch"]
))

print("\n=== ELIGIBLE SAMPLE LIST ===")
print(
    primary[
        [
            "count_matrix_label",
            "geo_accession",
            "primary_group",
            "age_months",
            "sex",
            "hospital_batch",
            "enrollment_year_batch",
        ]
    ].to_string(index=False)
)
