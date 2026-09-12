#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("results/DS001_GSE38900/signature_analysis")
TABLES = ROOT / "tables"
SENS = ROOT / "sensitivity"

TABLES.mkdir(parents=True, exist_ok=True)
SENS.mkdir(parents=True, exist_ok=True)

MASTER = (
    TABLES /
    "DS001_GSE38900_gene_evidence_master_with_response_patterns.tsv"
)

x = pd.read_csv(MASTER, sep="\t")

FDR = 0.05


# ---------------------------------------------------------
# 1. Primary selective-response definitions
#
# IMPORTANT:
# These are "selective", not "specific".
# Non-significant does not imply zero biological effect.
# ---------------------------------------------------------

x["influenza_selective_primary"] = (
    x["direct_significant"]
    &
    (
        x["pathogen_response_pattern"].isin([
            "influenza_specific_induction",
            "influenza_specific_suppression"
        ])
    )
)

x["RSV_selective_primary"] = (
    x["direct_significant"]
    &
    (
        x["pathogen_response_pattern"].isin([
            "RSV_specific_induction",
            "RSV_specific_suppression"
        ])
    )
)


# ---------------------------------------------------------
# 2. Opposite-direction class
# ---------------------------------------------------------

x["opposite_direction_primary"] = (
    x["direct_significant"]
    &
    x["pathogen_response_pattern"].isin([
        "opposite_influenza_up_RSV_down",
        "opposite_influenza_down_RSV_up"
    ])
)


# ---------------------------------------------------------
# 3. Shared replicated subclasses
# ---------------------------------------------------------

x["shared_replicated_core"] = (
    x["shared_replicated"]
    &
    ~x["direct_significant"]
)

x["shared_replicated_influenza_amplified"] = (
    x["shared_replicated"]
    &
    x["direct_significant"]
    &
    (x["flu_vs_rsv_logFC"] > 0)
)

x["shared_replicated_RSV_amplified"] = (
    x["shared_replicated"]
    &
    x["direct_significant"]
    &
    (x["flu_vs_rsv_logFC"] < 0)
)


# ---------------------------------------------------------
# 4. Stringent specificity-like sensitivity definitions
#
# Active pathogen must have a meaningful effect.
# Other pathogen must be close to zero.
#
# These are sensitivity analyses, NOT primary definitions.
# ---------------------------------------------------------

settings = [
    # active minimum, inactive maximum
    (0.25, 0.25),
    (0.50, 0.25),
    (0.50, 0.10),
    (1.00, 0.25),
]

summary_rows = []


def add_summary(label, mask, active_col, inactive_col):
    z = x.loc[mask].copy()

    summary_rows.append({
        "definition": label,
        "n_genes": len(z),
        "leading_edge_n":
            int(z["session17_leading_edge"].sum()),
        "leading_edge_fraction":
            float(z["session17_leading_edge"].mean())
            if len(z) else np.nan,
        "median_active_abs_logFC":
            float(z[active_col].abs().median())
            if len(z) else np.nan,
        "median_inactive_abs_logFC":
            float(z[inactive_col].abs().median())
            if len(z) else np.nan,
        "median_direct_abs_logFC":
            float(z["flu_vs_rsv_logFC"].abs().median())
            if len(z) else np.nan,
    })


for active_min, inactive_max in settings:

    # Influenza-selective-like
    flu_mask = (
        x["influenza_selective_primary"]
        &
        (x["flu_logFC"].abs() >= active_min)
        &
        (x["rsv6884_logFC"].abs() <= inactive_max)
    )

    label = (
        f"influenza_stringent_"
        f"activeGE{active_min}_inactiveLE{inactive_max}"
    )

    add_summary(
        label,
        flu_mask,
        "flu_logFC",
        "rsv6884_logFC"
    )

    x[label] = flu_mask

    # RSV-selective-like
    rsv_mask = (
        x["RSV_selective_primary"]
        &
        (x["rsv6884_logFC"].abs() >= active_min)
        &
        (x["flu_logFC"].abs() <= inactive_max)
    )

    label = (
        f"RSV_stringent_"
        f"activeGE{active_min}_inactiveLE{inactive_max}"
    )

    add_summary(
        label,
        rsv_mask,
        "rsv6884_logFC",
        "flu_logFC"
    )

    x[label] = rsv_mask


# ---------------------------------------------------------
# 5. Shared-signature effect-size sensitivity
# ---------------------------------------------------------

shared_thresholds = [0.0, 0.25, 0.50, 1.00]

for threshold in shared_thresholds:

    mask = (
        x["shared_replicated"]
        &
        (x["flu_logFC"].abs() >= threshold)
        &
        (x["rsv6884_logFC"].abs() >= threshold)
        &
        (x["rsv10558_logFC"].abs() >= threshold)
    )

    z = x.loc[mask]

    summary_rows.append({
        "definition":
            f"shared_replicated_all3_abslogFC_GE_{threshold}",
        "n_genes": len(z),
        "leading_edge_n":
            int(z["session17_leading_edge"].sum()),
        "leading_edge_fraction":
            float(z["session17_leading_edge"].mean())
            if len(z) else np.nan,
        "median_active_abs_logFC":
            float(z["flu_logFC"].abs().median())
            if len(z) else np.nan,
        "median_inactive_abs_logFC":
            float(z["rsv6884_logFC"].abs().median())
            if len(z) else np.nan,
        "median_direct_abs_logFC":
            float(z["flu_vs_rsv_logFC"].abs().median())
            if len(z) else np.nan,
    })


# ---------------------------------------------------------
# 6. Direct-difference effect-size sensitivity
# ---------------------------------------------------------

direct_thresholds = [0.0, 0.25, 0.50, 1.00]

for threshold in direct_thresholds:

    mask = (
        x["direct_significant"]
        &
        (x["flu_vs_rsv_logFC"].abs() >= threshold)
    )

    z = x.loc[mask]

    summary_rows.append({
        "definition":
            f"direct_FDR05_abslogFC_GE_{threshold}",
        "n_genes": len(z),
        "leading_edge_n":
            int(z["session17_leading_edge"].sum()),
        "leading_edge_fraction":
            float(z["session17_leading_edge"].mean())
            if len(z) else np.nan,
        "median_active_abs_logFC":
            np.nan,
        "median_inactive_abs_logFC":
            np.nan,
        "median_direct_abs_logFC":
            float(z["flu_vs_rsv_logFC"].abs().median())
            if len(z) else np.nan,
    })


# ---------------------------------------------------------
# 7. Primary class summary
# ---------------------------------------------------------

primary_rows = []

primary_classes = {
    "shared_replicated_core":
        x["shared_replicated_core"],

    "shared_replicated_influenza_amplified":
        x["shared_replicated_influenza_amplified"],

    "shared_replicated_RSV_amplified":
        x["shared_replicated_RSV_amplified"],

    "influenza_selective_primary":
        x["influenza_selective_primary"],

    "RSV_selective_primary":
        x["RSV_selective_primary"],

    "opposite_direction_primary":
        x["opposite_direction_primary"],
}

for name, mask in primary_classes.items():

    z = x.loc[mask]

    primary_rows.append({
        "signature": name,
        "n_genes": len(z),
        "leading_edge_n":
            int(z["session17_leading_edge"].sum()),
        "leading_edge_fraction":
            float(z["session17_leading_edge"].mean())
            if len(z) else np.nan,
        "median_flu_logFC":
            float(z["flu_logFC"].median())
            if len(z) else np.nan,
        "median_RSV_logFC":
            float(z["rsv6884_logFC"].median())
            if len(z) else np.nan,
        "median_direct_logFC":
            float(z["flu_vs_rsv_logFC"].median())
            if len(z) else np.nan,
    })


primary = pd.DataFrame(primary_rows)

primary.to_csv(
    TABLES /
    "DS001_GSE38900_primary_signature_summary.tsv",
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# 8. Sensitivity summary
# ---------------------------------------------------------

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    SENS /
    "DS001_GSE38900_signature_threshold_sensitivity.tsv",
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# 9. Updated evidence table
# ---------------------------------------------------------

x.to_csv(
    TABLES /
    "DS001_GSE38900_gene_evidence_master_final_classification.tsv",
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# 10. Export principal high-confidence signatures
# ---------------------------------------------------------

exports = {
    "DS001_GSE38900_signature_shared_core.tsv":
        x.loc[x["shared_replicated_core"]],

    "DS001_GSE38900_signature_shared_influenza_amplified.tsv":
        x.loc[x["shared_replicated_influenza_amplified"]],

    "DS001_GSE38900_signature_shared_RSV_amplified.tsv":
        x.loc[x["shared_replicated_RSV_amplified"]],

    "DS001_GSE38900_signature_influenza_selective.tsv":
        x.loc[x["influenza_selective_primary"]],

    "DS001_GSE38900_signature_RSV_selective.tsv":
        x.loc[x["RSV_selective_primary"]],

    "DS001_GSE38900_signature_opposite_direction.tsv":
        x.loc[x["opposite_direction_primary"]],
}

for filename, df in exports.items():
    df.to_csv(
        TABLES / filename,
        sep="\t",
        index=False
    )


print("\nPRIMARY SIGNATURE SUMMARY")
print(primary.to_string(index=False))

print("\nTHRESHOLD SENSITIVITY")
print(summary.to_string(index=False))

print("\nEXPORTED SIGNATURES")
for filename, df in exports.items():
    print(f"{filename}: {len(df)} genes")
