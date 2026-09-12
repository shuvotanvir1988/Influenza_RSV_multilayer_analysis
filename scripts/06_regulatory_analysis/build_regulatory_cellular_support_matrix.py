from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path.home() / "Influenza_RSV_Project"
BASE = ROOT / "results/regulatory_driver_analysis"

MASTER = BASE / "scrnaseq_validation/differential_activity/GSE283746_frozen37_celltype_validation_master_v1.0.tsv"
DRIVERS = BASE / "tables/REGULATORY_HIGH_CONFIDENCE_DRIVERS_v1.0.tsv"
RNA = BASE / "rnaseq_validation/GSE155925_frozen37_driver_validation_v1.0.tsv"

OUT = BASE / "tables"
OUT.mkdir(parents=True, exist_ok=True)

m = pd.read_csv(MASTER, sep="\t")
d = pd.read_csv(DRIVERS, sep="\t")
r = pd.read_csv(RNA, sep="\t")

celltypes = [
    "Naive CD4+ T",
    "Naive CD8+ T",
    "TrB",
    "B",
    "NK",
    "Memory CD4+ T",
    "CD14+ Monocyte",
    "Cytotoxic T",
    "Treg",
    "Memory B",
    "CD16+ Monocyte",
]

# Validate complete driver x cell-type grid.
expected = len(d) * len(celltypes)
if len(m) != expected:
    raise ValueError(f"Expected {expected} rows in cell-type master, found {len(m)}")

if set(m["TF"]) != set(d["TF"]):
    raise ValueError("TF set mismatch between frozen drivers and cell-type validation master.")

# Summaries by TF.
rows = []

for tf, z in m.groupby("TF", sort=False):
    z = z.set_index("cell_type").reindex(celltypes)

    if z["delta_activity"].isna().any():
        raise ValueError(f"Missing cell-type activity result for {tf}")

    concord = z["direction_concordant"].astype(bool)
    fdr = z["FDR_support"].astype(bool)
    nom = z["nominal_support"].astype(bool)

    # Direction + significance support.
    concord_fdr = concord & fdr
    concord_nom = concord & nom

    cd14 = z.loc["CD14+ Monocyte"]
    cd16 = z.loc["CD16+ Monocyte"]

    monocyte_concordant_n = int(cd14["direction_concordant"]) + int(cd16["direction_concordant"])
    monocyte_fdr_concordant_n = int(
        bool(cd14["direction_concordant"]) and bool(cd14["FDR_support"])
    ) + int(
        bool(cd16["direction_concordant"]) and bool(cd16["FDR_support"])
    )

    lymphoid_types = [
        "Naive CD4+ T",
        "Naive CD8+ T",
        "TrB",
        "B",
        "NK",
        "Memory CD4+ T",
        "Cytotoxic T",
        "Treg",
        "Memory B",
    ]

    lymph = z.loc[lymphoid_types]
    lymph_concord_n = int(lymph["direction_concordant"].astype(bool).sum())
    lymph_fdr_concord_n = int(
        (
            lymph["direction_concordant"].astype(bool)
            & lymph["FDR_support"].astype(bool)
        ).sum()
    )

    rows.append({
        "TF": tf,
        "celltypes_evaluable": len(z),
        "celltypes_direction_concordant": int(concord.sum()),
        "celltypes_nominal_concordant": int(concord_nom.sum()),
        "celltypes_FDR_concordant": int(concord_fdr.sum()),
        "celltypes_FDR_any_direction": int(fdr.sum()),
        "CD14_delta": float(cd14["delta_activity"]),
        "CD14_FDR": float(cd14["FDR"]),
        "CD14_direction_concordant": bool(cd14["direction_concordant"]),
        "CD14_FDR_support": bool(cd14["FDR_support"]),
        "CD16_delta": float(cd16["delta_activity"]),
        "CD16_FDR": float(cd16["FDR"]),
        "CD16_direction_concordant": bool(cd16["direction_concordant"]),
        "CD16_FDR_support": bool(cd16["FDR_support"]),
        "monocyte_direction_concordant_n": monocyte_concordant_n,
        "monocyte_FDR_concordant_n": monocyte_fdr_concordant_n,
        "lymphoid_direction_concordant_n": lymph_concord_n,
        "lymphoid_FDR_concordant_n": lymph_fdr_concord_n,
    })

support = pd.DataFrame(rows)

# Classification rules are descriptive and fixed here.
def classify(row):
    broad = (
        row["celltypes_direction_concordant"] >= 8
        and row["celltypes_FDR_concordant"] >= 4
    )
    monocyte_strong = (
        row["monocyte_direction_concordant_n"] == 2
        and row["monocyte_FDR_concordant_n"] >= 1
    )
    monocyte_dominant = (
        monocyte_strong
        and row["lymphoid_direction_concordant_n"] <= 4
    )
    lymphoid_skewed = (
        row["lymphoid_direction_concordant_n"] >= 6
        and row["monocyte_direction_concordant_n"] <= 1
    )
    divergent = row["celltypes_direction_concordant"] <= 3

    if broad:
        return "broadly_shared"
    if monocyte_dominant:
        return "monocyte_dominant"
    if lymphoid_skewed:
        return "lymphoid_skewed"
    if monocyte_strong:
        return "monocyte_supported_mixed"
    if divergent:
        return "context_divergent"
    return "intermediate_mixed"

support["cellular_regulatory_class"] = support.apply(classify, axis=1)

# Merge frozen bulk evidence.
e = d.merge(support, on="TF", how="left", validate="one_to_one")

rna_keep = r[[
    "TF",
    "rnaseq_delta_activity",
    "rnaseq_P",
    "rnaseq_FDR",
    "direction_concordant"
]].copy()

rna_keep = rna_keep.rename(columns={
    "direction_concordant": "GSE155925_direction_concordant"
})

e = e.merge(rna_keep, on="TF", how="left", validate="one_to_one")

# Explicit support tiers for the final regulatory evidence table.
def cellular_tier(row):
    if (
        row["CD16_direction_concordant"]
        and row["CD16_FDR_support"]
        and row["CD14_direction_concordant"]
        and row["CD14_FDR_support"]
    ):
        return "TierA_both_monocyte_FDR_concordant"
    if (
        row["CD16_direction_concordant"]
        and row["CD16_FDR_support"]
    ):
        return "TierB_CD16_FDR_concordant"
    if (
        row["CD14_direction_concordant"]
        and row["CD14_FDR_support"]
    ):
        return "TierC_CD14_FDR_concordant"
    if row["monocyte_direction_concordant_n"] == 2:
        return "TierD_both_monocyte_directional"
    return "TierE_weak_or_divergent"

e["cellular_support_tier"] = e.apply(cellular_tier, axis=1)

# Sort by strongest cellular support, then frozen evidence.
tier_order = {
    "TierA_both_monocyte_FDR_concordant": 0,
    "TierB_CD16_FDR_concordant": 1,
    "TierC_CD14_FDR_concordant": 2,
    "TierD_both_monocyte_directional": 3,
    "TierE_weak_or_divergent": 4,
}
e["_tier_order"] = e["cellular_support_tier"].map(tier_order)

e = e.sort_values(
    [
        "_tier_order",
        "celltypes_FDR_concordant",
        "celltypes_direction_concordant",
        "regulatory_driver_score"
    ],
    ascending=[True, False, False, False]
).drop(columns="_tier_order")

master_file = OUT / "REGULATORY_DRIVER_MULTILAYER_CELLULAR_SUPPORT_MASTER_v1.0.tsv"
e.to_csv(master_file, sep="\t", index=False)

class_counts = (
    e.groupby("cellular_regulatory_class")
     .size()
     .reset_index(name="n_drivers")
     .sort_values("n_drivers", ascending=False)
)

class_counts.to_csv(
    OUT / "REGULATORY_DRIVER_CELLULAR_CLASS_COUNTS_v1.0.tsv",
    sep="\t", index=False
)

tier_counts = (
    e.groupby("cellular_support_tier")
     .size()
     .reset_index(name="n_drivers")
)

tier_counts.to_csv(
    OUT / "REGULATORY_DRIVER_CELLULAR_SUPPORT_TIER_COUNTS_v1.0.tsv",
    sep="\t", index=False
)

# Compact cell-type support matrix for plotting.
mat = m.pivot(index="TF", columns="cell_type", values="delta_activity")
mat = mat.reindex(index=d["TF"], columns=celltypes)
mat.to_csv(
    OUT / "REGULATORY_DRIVER_37x11_CELLTYPE_DELTA_ACTIVITY_MATRIX_v1.0.tsv",
    sep="\t"
)

sig = m.copy()
sig["signed_support"] = np.where(
    sig["FDR_support"].astype(bool),
    np.where(sig["direction_concordant"].astype(bool), 1, -1),
    0
)
sigmat = sig.pivot(index="TF", columns="cell_type", values="signed_support")
sigmat = sigmat.reindex(index=d["TF"], columns=celltypes)
sigmat.to_csv(
    OUT / "REGULATORY_DRIVER_37x11_CELLTYPE_SIGNED_SUPPORT_MATRIX_v1.0.tsv",
    sep="\t"
)

print("=== CELLULAR REGULATORY CLASS COUNTS ===")
print(class_counts.to_string(index=False))

print("\n=== CELLULAR SUPPORT TIERS ===")
print(tier_counts.to_string(index=False))

show_cols = [
    "TF",
    "regulatory_class",
    "RSV_replication_tier",
    "cellular_regulatory_class",
    "cellular_support_tier",
    "celltypes_direction_concordant",
    "celltypes_FDR_concordant",
    "CD14_delta",
    "CD14_FDR",
    "CD16_delta",
    "CD16_FDR",
    "GSE155925_direction_concordant",
    "frozen170_targets",
]

print("\n=== MULTILAYER DRIVER SUPPORT ===")
print(e[show_cols].to_string(index=False))

print("\nWritten:")
print(master_file)
print("\nMULTILAYER CELLULAR SUPPORT MATRIX COMPLETE.")
