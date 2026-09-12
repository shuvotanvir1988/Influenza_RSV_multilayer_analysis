from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path.home() / "Influenza_RSV_Project"
BASE = ROOT / "results/regulatory_driver_analysis"

NET = BASE / "design/COLLECTRI_HUMAN_FROZEN_v1.0.tsv.gz"
REG = BASE / "tables/REGULATORY_ARCHITECTURE_CLASSIFICATION_FROZEN_v1.0.tsv"
GENES = ROOT / "results/DS001_GSE38900/signature_analysis/tables/DS001_GSE38900_gene_evidence_master_final_classification.tsv"

OUT = BASE / "tables"
OUT.mkdir(parents=True, exist_ok=True)

net = pd.read_csv(NET, sep="\t")
reg = pd.read_csv(REG, sep="\t")
genes = pd.read_csv(GENES, sep="\t")

gene_classes = {
    "shared_core": "shared_replicated_core",
    "shared_influenza_amplified": "shared_replicated_influenza_amplified",
    "shared_RSV_amplified": "shared_replicated_RSV_amplified",
    "influenza_selective": "influenza_selective_primary",
    "RSV_selective": "RSV_selective_primary",
    "opposite_direction": "opposite_direction_primary",
}

members = {}
for short, col in gene_classes.items():
    if col not in genes.columns:
        raise ValueError(f"Missing frozen gene-class column: {col}")
    vals = genes.loc[genes[col] == True, "gene_symbol"].dropna().astype(str)
    members[short] = set(vals)

# Long-form TF-target-class map
rows = []
for short, gset in members.items():
    z = net[net["target"].astype(str).isin(gset)].copy()
    z["frozen_gene_class"] = short
    rows.append(z)

long_map = pd.concat(rows, ignore_index=True)
long_map.to_csv(
    OUT / "REGULATORY_TF_TARGET_FROZEN_GENE_CLASS_MAP_v1.0.tsv.gz",
    sep="\t", index=False, compression="gzip"
)

# Summarize each TF against frozen classes
summary_rows = []

for tf, g in net.groupby("source", sort=True):
    row = {"TF": tf}

    for short, gset in members.items():
        x = g[g["target"].astype(str).isin(gset)].copy()

        row[f"{short}_targets"] = int(x["target"].nunique())
        row[f"{short}_positive_edges"] = int((x["weight"] > 0).sum())
        row[f"{short}_negative_edges"] = int((x["weight"] < 0).sum())

        # Signed-balance score: positive minus negative edge count
        row[f"{short}_signed_balance"] = (
            row[f"{short}_positive_edges"] -
            row[f"{short}_negative_edges"]
        )

        denom = len(gset)
        row[f"{short}_coverage_fraction"] = (
            row[f"{short}_targets"] / denom if denom else np.nan
        )

    row["frozen170_targets"] = (
        row["shared_core_targets"]
        + row["shared_influenza_amplified_targets"]
        + row["shared_RSV_amplified_targets"]
    )
    row["frozen170_coverage_fraction"] = row["frozen170_targets"] / 170.0

    summary_rows.append(row)

tfmap = pd.DataFrame(summary_rows)

# Merge with frozen regulatory evidence
evidence = reg.merge(tfmap, on="TF", how="left")

# Primary evidence score is descriptive, not inferential.
# It rewards replicated regulatory activity and direct connectivity to the frozen 170-gene program.
tier_points = {
    "Tier1_significant_both_same_direction": 3,
    "Tier2_GPL6884_significant_same_direction": 2,
    "Tier3_GPL6884_significant_discordant": 0,
    "Tier4_replication_only_same_direction": 1,
    "Tier5_no_strong_replication": 0,
    "not_evaluable": 0,
}
evidence["replication_points"] = evidence["RSV_replication_tier"].map(tier_points).fillna(0)

class_points = {
    "shared_concordant": 3,
    "shared_influenza_amplified": 3,
    "shared_RSV_amplified": 3,
    "influenza_selective_supported": 2,
    "RSV_selective_supported": 2,
    "discordant_opposite": 2,
    "influenza_selective_unresolved": 1,
    "RSV_selective_unresolved": 1,
    "pathogen_differential_without_control_significance": 1,
    "no_primary_regulatory_signal": 0,
}
evidence["class_points"] = evidence["regulatory_class"].map(class_points).fillna(0)

# Connectivity tiers to the frozen 170-gene architecture
def connectivity_tier(n):
    if n >= 20:
        return "very_high"
    if n >= 10:
        return "high"
    if n >= 5:
        return "moderate"
    if n >= 1:
        return "low"
    return "none"

evidence["frozen170_connectivity_tier"] = evidence["frozen170_targets"].apply(connectivity_tier)

# Composite descriptive score; not used for statistical significance.
evidence["regulatory_driver_score"] = (
    evidence["replication_points"]
    + evidence["class_points"]
    + np.minimum(evidence["frozen170_targets"], 20) / 5.0
    + evidence["min10_all_primary"].astype(int)
)

# High-confidence candidate definition
evidence["high_confidence_regulatory_driver"] = (
    evidence["high_conf_shared"]
    & (evidence["frozen170_targets"] >= 5)
    & evidence["min10_all_primary"]
)

# Pathogen-biased candidates with direct connectivity
evidence["high_confidence_pathogen_biased_driver"] = (
    evidence["high_conf_pathogen_bias"]
    & (evidence["frozen170_targets"] >= 3)
    & evidence["min10_all_primary"]
)

evidence = evidence.sort_values(
    ["high_confidence_regulatory_driver",
     "regulatory_driver_score",
     "frozen170_targets",
     "rsv_FDR",
     "flu_FDR"],
    ascending=[False, False, False, True, True]
)

master = OUT / "REGULATORY_DRIVER_EVIDENCE_MASTER_v1.0.tsv"
evidence.to_csv(master, sep="\t", index=False)

hc = evidence[evidence["high_confidence_regulatory_driver"]].copy()
hc.to_csv(
    OUT / "REGULATORY_HIGH_CONFIDENCE_DRIVERS_v1.0.tsv",
    sep="\t", index=False
)

biased = evidence[evidence["high_confidence_pathogen_biased_driver"]].copy()
biased.to_csv(
    OUT / "REGULATORY_HIGH_CONFIDENCE_PATHOGEN_BIASED_DRIVERS_v1.0.tsv",
    sep="\t", index=False
)

# Summaries
summary = pd.DataFrame([{
    "regulators_total": len(evidence),
    "high_confidence_shared_drivers": int(evidence["high_confidence_regulatory_driver"].sum()),
    "high_confidence_pathogen_biased_drivers": int(evidence["high_confidence_pathogen_biased_driver"].sum()),
    "TFs_targeting_at_least_1_frozen170_gene": int((evidence["frozen170_targets"] >= 1).sum()),
    "TFs_targeting_at_least_5_frozen170_genes": int((evidence["frozen170_targets"] >= 5).sum()),
    "TFs_targeting_at_least_10_frozen170_genes": int((evidence["frozen170_targets"] >= 10).sum()),
    "TFs_targeting_at_least_20_frozen170_genes": int((evidence["frozen170_targets"] >= 20).sum()),
}])

summary.to_csv(
    OUT / "REGULATORY_DRIVER_EVIDENCE_SUMMARY_v1.0.tsv",
    sep="\t", index=False
)

print("=== FROZEN GENE CLASS SIZES ===")
for k, v in members.items():
    print(k, len(v))

print("\n=== DRIVER SUMMARY ===")
print(summary.to_string(index=False))

print("\n=== TOP HIGH-CONFIDENCE REGULATORY DRIVERS ===")
cols = [
    "TF",
    "regulatory_class",
    "RSV_replication_tier",
    "flu_delta",
    "flu_FDR",
    "rsv_delta",
    "rsv_FDR",
    "flu_vs_rsv_delta",
    "flu_vs_rsv_FDR",
    "frozen170_targets",
    "shared_core_targets",
    "shared_influenza_amplified_targets",
    "shared_RSV_amplified_targets",
    "min10_all_primary",
    "regulatory_driver_score"
]
print(hc[cols].head(50).to_string(index=False))

print("\n=== TOP PATHOGEN-BIASED DRIVERS ===")
print(biased[cols].head(40).to_string(index=False))

print("\nWritten:")
print(master)
print("\nREGULATORY DRIVER EVIDENCE MAPPING COMPLETE.")
