from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path.home() / "Influenza_RSV_Project"
BASE = ROOT / "results/regulatory_driver_analysis"

FLU = BASE / "bulk/differential_activity/GPL6884_InfluenzaA_vs_control_TF_activity_limma_v1.0.tsv"
RSV = BASE / "bulk/differential_activity/GPL6884_RSV_acute_vs_control_TF_activity_limma_v1.0.tsv"
DIR = BASE / "bulk/differential_activity/GPL6884_InfluenzaA_vs_RSVacute_TF_activity_limma_v1.0.tsv"
REP = BASE / "replication/RSV_TF_crossplatform_replication_master_v1.0.tsv"

OUT = BASE / "tables"
OUT.mkdir(parents=True, exist_ok=True)

flu = pd.read_csv(FLU, sep="\t")
rsv = pd.read_csv(RSV, sep="\t")
direct = pd.read_csv(DIR, sep="\t")
rep = pd.read_csv(REP, sep="\t")

keep_flu = flu[["TF","delta_activity","P_value","FDR","measurable_targets","eligible_min10"]].copy()
keep_flu = keep_flu.rename(columns={
    "delta_activity":"flu_delta",
    "P_value":"flu_P",
    "FDR":"flu_FDR",
    "measurable_targets":"flu_targets",
    "eligible_min10":"flu_min10"
})

keep_rsv = rsv[["TF","delta_activity","P_value","FDR","measurable_targets","eligible_min10"]].copy()
keep_rsv = keep_rsv.rename(columns={
    "delta_activity":"rsv_delta",
    "P_value":"rsv_P",
    "FDR":"rsv_FDR",
    "measurable_targets":"rsv_targets",
    "eligible_min10":"rsv_min10"
})

keep_dir = direct[["TF","delta_activity","P_value","FDR","measurable_targets","eligible_min10"]].copy()
keep_dir = keep_dir.rename(columns={
    "delta_activity":"flu_vs_rsv_delta",
    "P_value":"flu_vs_rsv_P",
    "FDR":"flu_vs_rsv_FDR",
    "measurable_targets":"direct_targets",
    "eligible_min10":"direct_min10"
})

df = keep_flu.merge(keep_rsv,on="TF",how="inner").merge(keep_dir,on="TF",how="inner")

rep_cols = [
    "TF",
    "delta_activity_10558",
    "FDR_10558",
    "same_direction",
    "significant_both",
    "replicated_significant_same_direction",
    "measurable_targets_10558"
]
rep2 = rep[rep_cols].copy()
df = df.merge(rep2,on="TF",how="left")

# Primary significance flags
df["flu_sig"] = df["flu_FDR"] < 0.05
df["rsv_sig"] = df["rsv_FDR"] < 0.05
df["direct_sig"] = df["flu_vs_rsv_FDR"] < 0.05

df["same_infection_direction"] = np.sign(df["flu_delta"]) == np.sign(df["rsv_delta"])
df["opposite_infection_direction"] = np.sign(df["flu_delta"]) == -np.sign(df["rsv_delta"])

# Deterministic primary regulatory class
def classify(row):
    flu_sig = row["flu_sig"]
    rsv_sig = row["rsv_sig"]
    direct_sig = row["direct_sig"]
    same = row["same_infection_direction"]
    opp = row["opposite_infection_direction"]
    ddir = row["flu_vs_rsv_delta"]

    # Both infections significant
    if flu_sig and rsv_sig:
        if opp:
            return "discordant_opposite"
        if same:
            if direct_sig:
                if ddir > 0:
                    return "shared_influenza_amplified"
                elif ddir < 0:
                    return "shared_RSV_amplified"
            return "shared_concordant"
        return "shared_unresolved"

    # Influenza only
    if flu_sig and not rsv_sig:
        if direct_sig and ddir > 0:
            return "influenza_selective_supported"
        return "influenza_selective_unresolved"

    # RSV only
    if rsv_sig and not flu_sig:
        if direct_sig and ddir < 0:
            return "RSV_selective_supported"
        return "RSV_selective_unresolved"

    # Neither infection individually significant
    if not flu_sig and not rsv_sig:
        if direct_sig:
            return "pathogen_differential_without_control_significance"
        return "no_primary_regulatory_signal"

    return "unclassified"

df["regulatory_class"] = df.apply(classify, axis=1)

# RSV replication tier
def rep_tier(row):
    if pd.isna(row["delta_activity_10558"]):
        return "not_evaluable"
    same = bool(row["same_direction"])
    sig6884 = row["rsv_FDR"] < 0.05
    sig10558 = row["FDR_10558"] < 0.05
    if sig6884 and sig10558 and same:
        return "Tier1_significant_both_same_direction"
    if sig6884 and same:
        return "Tier2_GPL6884_significant_same_direction"
    if sig6884 and not same:
        return "Tier3_GPL6884_significant_discordant"
    if (not sig6884) and sig10558 and same:
        return "Tier4_replication_only_same_direction"
    return "Tier5_no_strong_replication"

df["RSV_replication_tier"] = df.apply(rep_tier, axis=1)

# Sensitivity flag
df["min10_all_primary"] = (
    (df["flu_targets"] >= 10) &
    (df["rsv_targets"] >= 10) &
    (df["direct_targets"] >= 10)
)

# High-confidence flags
df["high_conf_shared"] = (
    df["regulatory_class"].isin([
        "shared_concordant",
        "shared_influenza_amplified",
        "shared_RSV_amplified"
    ])
    & df["RSV_replication_tier"].isin([
        "Tier1_significant_both_same_direction",
        "Tier2_GPL6884_significant_same_direction"
    ])
)

df["high_conf_pathogen_bias"] = (
    df["regulatory_class"].isin([
        "shared_influenza_amplified",
        "shared_RSV_amplified",
        "influenza_selective_supported",
        "RSV_selective_supported",
        "discordant_opposite"
    ])
)

# Save master table
master_file = OUT / "REGULATORY_ARCHITECTURE_CLASSIFICATION_FROZEN_v1.0.tsv"
df.to_csv(master_file, sep="\t", index=False)

# Summaries
class_summary = (
    df.groupby("regulatory_class", dropna=False)
      .agg(
          n_TFs=("TF","size"),
          min10_all_primary=("min10_all_primary","sum"),
          Tier1_RSV=("RSV_replication_tier",lambda x:(x=="Tier1_significant_both_same_direction").sum()),
          Tier2_RSV=("RSV_replication_tier",lambda x:(x=="Tier2_GPL6884_significant_same_direction").sum()),
          high_conf_shared=("high_conf_shared","sum")
      )
      .reset_index()
      .sort_values("n_TFs", ascending=False)
)

class_summary_file = OUT / "REGULATORY_ARCHITECTURE_CLASS_COUNTS_v1.0.tsv"
class_summary.to_csv(class_summary_file, sep="\t", index=False)

tier_summary = (
    df.groupby("RSV_replication_tier")
      .size()
      .reset_index(name="n_TFs")
      .sort_values("n_TFs", ascending=False)
)
tier_summary_file = OUT / "RSV_REGULATORY_REPLICATION_TIER_COUNTS_v1.0.tsv"
tier_summary.to_csv(tier_summary_file, sep="\t", index=False)

# Leading high-confidence tables
shared_hc = df[df["high_conf_shared"]].copy()
shared_hc = shared_hc.sort_values(["rsv_FDR","flu_FDR","FDR_10558"])
shared_hc.to_csv(
    OUT / "REGULATORY_HIGH_CONFIDENCE_SHARED_v1.0.tsv",
    sep="\t", index=False
)

bias_hc = df[df["high_conf_pathogen_bias"]].copy()
bias_hc = bias_hc.sort_values(["flu_vs_rsv_FDR","flu_FDR","rsv_FDR"])
bias_hc.to_csv(
    OUT / "REGULATORY_PATHOGEN_BIASED_CANDIDATES_v1.0.tsv",
    sep="\t", index=False
)

print("=== REGULATORY CLASS COUNTS ===")
print(class_summary.to_string(index=False))

print("\n=== RSV REPLICATION TIERS ===")
print(tier_summary.to_string(index=False))

print("\n=== HIGH-CONFIDENCE SHARED TFs ===")
show = shared_hc[[
    "TF","regulatory_class",
    "flu_delta","flu_FDR",
    "rsv_delta","rsv_FDR",
    "flu_vs_rsv_delta","flu_vs_rsv_FDR",
    "delta_activity_10558","FDR_10558",
    "RSV_replication_tier",
    "min10_all_primary"
]]
print(show.head(50).to_string(index=False))

print("\nWritten:")
print(master_file)
print(class_summary_file)
print(tier_summary_file)
print("\nREGULATORY ARCHITECTURE CLASSIFICATION COMPLETE.")
