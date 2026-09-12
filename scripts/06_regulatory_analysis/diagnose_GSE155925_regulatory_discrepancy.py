from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr

ROOT = Path.home() / "Influenza_RSV_Project"
BASE = ROOT / "results/regulatory_driver_analysis"

MICRO = BASE / "bulk/differential_activity/GPL6884_RSV_acute_vs_control_TF_activity_limma_v1.0.tsv"
RNA = BASE / "rnaseq_validation/GSE155925_all_TF_activity_limma_v1.0.tsv"
DRIVERS = BASE / "tables/REGULATORY_HIGH_CONFIDENCE_DRIVERS_v1.0.tsv"

OUT = BASE / "rnaseq_validation/diagnostics"
OUT.mkdir(parents=True, exist_ok=True)

m = pd.read_csv(MICRO, sep="\t")
r = pd.read_csv(RNA, sep="\t")
d = pd.read_csv(DRIVERS, sep="\t")

x = m.merge(r, on="TF", how="inner")

# Global evaluable comparison
pear = pearsonr(x["delta_activity"], x["rnaseq_delta_activity"])
spear = spearmanr(x["delta_activity"], x["rnaseq_delta_activity"])
same = np.sign(x["delta_activity"]) == np.sign(x["rnaseq_delta_activity"])

# >=10 target sensitivity
x10 = x[x["measurable_targets"] >= 10].copy()
pear10 = pearsonr(x10["delta_activity"], x10["rnaseq_delta_activity"])
spear10 = spearmanr(x10["delta_activity"], x10["rnaseq_delta_activity"])
same10 = np.sign(x10["delta_activity"]) == np.sign(x10["rnaseq_delta_activity"])

# Frozen 37
f37 = d[["TF","RSV_replication_tier","frozen170_targets","rsv_delta"]].merge(
    r, on="TF", how="left"
)
f37["direction_concordant"] = np.sign(f37["rsv_delta"]) == np.sign(f37["rnaseq_delta_activity"])

# Tier1 subset
tier1 = f37[f37["RSV_replication_tier"].astype(str).str.startswith("Tier1_")].copy()

summary = pd.DataFrame([
    {
        "analysis":"all_evaluable_TFs",
        "n":len(x),
        "pearson_r":pear.statistic,
        "pearson_p":pear.pvalue,
        "spearman_rho":spear.statistic,
        "spearman_p":spear.pvalue,
        "direction_concordant_n":int(same.sum()),
        "direction_concordant_pct":float(same.mean()*100),
    },
    {
        "analysis":"min10_TFs",
        "n":len(x10),
        "pearson_r":pear10.statistic,
        "pearson_p":pear10.pvalue,
        "spearman_rho":spear10.statistic,
        "spearman_p":spear10.pvalue,
        "direction_concordant_n":int(same10.sum()),
        "direction_concordant_pct":float(same10.mean()*100),
    },
    {
        "analysis":"frozen37",
        "n":len(f37),
        "pearson_r":pearsonr(f37["rsv_delta"],f37["rnaseq_delta_activity"]).statistic,
        "pearson_p":pearsonr(f37["rsv_delta"],f37["rnaseq_delta_activity"]).pvalue,
        "spearman_rho":spearmanr(f37["rsv_delta"],f37["rnaseq_delta_activity"]).statistic,
        "spearman_p":spearmanr(f37["rsv_delta"],f37["rnaseq_delta_activity"]).pvalue,
        "direction_concordant_n":int(f37["direction_concordant"].sum()),
        "direction_concordant_pct":float(f37["direction_concordant"].mean()*100),
    },
    {
        "analysis":"frozen37_Tier1",
        "n":len(tier1),
        "pearson_r":pearsonr(tier1["rsv_delta"],tier1["rnaseq_delta_activity"]).statistic,
        "pearson_p":pearsonr(tier1["rsv_delta"],tier1["rnaseq_delta_activity"]).pvalue,
        "spearman_rho":spearmanr(tier1["rsv_delta"],tier1["rnaseq_delta_activity"]).statistic,
        "spearman_p":spearmanr(tier1["rsv_delta"],tier1["rnaseq_delta_activity"]).pvalue,
        "direction_concordant_n":int(tier1["direction_concordant"].sum()),
        "direction_concordant_pct":float(tier1["direction_concordant"].mean()*100),
    },
])

summary.to_csv(
    OUT / "GSE155925_regulatory_discrepancy_correlation_summary_v1.0.tsv",
    sep="\t", index=False
)

x.to_csv(
    OUT / "GSE155925_vs_GPL6884_allTF_effect_comparison_v1.0.tsv",
    sep="\t", index=False
)

f37.to_csv(
    OUT / "GSE155925_frozen37_effect_comparison_v1.0.tsv",
    sep="\t", index=False
)

print("=== GSE155925 REGULATORY DISCREPANCY SUMMARY ===")
print(summary.to_string(index=False))
