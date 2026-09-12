from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path.home() / "Influenza_RSV_Project"

MASTER = ROOT / (
    "results/influenza_multilayer_integration/gene_level/"
    "SESSION27_INFLUENZA_OBJECTIVE_CANDIDATE_CLASSIFICATION_v1.0.tsv"
)

LE = ROOT / (
    "results/figure3_influenza_centered/frozen_submission_v4/tables/"
    "Figure3_SUBMISSION_v4_leading_edge_recurrence.tsv"
)

OUT = ROOT / "results/session38_additional_analysis/tables"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(MASTER, sep="\t")
le = pd.read_csv(LE, sep="\t")

assert len(df) == 170
assert df["gene_symbol"].nunique() == 170

# ------------------------------------------------------------
# RNA SCORE
# ------------------------------------------------------------
def score_rna(row):
    fdr_n = int(row["rnaseq_n_seasons_FDR05_direction_match"])
    dir_n = int(row["n_seasons_direction_match"])

    if fdr_n >= 4:
        return 4.0
    elif fdr_n >= 3:
        return 3.0
    elif dir_n >= 3:
        return 2.0
    elif dir_n >= 1:
        return 1.0
    return 0.0

df["score_rna"] = df.apply(score_rna, axis=1)

# ------------------------------------------------------------
# PROTEOMICS SCORE
# ------------------------------------------------------------
def score_protein(row):
    if not bool(row["proteomics_evaluable"]):
        return np.nan

    if bool(row["proteomics_strong_support"]):
        return 3.0

    if bool(row["proteomics_nominal_support"]):
        return 2.0

    if bool(row["proteomics_direction_only"]):
        return 1.0

    return 0.0

df["score_proteomics"] = df.apply(score_protein, axis=1)

# ------------------------------------------------------------
# CRISPR SCORE
# ------------------------------------------------------------
def score_crispr(row):
    if not bool(row["crispr_iav_evaluable"]):
        return np.nan

    if bool(row["crispr_high_confidence"]):
        return 3.0

    if bool(row["crispr_support"]):
        return 2.0

    return 0.0

df["score_crispr"] = df.apply(score_crispr, axis=1)

# ------------------------------------------------------------
# REGULATORY SCORE
# ------------------------------------------------------------
def score_regulatory(n):
    n = int(n)

    if n >= 6:
        return 2.0
    elif n >= 3:
        return 1.0
    elif n >= 1:
        return 0.5
    return 0.0

df["score_regulatory"] = (
    df["influenza_sig_driver_count"]
    .fillna(0)
    .apply(score_regulatory)
)

# ------------------------------------------------------------
# LEADING-EDGE RECURRENCE
# ------------------------------------------------------------
le_genes = set(le["leading_edge_gene"].astype(str))

df["leading_edge_recurrent"] = (
    df["gene_symbol"].astype(str).isin(le_genes)
)

df["score_leading_edge"] = np.where(
    df["leading_edge_recurrent"], 2.0, 0.0
)

# ------------------------------------------------------------
# INFLUENZA-AMPLIFIED ARCHITECTURE
# ------------------------------------------------------------
df["influenza_amplified"] = (
    df["proteomics_architecture_class"]
    .astype(str)
    .str.lower()
    .eq("influenza_amplified")
)

df["score_architecture"] = np.where(
    df["influenza_amplified"], 1.0, 0.0
)

# ------------------------------------------------------------
# RAW SCORE
# Missing proteomics/CRISPR contribute no observed points.
# Missingness is handled explicitly in denominator below.
# ------------------------------------------------------------
score_cols = [
    "score_rna",
    "score_proteomics",
    "score_crispr",
    "score_regulatory",
    "score_leading_edge",
    "score_architecture",
]

df["evidence_score_raw"] = (
    df[score_cols]
    .fillna(0)
    .sum(axis=1)
)

# ------------------------------------------------------------
# AVAILABLE MAXIMUM
# RNA 4 + regulatory 2 + LE 2 + architecture 1 = 9 always.
# Protein adds 3 if evaluable.
# CRISPR adds 3 if evaluable.
# ------------------------------------------------------------
df["evidence_score_available_max"] = (
    9.0
    + np.where(df["proteomics_evaluable"], 3.0, 0.0)
    + np.where(df["crispr_iav_evaluable"], 3.0, 0.0)
)

df["evidence_score_coverage_adjusted"] = (
    100.0
    * df["evidence_score_raw"]
    / df["evidence_score_available_max"]
)

# ------------------------------------------------------------
# BIOLOGICAL AXES
# ------------------------------------------------------------
response_support = (
    df["independent_rnaseq_supported"].astype(bool)
    & (
        df["independent_proteomics_supported"].astype(bool)
        | (df["score_regulatory"] > 0)
        | df["leading_edge_recurrent"]
    )
)

dependency_support = df["crispr_support"].astype(bool)

df["response_priority"] = response_support
df["dependency_priority"] = dependency_support

df["response_dependency_bridge"] = (
    df["independent_rnaseq_supported"].astype(bool)
    & dependency_support
)

# ------------------------------------------------------------
# RANK
# ------------------------------------------------------------
df["abs_discovery_iav_logFC"] = df["discovery_iav_logFC"].abs()

df = df.sort_values(
    by=[
        "independent_modalities_supported_n",
        "evidence_score_raw",
        "evidence_score_coverage_adjusted",
        "score_rna",
        "abs_discovery_iav_logFC",
        "gene_symbol",
    ],
    ascending=[False, False, False, False, False, True],
).reset_index(drop=True)

df.insert(0, "session38_rank", np.arange(1, len(df) + 1))

# ------------------------------------------------------------
# OUTPUTS
# ------------------------------------------------------------
master_out = OUT / "Session38A_EVIDENCE_WEIGHTED_GENE_PRIORITY_MASTER_v1.0.tsv"
df.to_csv(master_out, sep="\t", index=False)

top_out = OUT / "Session38A_TOP30_EVIDENCE_WEIGHTED_GENES_v1.0.tsv"
df.head(30).to_csv(top_out, sep="\t", index=False)

bridge_out = OUT / "Session38A_RESPONSE_DEPENDENCY_BRIDGE_GENES_v1.0.tsv"
df[df["response_dependency_bridge"]].to_csv(
    bridge_out, sep="\t", index=False
)

summary = pd.DataFrame({
    "metric": [
        "genes_total",
        "proteomics_evaluable",
        "crispr_evaluable",
        "leading_edge_recurrent",
        "influenza_amplified",
        "response_priority",
        "dependency_priority",
        "response_dependency_bridge",
        "raw_score_ge_8",
        "raw_score_ge_10",
    ],
    "value": [
        len(df),
        int(df["proteomics_evaluable"].sum()),
        int(df["crispr_iav_evaluable"].sum()),
        int(df["leading_edge_recurrent"].sum()),
        int(df["influenza_amplified"].sum()),
        int(df["response_priority"].sum()),
        int(df["dependency_priority"].sum()),
        int(df["response_dependency_bridge"].sum()),
        int((df["evidence_score_raw"] >= 8).sum()),
        int((df["evidence_score_raw"] >= 10).sum()),
    ],
})

summary_out = OUT / "Session38A_EVIDENCE_WEIGHTED_SUMMARY_v1.0.tsv"
summary.to_csv(summary_out, sep="\t", index=False)

print("=== SESSION 38A COMPLETE ===")
print(f"Genes: {len(df)}")
print(f"Proteomics evaluable: {int(df['proteomics_evaluable'].sum())}")
print(f"CRISPR evaluable: {int(df['crispr_iav_evaluable'].sum())}")
print(f"Leading-edge recurrent: {int(df['leading_edge_recurrent'].sum())}")
print(f"Influenza-amplified: {int(df['influenza_amplified'].sum())}")
print(f"Response-priority: {int(df['response_priority'].sum())}")
print(f"Dependency-priority: {int(df['dependency_priority'].sum())}")
print(f"Response-dependency bridges: {int(df['response_dependency_bridge'].sum())}")

print()
print("=== TOP 20 ===")

cols = [
    "session38_rank",
    "gene_symbol",
    "proteomics_architecture_class",
    "independent_modalities_supported_n",
    "score_rna",
    "score_proteomics",
    "score_crispr",
    "score_regulatory",
    "score_leading_edge",
    "score_architecture",
    "evidence_score_raw",
    "evidence_score_available_max",
    "evidence_score_coverage_adjusted",
    "candidate_tier",
    "response_dependency_bridge",
]

print(df[cols].head(20).to_string(index=False))

print()
print("=== RESPONSE-DEPENDENCY BRIDGES ===")

bridges = df[df["response_dependency_bridge"]]

if len(bridges):
    print(bridges[cols].to_string(index=False))
else:
    print("NONE")

print()
print("Written:")
print(master_out)
print(top_out)
print(bridge_out)
print(summary_out)
