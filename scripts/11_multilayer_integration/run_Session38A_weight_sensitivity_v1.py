from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

ROOT = Path.home() / "Influenza_RSV_Project"

INFILE = ROOT / (
    "results/session38_additional_analysis/tables/"
    "Session38A_EVIDENCE_WEIGHTED_GENE_PRIORITY_MASTER_v1.0.tsv"
)

OUT = ROOT / "results/session38_additional_analysis/tables"
df = pd.read_csv(INFILE, sep="\t")

assert len(df) == 170
assert df["gene_symbol"].nunique() == 170

# Normalize existing component scores to 0-1.
components = pd.DataFrame({
    "gene_symbol": df["gene_symbol"],
    "rna": df["score_rna"] / 4.0,
    "protein": df["score_proteomics"] / 3.0,
    "crispr": df["score_crispr"] / 3.0,
    "regulatory": df["score_regulatory"] / 2.0,
    "leading_edge": df["score_leading_edge"] / 2.0,
    "architecture": df["score_architecture"] / 1.0,
})

# Missing independent assays remain missing for the adjusted denominator.
models = {
    "primary": {
        "rna": 4, "protein": 3, "crispr": 3,
        "regulatory": 2, "leading_edge": 2, "architecture": 1
    },
    "independent_heavy": {
        "rna": 4, "protein": 4, "crispr": 4,
        "regulatory": 1, "leading_edge": 1, "architecture": 0.5
    },
    "functional_heavy": {
        "rna": 3, "protein": 2, "crispr": 5,
        "regulatory": 1, "leading_edge": 1, "architecture": 1
    },
    "response_heavy": {
        "rna": 5, "protein": 4, "crispr": 2,
        "regulatory": 2, "leading_edge": 2, "architecture": 1
    },
}

always_available = {"rna", "regulatory", "leading_edge", "architecture"}

result = df[[
    "gene_symbol",
    "independent_modalities_supported_n",
    "discovery_iav_logFC",
]].copy()

for model_name, weights in models.items():

    raw = np.zeros(len(df), dtype=float)
    available_max = np.zeros(len(df), dtype=float)

    for comp, weight in weights.items():
        vals = components[comp]

        if comp in always_available:
            raw += vals.fillna(0).to_numpy() * weight
            available_max += weight

        elif comp == "protein":
            evaluable = df["proteomics_evaluable"].astype(bool).to_numpy()
            raw += vals.fillna(0).to_numpy() * weight
            available_max += evaluable * weight

        elif comp == "crispr":
            evaluable = df["crispr_iav_evaluable"].astype(bool).to_numpy()
            raw += vals.fillna(0).to_numpy() * weight
            available_max += evaluable * weight

    adjusted = 100 * raw / available_max

    result[f"{model_name}_raw"] = raw
    result[f"{model_name}_adjusted"] = adjusted

    # Preserve preregistered rank philosophy:
    # modality breadth first, then weighted raw score.
    temp = result[[
        "gene_symbol",
        "independent_modalities_supported_n",
        "discovery_iav_logFC",
    ]].copy()

    temp["raw"] = raw
    temp["adjusted"] = adjusted
    temp["abs_effect"] = temp["discovery_iav_logFC"].abs()

    temp = temp.sort_values(
        [
            "independent_modalities_supported_n",
            "raw",
            "adjusted",
            "abs_effect",
            "gene_symbol",
        ],
        ascending=[False, False, False, False, True]
    ).reset_index(drop=True)

    rank_map = {
        gene: rank
        for rank, gene in enumerate(temp["gene_symbol"], start=1)
    }

    result[f"{model_name}_rank"] = (
        result["gene_symbol"].map(rank_map)
    )

# ------------------------------------------------------------
# Rank robustness
# ------------------------------------------------------------
rank_cols = [f"{m}_rank" for m in models]

result["rank_min"] = result[rank_cols].min(axis=1)
result["rank_max"] = result[rank_cols].max(axis=1)
result["rank_range"] = result["rank_max"] - result["rank_min"]
result["rank_median"] = result[rank_cols].median(axis=1)

result["top10_models_n"] = sum(
    result[c] <= 10 for c in rank_cols
)

result["top20_models_n"] = sum(
    result[c] <= 20 for c in rank_cols
)

result["top30_models_n"] = sum(
    result[c] <= 30 for c in rank_cols
)

# ------------------------------------------------------------
# Spearman correlations
# ------------------------------------------------------------
corr_rows = []

for m1 in models:
    for m2 in models:
        if m1 >= m2:
            continue

        rho, p = spearmanr(
            result[f"{m1}_rank"],
            result[f"{m2}_rank"]
        )

        corr_rows.append({
            "model_1": m1,
            "model_2": m2,
            "spearman_rho": rho,
            "p_value": p,
        })

corr = pd.DataFrame(corr_rows)

# ------------------------------------------------------------
# Consensus ranking
# ------------------------------------------------------------
result = result.sort_values(
    [
        "independent_modalities_supported_n",
        "top10_models_n",
        "top20_models_n",
        "rank_median",
        "rank_range",
        "gene_symbol",
    ],
    ascending=[False, False, False, True, True, True]
).reset_index(drop=True)

result.insert(
    0,
    "consensus_rank",
    np.arange(1, len(result) + 1)
)

result.to_csv(
    OUT / "Session38A_WEIGHT_SENSITIVITY_GENE_RANKS_v1.0.tsv",
    sep="\t",
    index=False
)

corr.to_csv(
    OUT / "Session38A_WEIGHT_SENSITIVITY_CORRELATIONS_v1.0.tsv",
    sep="\t",
    index=False
)

print("=== SESSION 38A WEIGHT SENSITIVITY ===")
print()

print("Spearman rank correlations:")
print(corr.to_string(index=False))

print()
print("=== CONSENSUS TOP 30 ===")

show = [
    "consensus_rank",
    "gene_symbol",
    "independent_modalities_supported_n",
    "primary_rank",
    "independent_heavy_rank",
    "functional_heavy_rank",
    "response_heavy_rank",
    "rank_min",
    "rank_max",
    "rank_range",
    "top10_models_n",
    "top20_models_n",
]

print(result[show].head(30).to_string(index=False))

print()
print("=== GENES TOP-10 IN ALL FOUR MODELS ===")

stable10 = result[result["top10_models_n"] == 4]

if len(stable10):
    print(stable10[show].to_string(index=False))
else:
    print("NONE")

print()
print("=== GENES TOP-20 IN ALL FOUR MODELS ===")

stable20 = result[result["top20_models_n"] == 4]

if len(stable20):
    print(stable20[show].to_string(index=False))
else:
    print("NONE")

print()
print("=== BRIDGE GENE ROBUSTNESS ===")

bridges = {
    "KPNB1", "HERC5", "OTOF",
    "CHMP5", "TOP2A", "FCGR1B", "HIST2H2AC"
}

x = result[result["gene_symbol"].isin(bridges)]

print(
    x[show]
    .sort_values("consensus_rank")
    .to_string(index=False)
)

print()
print("Written:")
print(
    OUT /
    "Session38A_WEIGHT_SENSITIVITY_GENE_RANKS_v1.0.tsv"
)
print(
    OUT /
    "Session38A_WEIGHT_SENSITIVITY_CORRELATIONS_v1.0.tsv"
)
