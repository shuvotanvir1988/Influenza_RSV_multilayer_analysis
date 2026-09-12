from pathlib import Path
import csv
import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

SEED = 20260812
N_PERM = 50000

ROOT = Path(
    "results/proteomics_validation/Influenza-P01"
)

DRIVERS = Path(
    "results/regulatory_driver_analysis/rnaseq_validation/"
    "GSE155925_frozen37_driver_validation_v1.0.tsv"
)

NETWORK = Path(
    "results/regulatory_driver_analysis/design/"
    "COLLECTRI_HUMAN_FROZEN_v1.0.tsv.gz"
)

PRIMARY_RANK = ROOT / (
    "pathway_validation/"
    "Influenza-P01_proteome_gene_level_median_t_rank_v1.0.tsv"
)

NONICU = ROOT / (
    "differential_proteomics/"
    "Influenza-P01_nonICU_infection_limma_v1.0.tsv"
)

ANNOT = ROOT / (
    "analysis_inputs/"
    "Influenza-P01_human_assay_annotation_v1.0.tsv"
)

OUT = ROOT / "regulatory_validation"
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. READ INPUTS
# ============================================================

drivers = pd.read_csv(
    DRIVERS,
    sep="\t"
)

net = pd.read_csv(
    NETWORK,
    sep="\t",
    compression="gzip"
)

rank = pd.read_csv(
    PRIMARY_RANK,
    sep="\t"
)

nonicu = pd.read_csv(
    NONICU,
    sep="\t"
)

ann = pd.read_csv(
    ANNOT,
    sep="\t"
)

assert len(drivers) == 37
assert drivers["TF"].nunique() == 37
assert not net.duplicated(
    subset=["source", "target"]
).any()

print("=== INPUT QC ===")
print("Frozen drivers:", len(drivers))
print("CollecTRI edges:", len(net))
print(
    "Primary proteomic ranked genes:",
    rank["gene_symbol"].nunique()
)

# ============================================================
# 2. PRIMARY GENE-LEVEL T STATISTICS
# ============================================================

primary_gene = (
    rank[
        ["gene_symbol", "protein_median_t"]
    ]
    .drop_duplicates("gene_symbol")
    .copy()
)

primary_gene["gene_symbol"] = (
    primary_gene["gene_symbol"]
    .astype(str)
    .str.strip()
)

primary_t = dict(
    zip(
        primary_gene["gene_symbol"],
        primary_gene["protein_median_t"]
    )
)

# ============================================================
# 3. BUILD NON-ICU GENE-LEVEL RANK USING SAME FROZEN RULE
# ============================================================

# The non-ICU limma result already contains the full assay annotation,
# including EntrezGeneSymbol. Do not merge annotation a second time.
if "EntrezGeneSymbol" not in nonicu.columns:
    raise SystemExit(
        "STOP: EntrezGeneSymbol missing from non-ICU limma result."
    )

x = nonicu.copy()

expanded = []

for _, r in x.iterrows():

    raw = str(r["EntrezGeneSymbol"]).strip()

    if not raw or raw.lower() == "nan":
        continue

    for symbol in raw.split("|"):

        symbol = symbol.strip()

        if symbol:

            expanded.append(
                {
                    "gene_symbol": symbol,
                    "t": r["t"],
                    "SeqId": r["SeqId"],
                }
            )

expanded = pd.DataFrame(expanded)

nonicu_gene = (
    expanded
    .groupby("gene_symbol", as_index=False)
    .agg(
        protein_median_t=("t", "median"),
        n_assays=("SeqId", "nunique"),
    )
)

nonicu_t = dict(
    zip(
        nonicu_gene["gene_symbol"],
        nonicu_gene["protein_median_t"]
    )
)

# ============================================================
# 4. SIGNED TARGET-PROGRAM TEST
# ============================================================

rng = np.random.default_rng(SEED)

primary_universe = np.asarray(
    primary_gene["protein_median_t"],
    dtype=float
)

def permutation_test(
    target_t,
    weights,
    universe,
    n_perm=N_PERM,
):

    target_t = np.asarray(
        target_t,
        dtype=float
    )

    weights = np.asarray(
        weights,
        dtype=float
    )

    observed = np.mean(
        target_t * weights
    )

    n = len(weights)

    extreme = 0

    # Batched to avoid a large memory footprint.
    batch_size = 1000

    completed = 0

    while completed < n_perm:

        b = min(
            batch_size,
            n_perm - completed
        )

        null_scores = np.empty(
            b,
            dtype=float
        )

        for i in range(b):

            sampled = rng.choice(
                universe,
                size=n,
                replace=False
            )

            null_scores[i] = np.mean(
                sampled * weights
            )

        extreme += int(
            np.sum(
                np.abs(null_scores) >=
                abs(observed)
            )
        )

        completed += b

    p = (
        extreme + 1
    ) / (
        n_perm + 1
    )

    return observed, p

rows = []

for _, d in drivers.iterrows():

    tf = d["TF"]

    z = net.loc[
        net["source"] == tf,
        ["source", "target", "weight"]
    ].copy()

    z = z.loc[
        z["target"].isin(primary_t)
    ].copy()

    n_targets = z["target"].nunique()

    if n_targets < 10:

        rows.append(
            {
                "TF": tf,
                "n_measurable_targets": n_targets,
                "proteomic_program_score": np.nan,
                "permutation_P": np.nan,
            }
        )

        continue

    target_stats = np.asarray(
        [
            primary_t[g]
            for g in z["target"]
        ],
        dtype=float
    )

    weights = np.asarray(
        z["weight"],
        dtype=float
    )

    score, p = permutation_test(
        target_stats,
        weights,
        primary_universe
    )

    rows.append(
        {
            "TF": tf,
            "n_measurable_targets": n_targets,
            "n_positive_edges":
                int((weights > 0).sum()),
            "n_negative_edges":
                int((weights < 0).sum()),
            "proteomic_program_score": score,
            "permutation_P": p,
        }
    )

result = pd.DataFrame(rows)

# ============================================================
# 5. MERGE FROZEN TRANSCRIPTOMIC EXPECTATION
# ============================================================

keep_cols = [
    "TF",
    "flu_delta",
    "flu_FDR",
    "regulatory_class",
    "RSV_replication_tier",
    "regulatory_driver_score",
    "high_confidence_regulatory_driver",
    "high_confidence_pathogen_biased_driver",
]

result = result.merge(
    drivers[keep_cols],
    on="TF",
    how="left",
    validate="one_to_one"
)

# BH correction
pvals = result["permutation_P"].to_numpy()

order = np.argsort(pvals)
rank_order = np.empty_like(order)
rank_order[order] = np.arange(1, len(pvals) + 1)

q = pvals * len(pvals) / rank_order

# enforce monotonic BH in sorted order
sorted_q = q[order]
sorted_q = np.minimum.accumulate(
    sorted_q[::-1]
)[::-1]

q_final = np.empty_like(sorted_q)
q_final[order] = np.minimum(
    sorted_q,
    1.0
)

result["proteomic_FDR"] = q_final

result["direction_concordant"] = (
    np.sign(
        result["proteomic_program_score"]
    )
    ==
    np.sign(
        result["flu_delta"]
    )
)

def classify(r):

    concordant = bool(
        r["direction_concordant"]
    )

    fdr = r["proteomic_FDR"]
    p = r["permutation_P"]

    if concordant and fdr < 0.05:
        return "PROTEOMIC_REGULATORY_CONCORDANT_FDR"

    if concordant and p < 0.05:
        return "PROTEOMIC_REGULATORY_CONCORDANT_NOMINAL"

    if concordant:
        return "PROTEOMIC_REGULATORY_CONCORDANT_UNSUPPORTED"

    if fdr < 0.05:
        return "PROTEOMIC_REGULATORY_DISCORDANT_FDR"

    if p < 0.05:
        return "PROTEOMIC_REGULATORY_DISCORDANT_NOMINAL"

    return "PROTEOMIC_REGULATORY_DISCORDANT_UNSUPPORTED"

result["validation_category"] = (
    result.apply(
        classify,
        axis=1
    )
)

# ============================================================
# 6. FAMILY-LEVEL VALIDATION
# ============================================================

k = int(
    result["direction_concordant"].sum()
)

n = len(result)

binom = binomtest(
    k,
    n,
    p=0.5,
    alternative="greater"
)

rho = spearmanr(
    result["flu_delta"],
    result["proteomic_program_score"]
)

vc = result[
    "validation_category"
].value_counts()

summary = pd.DataFrame(
    [
        {
            "n_frozen_drivers": n,
            "n_testable": n,
            "n_direction_concordant": k,
            "direction_concordance_fraction":
                k / n,
            "direction_binomial_p":
                binom.pvalue,
            "spearman_rho":
                rho.statistic,
            "spearman_p":
                rho.pvalue,
            "concordant_FDR":
                int(
                    vc.get(
                        "PROTEOMIC_REGULATORY_CONCORDANT_FDR",
                        0
                    )
                ),
            "concordant_nominal":
                int(
                    vc.get(
                        "PROTEOMIC_REGULATORY_CONCORDANT_NOMINAL",
                        0
                    )
                ),
            "discordant_FDR":
                int(
                    vc.get(
                        "PROTEOMIC_REGULATORY_DISCORDANT_FDR",
                        0
                    )
                ),
            "discordant_nominal":
                int(
                    vc.get(
                        "PROTEOMIC_REGULATORY_DISCORDANT_NOMINAL",
                        0
                    )
                ),
        }
    ]
)

# ============================================================
# 7. NON-ICU SENSITIVITY
# ============================================================

nonicu_rows = []

nonicu_universe = np.asarray(
    nonicu_gene["protein_median_t"],
    dtype=float
)

# Separate deterministic RNG stream.
rng = np.random.default_rng(
    SEED + 1
)

for _, d in drivers.iterrows():

    tf = d["TF"]

    z = net.loc[
        net["source"] == tf,
        ["target", "weight"]
    ].copy()

    z = z.loc[
        z["target"].isin(nonicu_t)
    ].copy()

    if z["target"].nunique() < 10:
        continue

    target_stats = np.asarray(
        [
            nonicu_t[g]
            for g in z["target"]
        ],
        dtype=float
    )

    weights = np.asarray(
        z["weight"],
        dtype=float
    )

    score, p = permutation_test(
        target_stats,
        weights,
        nonicu_universe
    )

    nonicu_rows.append(
        {
            "TF": tf,
            "nonICU_program_score":
                score,
            "nonICU_permutation_P":
                p,
            "nonICU_direction_concordant":
                np.sign(score) ==
                np.sign(d["flu_delta"]),
        }
    )

nonicu_res = pd.DataFrame(
    nonicu_rows
)

result = result.merge(
    nonicu_res,
    on="TF",
    how="left",
    validate="one_to_one"
)

nk = int(
    result[
        "nonICU_direction_concordant"
    ].sum()
)

nbinom = binomtest(
    nk,
    len(result),
    p=0.5,
    alternative="greater"
)

nrho = spearmanr(
    result["flu_delta"],
    result["nonICU_program_score"]
)

sensitivity = pd.DataFrame(
    [
        {
            "analysis":
                "nonICU_frozen37",
            "n_drivers":
                len(result),
            "n_direction_concordant":
                nk,
            "direction_concordance_fraction":
                nk / len(result),
            "direction_binomial_p":
                nbinom.pvalue,
            "spearman_rho":
                nrho.statistic,
            "spearman_p":
                nrho.pvalue,
        }
    ]
)

# ============================================================
# 8. WRITE OUTPUTS
# ============================================================

result = result.sort_values(
    [
        "proteomic_FDR",
        "permutation_P",
        "TF",
    ]
)

result_file = OUT / (
    "Influenza-P01_FROZEN37_regulatory_program_validation_v1.0.tsv"
)

summary_file = OUT / (
    "Influenza-P01_FROZEN37_regulatory_program_summary_v1.0.tsv"
)

sensitivity_file = OUT / (
    "Influenza-P01_FROZEN37_regulatory_nonICU_sensitivity_v1.0.tsv"
)

result.to_csv(
    result_file,
    sep="\t",
    index=False
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False
)

sensitivity.to_csv(
    sensitivity_file,
    sep="\t",
    index=False
)

# ============================================================
# 9. COMPACT REPORT
# ============================================================

print("\n=== PRIMARY FROZEN37 REGULATORY VALIDATION ===")
print(
    summary.to_string(
        index=False
    )
)

print("\n=== VALIDATION CATEGORIES ===")
print(
    result[
        "validation_category"
    ]
    .value_counts()
    .to_string()
)

print("\n=== FDR-SUPPORTED CONCORDANT DRIVERS ===")

x = result.loc[
    result["validation_category"] ==
    "PROTEOMIC_REGULATORY_CONCORDANT_FDR",
    [
        "TF",
        "flu_delta",
        "proteomic_program_score",
        "n_measurable_targets",
        "permutation_P",
        "proteomic_FDR",
        "regulatory_class",
    ]
]

print(
    x.to_string(index=False)
    if len(x)
    else "None"
)

print("\n=== FDR-SUPPORTED DISCORDANT DRIVERS ===")

x = result.loc[
    result["validation_category"] ==
    "PROTEOMIC_REGULATORY_DISCORDANT_FDR",
    [
        "TF",
        "flu_delta",
        "proteomic_program_score",
        "n_measurable_targets",
        "permutation_P",
        "proteomic_FDR",
        "regulatory_class",
    ]
]

print(
    x.to_string(index=False)
    if len(x)
    else "None"
)

print("\n=== NON-ICU SENSITIVITY ===")
print(
    sensitivity.to_string(
        index=False
    )
)

print("\nWritten:")
print(result_file)
print(summary_file)
print(sensitivity_file)

print(
    "\nSTATUS: FROZEN37 REGULATORY PROGRAM VALIDATION COMPLETE"
)
