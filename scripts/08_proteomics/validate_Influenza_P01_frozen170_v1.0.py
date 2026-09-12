from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

ROOT = Path("results/proteomics_validation/Influenza-P01")

TARGET = Path(
    "results/proteomics_validation/frozen_targets/"
    "PROTEOMICS_TARGETS_170_FROZEN_v1.0.tsv"
)

MAP = ROOT / (
    "frozen_mapping/"
    "Influenza-P01_FROZEN_ASSAY_MAPPING_v1.0.tsv"
)

PRIMARY = ROOT / (
    "differential_proteomics/"
    "Influenza-P01_primary_adjusted_infection_limma_v1.0.tsv"
)

UNADJ = ROOT / (
    "differential_proteomics/"
    "Influenza-P01_unadjusted_infection_limma_v1.0.tsv"
)

NONICU = ROOT / (
    "differential_proteomics/"
    "Influenza-P01_nonICU_infection_limma_v1.0.tsv"
)

OUT = ROOT / "validation"
OUT.mkdir(parents=True, exist_ok=True)

target = pd.read_csv(TARGET, sep="\t")
mapping = pd.read_csv(MAP, sep="\t")
primary = pd.read_csv(PRIMARY, sep="\t")
unadj = pd.read_csv(UNADJ, sep="\t")
nonicu = pd.read_csv(NONICU, sep="\t")

for x in [mapping, primary, unadj, nonicu]:
    x["row_index"] = pd.to_numeric(
        x["row_index"],
        errors="coerce"
    ).astype("Int64")

meas = mapping.loc[
    mapping["measurable"].astype(str).str.lower().isin(
        ["true", "1"]
    )
].copy()

assert meas["gene_symbol"].nunique() == 84

tx = target[
    [
        "gene_symbol",
        "proteomics_architecture_class",
        "flu_logFC",
        "flu_FDR",
    ]
].copy()

tx = tx.rename(
    columns={
        "proteomics_architecture_class": "architecture_class",
        "flu_logFC": "transcript_logFC",
        "flu_FDR": "transcript_FDR",
    }
)

assay = meas.merge(
    primary[
        [
            "row_index",
            "logFC",
            "P.Value",
            "adj.P.Val",
            "standardized_effect",
            "residual_sigma",
        ]
    ],
    on="row_index",
    how="left",
    validate="many_to_one"
)

assay = assay.merge(
    tx,
    on=["gene_symbol", "architecture_class"],
    how="left",
    validate="many_to_one"
)

assert assay["logFC"].notna().all()
assert assay["transcript_logFC"].notna().all()

assay["transcript_direction"] = np.sign(
    assay["transcript_logFC"]
)

assay["protein_direction"] = np.sign(
    assay["standardized_effect"]
)

assay["assay_direction_concordant"] = (
    assay["transcript_direction"] ==
    assay["protein_direction"]
)

gene_rows = []

for gene, z in assay.groupby("gene_symbol"):

    z = z.copy()

    tx_effect = z["transcript_logFC"].iloc[0]
    cls = z["architecture_class"].iloc[0]

    median_std = z["standardized_effect"].median()
    median_logfc = z["logFC"].median()

    tx_dir = np.sign(tx_effect)
    protein_dir = np.sign(median_std)

    gene_concordant = bool(tx_dir == protein_dir)

    unanimous_concordant = bool(
        z["assay_direction_concordant"].all()
    )

    all_fdr = bool(
        (z["adj.P.Val"] < 0.05).all()
    )

    all_nominal = bool(
        (z["P.Value"] < 0.05).all()
    )

    all_opposite = bool(
        (~z["assay_direction_concordant"]).all()
    )

    if unanimous_concordant and all_fdr:
        category = "PROTEIN_CONCORDANT_FDR"
    elif unanimous_concordant and all_nominal:
        category = "PROTEIN_CONCORDANT_NOMINAL"
    elif gene_concordant:
        category = "PROTEIN_CONCORDANT_UNSUPPORTED"
    elif all_opposite and all_fdr:
        category = "PROTEIN_DISCORDANT_FDR"
    elif all_opposite and all_nominal:
        category = "PROTEIN_DISCORDANT_NOMINAL"
    else:
        category = "PROTEIN_DISCORDANT_UNSUPPORTED"

    gene_rows.append({
        "gene_symbol": gene,
        "architecture_class": cls,
        "mapping_method": z["mapping_method"].iloc[0],
        "n_assays": len(z),
        "transcript_logFC": tx_effect,
        "protein_median_logFC": median_logfc,
        "protein_median_standardized_effect": median_std,
        "gene_direction_concordant": gene_concordant,
        "all_assays_direction_concordant": unanimous_concordant,
        "fraction_assays_direction_concordant":
            z["assay_direction_concordant"].mean(),
        "all_assays_FDR05": all_fdr,
        "all_assays_P05": all_nominal,
        "min_assay_FDR": z["adj.P.Val"].min(),
        "max_assay_FDR": z["adj.P.Val"].max(),
        "validation_category": category,
    })

gene = pd.DataFrame(gene_rows)

assert len(gene) == 84

def summarize(z, label):
    n = len(z)
    k = int(z["gene_direction_concordant"].sum())

    b = binomtest(
        k,
        n,
        p=0.5,
        alternative="greater"
    )

    rho, rho_p = spearmanr(
        z["transcript_logFC"],
        z["protein_median_standardized_effect"]
    )

    vc = z["validation_category"].value_counts()

    return {
        "group": label,
        "n_genes": n,
        "n_direction_concordant": k,
        "direction_concordance_fraction": k / n,
        "direction_binomial_p": b.pvalue,
        "spearman_rho": rho,
        "spearman_p": rho_p,
        "concordant_FDR":
            int(vc.get("PROTEIN_CONCORDANT_FDR", 0)),
        "concordant_nominal":
            int(vc.get("PROTEIN_CONCORDANT_NOMINAL", 0)),
        "discordant_FDR":
            int(vc.get("PROTEIN_DISCORDANT_FDR", 0)),
        "discordant_nominal":
            int(vc.get("PROTEIN_DISCORDANT_NOMINAL", 0)),
    }

summaries = [
    summarize(gene, "ALL_84")
]

for cls in [
    "shared_core",
    "influenza_amplified",
    "RSV_amplified",
]:
    summaries.append(
        summarize(
            gene.loc[
                gene["architecture_class"] == cls
            ],
            cls
        )
    )

summaries.append(
    summarize(
        gene.loc[
            gene["mapping_method"] == "exact_gene_symbol"
        ],
        "EXACT_SYMBOL_ONLY"
    )
)

summary = pd.DataFrame(summaries)

# adjusted vs unadjusted
sens = primary[
    ["row_index", "logFC"]
].rename(
    columns={"logFC": "adjusted_logFC"}
).merge(
    unadj[
        ["row_index", "logFC"]
    ].rename(
        columns={"logFC": "unadjusted_logFC"}
    ),
    on="row_index"
)

sens_rho = spearmanr(
    sens["adjusted_logFC"],
    sens["unadjusted_logFC"]
)

# non-ICU frozen-gene sensitivity
ni = meas.merge(
    nonicu[
        [
            "row_index",
            "standardized_effect",
        ]
    ],
    on="row_index",
    how="left"
).merge(
    tx[
        [
            "gene_symbol",
            "architecture_class",
            "transcript_logFC",
        ]
    ],
    on=["gene_symbol", "architecture_class"],
    how="left"
)

ni_gene = (
    ni.groupby(
        [
            "gene_symbol",
            "architecture_class",
            "transcript_logFC",
        ],
        as_index=False
    )
    .agg(
        protein_median_standardized_effect=(
            "standardized_effect",
            "median"
        )
    )
)

ni_gene["direction_concordant"] = (
    np.sign(ni_gene["transcript_logFC"]) ==
    np.sign(
        ni_gene["protein_median_standardized_effect"]
    )
)

ni_k = int(
    ni_gene["direction_concordant"].sum()
)

ni_n = len(ni_gene)

ni_binom = binomtest(
    ni_k,
    ni_n,
    p=0.5,
    alternative="greater"
)

ni_rho = spearmanr(
    ni_gene["transcript_logFC"],
    ni_gene["protein_median_standardized_effect"]
)

sensitivity = pd.DataFrame([
    {
        "analysis": "adjusted_vs_unadjusted_all_assays",
        "n": len(sens),
        "spearman_rho": sens_rho.statistic,
        "p_value": sens_rho.pvalue,
    },
    {
        "analysis": "nonICU_frozen84",
        "n": ni_n,
        "spearman_rho": ni_rho.statistic,
        "p_value": ni_rho.pvalue,
        "n_direction_concordant": ni_k,
        "direction_concordance_fraction": ni_k / ni_n,
        "direction_binomial_p": ni_binom.pvalue,
    }
])

assay_file = OUT / (
    "Influenza-P01_frozen84_assay_level_validation_v1.0.tsv"
)

gene_file = OUT / (
    "Influenza-P01_frozen84_gene_level_validation_v1.0.tsv"
)

summary_file = OUT / (
    "Influenza-P01_frozen84_validation_summary_v1.0.tsv"
)

sensitivity_file = OUT / (
    "Influenza-P01_validation_sensitivity_summary_v1.0.tsv"
)

assay.to_csv(
    assay_file,
    sep="\t",
    index=False
)

gene.to_csv(
    gene_file,
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

print("=== PRIMARY FROZEN PROTEOMIC VALIDATION ===")
print(summary.to_string(index=False))

print("\n=== VALIDATION CATEGORY COUNTS ===")
print(
    gene["validation_category"]
    .value_counts()
    .to_string()
)

print("\n=== SENSITIVITY ===")
print(
    sensitivity.to_string(index=False)
)

print("\n=== FDR-SUPPORTED CONCORDANT GENES ===")

x = gene.loc[
    gene["validation_category"] ==
    "PROTEIN_CONCORDANT_FDR"
].sort_values(
    "protein_median_standardized_effect",
    key=lambda s: s.abs(),
    ascending=False
)

if len(x):
    print(
        x[
            [
                "gene_symbol",
                "architecture_class",
                "transcript_logFC",
                "protein_median_standardized_effect",
                "n_assays",
            ]
        ].to_string(index=False)
    )
else:
    print("None")

print("\n=== FDR-SUPPORTED DISCORDANT GENES ===")

x = gene.loc[
    gene["validation_category"] ==
    "PROTEIN_DISCORDANT_FDR"
]

if len(x):
    print(
        x[
            [
                "gene_symbol",
                "architecture_class",
                "transcript_logFC",
                "protein_median_standardized_effect",
                "n_assays",
            ]
        ].to_string(index=False)
    )
else:
    print("None")

print("\nSTATUS: FROZEN PROTEOMIC VALIDATION COMPLETE")
