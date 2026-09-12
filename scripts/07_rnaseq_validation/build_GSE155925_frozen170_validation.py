import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr, binomtest

FROZEN = Path(
    "results/DS001_GSE38900/integrated_architecture/tables/"
    "Figure7A_high_confidence_gene_architecture.tsv"
)

RNASEQ = Path(
    "results/rnaseq_validation/GSE155925_DE/"
    "GSE155925_RSV_vs_negative_DESeq2_adjusted.tsv"
)

OUTDIR = Path("results/rnaseq_validation/validation")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Frozen Session 22 architecture
# ------------------------------------------------------------
frozen = pd.read_csv(FROZEN, sep="\t")

print("=== FROZEN ARCHITECTURE ===")
print("Rows:", len(frozen))
print(frozen["Figure7A_class"].value_counts())

assert len(frozen) == 170
assert frozen["gene_symbol"].is_unique

# ------------------------------------------------------------
# RNA-seq DE
# gene_id format: SYMBOL:ENSG...
# ------------------------------------------------------------
rna = pd.read_csv(RNASEQ, sep="\t")

rna["gene_symbol"] = (
    rna["gene_id"]
    .astype(str)
    .str.split(":", n=1)
    .str[0]
)

print("\n=== RNA-seq ===")
print("Rows:", len(rna))
print("Unique symbols:", rna["gene_symbol"].nunique())
print("Duplicated symbols:", rna["gene_symbol"].duplicated().sum())

# If duplicate symbols occur, retain row with smallest padj;
# if padj missing, use pvalue.
rna["_rank_p"] = rna["padj"].fillna(rna["pvalue"])

rna = (
    rna.sort_values(["gene_symbol", "_rank_p"])
       .drop_duplicates("gene_symbol", keep="first")
       .drop(columns="_rank_p")
)

# ------------------------------------------------------------
# Merge WITHOUT redefining frozen classes
# ------------------------------------------------------------
v = frozen.merge(
    rna[
        [
            "gene_symbol",
            "gene_id",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
        ]
    ],
    on="gene_symbol",
    how="left",
    validate="one_to_one"
)

v = v.rename(columns={
    "log2FoldChange": "rnaseq_RSV_log2FC",
    "lfcSE": "rnaseq_lfcSE",
    "stat": "rnaseq_stat",
    "pvalue": "rnaseq_pvalue",
    "padj": "rnaseq_padj",
    "baseMean": "rnaseq_baseMean",
    "gene_id": "rnaseq_gene_id",
})

v["rnaseq_evaluable"] = v["rnaseq_RSV_log2FC"].notna()

# Frozen RSV reference direction:
# GPL6884 is the larger discovery platform.
v["direction_concordant_GPL6884"] = np.where(
    v["rnaseq_evaluable"],
    np.sign(v["rnaseq_RSV_log2FC"]) == np.sign(v["rsv6884_logFC"]),
    np.nan
)

# Stronger criterion: RNA-seq agrees with BOTH frozen RSV platforms.
v["direction_concordant_both_RSV_arrays"] = np.where(
    v["rnaseq_evaluable"],
    (
        (np.sign(v["rnaseq_RSV_log2FC"]) == np.sign(v["rsv6884_logFC"])) &
        (np.sign(v["rnaseq_RSV_log2FC"]) == np.sign(v["rsv10558_logFC"]))
    ),
    np.nan
)

v["rnaseq_FDR05"] = (
    v["rnaseq_padj"].notna() &
    (v["rnaseq_padj"] < 0.05)
)

v["rnaseq_FDR05_absLFC05"] = (
    v["rnaseq_FDR05"] &
    (v["rnaseq_RSV_log2FC"].abs() >= 0.5)
)

v["rnaseq_FDR05_same_direction"] = (
    v["rnaseq_FDR05"] &
    (v["direction_concordant_both_RSV_arrays"] == True)
)

outfile = OUTDIR / "GSE155925_frozen170_validation_master.tsv"
v.to_csv(outfile, sep="\t", index=False)

# ------------------------------------------------------------
# Summary function
# ------------------------------------------------------------
def summarize(d, label):
    e = d[d["rnaseq_evaluable"]].copy()

    n_total = len(d)
    n_eval = len(e)

    n_conc = int(
        (e["direction_concordant_both_RSV_arrays"] == True).sum()
    )

    n_sig = int(e["rnaseq_FDR05"].sum())
    n_sig05 = int(e["rnaseq_FDR05_absLFC05"].sum())
    n_sig_conc = int(e["rnaseq_FDR05_same_direction"].sum())

    pearson = np.nan
    spearman = np.nan
    pearson_p = np.nan
    spearman_p = np.nan

    if n_eval >= 3:
        pearson, pearson_p = pearsonr(
            e["rsv6884_logFC"],
            e["rnaseq_RSV_log2FC"]
        )
        spearman, spearman_p = spearmanr(
            e["rsv6884_logFC"],
            e["rnaseq_RSV_log2FC"]
        )

    binom_p = (
        binomtest(n_conc, n_eval, 0.5, alternative="greater").pvalue
        if n_eval > 0 else np.nan
    )

    return {
        "class": label,
        "frozen_genes": n_total,
        "rnaseq_evaluable": n_eval,
        "direction_concordant": n_conc,
        "direction_concordance_fraction":
            n_conc / n_eval if n_eval else np.nan,
        "direction_binomial_p": binom_p,
        "rnaseq_FDR05": n_sig,
        "rnaseq_FDR05_fraction":
            n_sig / n_eval if n_eval else np.nan,
        "rnaseq_FDR05_absLFC05": n_sig05,
        "rnaseq_FDR05_same_direction": n_sig_conc,
        "pearson_r_GPL6884_vs_RNAseq": pearson,
        "pearson_p": pearson_p,
        "spearman_rho_GPL6884_vs_RNAseq": spearman,
        "spearman_p": spearman_p,
        "median_microarray_RSV_logFC":
            e["rsv6884_logFC"].median() if n_eval else np.nan,
        "median_RNAseq_RSV_log2FC":
            e["rnaseq_RSV_log2FC"].median() if n_eval else np.nan,
    }

summary = []

summary.append(summarize(v, "ALL_170"))

for cls, d in v.groupby("Figure7A_class", sort=False):
    summary.append(summarize(d, cls))

s = pd.DataFrame(summary)

summary_file = OUTDIR / "GSE155925_frozen170_validation_summary.tsv"
s.to_csv(summary_file, sep="\t", index=False)

# ------------------------------------------------------------
# Console audit
# ------------------------------------------------------------
print("\n=== VALIDATION COVERAGE ===")
print("Frozen genes:", len(v))
print("RNA-seq evaluable:", int(v["rnaseq_evaluable"].sum()))
print("Missing:", int((~v["rnaseq_evaluable"]).sum()))

print("\n=== VALIDATION SUMMARY ===")
print(s.to_string(index=False))

print("\n=== MISSING FROZEN GENES ===")
missing = v.loc[~v["rnaseq_evaluable"], ["gene_symbol", "Figure7A_class"]]
if len(missing):
    print(missing.to_string(index=False))
else:
    print("NONE")

print("\n=== DISCORDANT GENES ===")
discord = v.loc[
    (v["rnaseq_evaluable"]) &
    (v["direction_concordant_both_RSV_arrays"] == False),
    [
        "gene_symbol",
        "Figure7A_class",
        "rsv6884_logFC",
        "rsv10558_logFC",
        "rnaseq_RSV_log2FC",
        "rnaseq_padj",
    ]
].sort_values("rnaseq_padj")

print(discord.to_string(index=False))

print("\nWritten:")
print(outfile)
print(summary_file)
