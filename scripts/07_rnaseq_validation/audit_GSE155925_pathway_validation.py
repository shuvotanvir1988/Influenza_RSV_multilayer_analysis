import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

INFILE = Path(
    "results/rnaseq_validation/validation/pathway_validation/"
    "GSE155925_preregistered_pathway_GSEA_M3.tsv"
)

OUTDIR = Path(
    "results/rnaseq_validation/validation/pathway_validation"
)

d = pd.read_csv(INFILE, sep="\t")

# ------------------------------------------------------------
# Avoid double-counting pathways assigned to >1 domain
# for GLOBAL pathway-level calculations.
# ------------------------------------------------------------

global_d = (
    d.sort_values(["pathway", "preregistered_domain"])
     .drop_duplicates("pathway")
     .copy()
)

global_d = global_d[global_d["NES"].notna()].copy()

print("========================================")
print("1. GLOBAL UNIQUE-PATHWAY AUDIT")
print("========================================")

print("Evaluable unique pathways:", len(global_d))

for micro_col, label in [
    ("NES_RSV_GPL10558", "GPL10558"),
    ("NES_RSV_GPL6884", "GPL6884"),
]:

    x = global_d[[micro_col, "NES"]].dropna()

    pr, pp = pearsonr(x[micro_col], x["NES"])
    sr, sp = spearmanr(x[micro_col], x["NES"])

    concord = (
        np.sign(x[micro_col]) ==
        np.sign(x["NES"])
    )

    print(f"\n{label}")
    print(
        "Direction concordance:",
        int(concord.sum()), "/", len(concord),
        f"({concord.mean():.3f})"
    )
    print(
        f"Pearson r={pr:.3f}, p={pp:.4g}"
    )
    print(
        f"Spearman rho={sr:.3f}, p={sp:.4g}"
    )

# ------------------------------------------------------------
# Original microarray agreement
# ------------------------------------------------------------

global_d["microarray_direction_agreement"] = (
    np.sign(global_d["NES_RSV_GPL10558"]) ==
    np.sign(global_d["NES_RSV_GPL6884"])
)

agree = global_d[
    global_d["microarray_direction_agreement"]
].copy()

agree["RNAseq_matches_microarray_direction"] = (
    np.sign(agree["NES"]) ==
    np.sign(agree["NES_RSV_GPL6884"])
)

print("\n========================================")
print("2. PATHWAYS WHERE MICROARRAYS AGREED")
print("========================================")

print(
    "Microarray-concordant pathways:",
    len(agree)
)

print(
    "RNA-seq same direction:",
    int(agree["RNAseq_matches_microarray_direction"].sum()),
    "/",
    len(agree),
    f"({agree['RNAseq_matches_microarray_direction'].mean():.3f})"
)

print(
    "RNA-seq significant AND same direction:",
    int(
        (
            agree["RNAseq_matches_microarray_direction"] &
            (agree["padj"] < 0.05)
        ).sum()
    )
)

# ------------------------------------------------------------
# Frozen RSV-replicated pathways
# ------------------------------------------------------------

rep = global_d[
    global_d["RSV_replicated"] == True
].copy()

rep["RNAseq_matches_GPL10558"] = (
    np.sign(rep["NES"]) ==
    np.sign(rep["NES_RSV_GPL10558"])
)

rep["RNAseq_matches_GPL6884"] = (
    np.sign(rep["NES"]) ==
    np.sign(rep["NES_RSV_GPL6884"])
)

rep["RNAseq_matches_both"] = (
    rep["RNAseq_matches_GPL10558"] &
    rep["RNAseq_matches_GPL6884"]
)

rep["RNAseq_FDR05"] = rep["padj"] < 0.05

print("\n========================================")
print("3. FROZEN RSV-REPLICATED PATHWAYS")
print("========================================")

print("Evaluable:", len(rep))

print(
    "Direction matches BOTH:",
    int(rep["RNAseq_matches_both"].sum()),
    "/",
    len(rep),
    f"({rep['RNAseq_matches_both'].mean():.3f})"
)

print(
    "Direction matches BOTH + RNAseq FDR<0.05:",
    int(
        (
            rep["RNAseq_matches_both"] &
            rep["RNAseq_FDR05"]
        ).sum()
    ),
    "/",
    len(rep)
)

# Correlation restricted to replicated pathways
if len(rep) >= 3:
    for micro_col, label in [
        ("NES_RSV_GPL10558", "GPL10558"),
        ("NES_RSV_GPL6884", "GPL6884"),
    ]:
        x = rep[[micro_col, "NES"]].dropna()

        if len(x) >= 3:
            pr, pp = pearsonr(
                x[micro_col], x["NES"]
            )
            sr, sp = spearmanr(
                x[micro_col], x["NES"]
            )

            print(
                f"{label}: Pearson={pr:.3f} "
                f"(p={pp:.4g}); "
                f"Spearman={sr:.3f} "
                f"(p={sp:.4g})"
            )

# ------------------------------------------------------------
# Domain-level audit
# Keep duplicated pathway-domain assignments because
# domain membership itself is the quantity being summarized.
# ------------------------------------------------------------

domain_rows = d[d["NES"].notna()].copy()

domain_rows["microarray_agreed"] = (
    np.sign(domain_rows["NES_RSV_GPL10558"]) ==
    np.sign(domain_rows["NES_RSV_GPL6884"])
)

domain_rows["RNAseq_matches_both"] = (
    (
        np.sign(domain_rows["NES"]) ==
        np.sign(domain_rows["NES_RSV_GPL10558"])
    )
    &
    (
        np.sign(domain_rows["NES"]) ==
        np.sign(domain_rows["NES_RSV_GPL6884"])
    )
)

domain_rows["RNAseq_FDR05"] = (
    domain_rows["padj"] < 0.05
)

domain_rows["RNAseq_FDR05_and_matches_both"] = (
    domain_rows["RNAseq_FDR05"] &
    domain_rows["RNAseq_matches_both"]
)

domain_rows["replicated_and_RNAseq_match"] = (
    (domain_rows["RSV_replicated"] == True) &
    domain_rows["RNAseq_matches_both"]
)

domain_rows["replicated_and_RNAseq_sig_match"] = (
    (domain_rows["RSV_replicated"] == True) &
    domain_rows["RNAseq_matches_both"] &
    domain_rows["RNAseq_FDR05"]
)

domain = (
    domain_rows
    .groupby("preregistered_domain")
    .agg(
        evaluable=("pathway", "size"),
        original_RSV_replicated=("RSV_replicated", "sum"),
        microarray_agreed=("microarray_agreed", "sum"),
        RNAseq_direction_matches_both=("RNAseq_matches_both", "sum"),
        RNAseq_FDR05=("RNAseq_FDR05", "sum"),
        RNAseq_FDR05_and_matches_both=(
            "RNAseq_FDR05_and_matches_both", "sum"
        ),
        replicated_RNAseq_direction_match=(
            "replicated_and_RNAseq_match", "sum"
        ),
        replicated_RNAseq_sig_direction_match=(
            "replicated_and_RNAseq_sig_match", "sum"
        ),
        median_RNAseq_NES=("NES", "median")
    )
    .reset_index()
)

domain["RNAseq_match_fraction"] = (
    domain["RNAseq_direction_matches_both"] /
    domain["evaluable"]
)

print("\n========================================")
print("4. DOMAIN-LEVEL AUDIT")
print("========================================")

print(domain.to_string(index=False))

# ------------------------------------------------------------
# Interferon detail
# ------------------------------------------------------------

interferon = domain_rows[
    domain_rows["preregistered_domain"] ==
    "Interferon"
].copy()

print("\n========================================")
print("5. INTERFERON DETAIL")
print("========================================")

cols = [
    "pathway",
    "NES_RSV_GPL10558",
    "padj_RSV_GPL10558",
    "NES_RSV_GPL6884",
    "padj_RSV_GPL6884",
    "RSV_replicated",
    "NES",
    "padj",
    "RNAseq_matches_both"
]

print(
    interferon[cols]
    .sort_values("padj")
    .to_string(index=False)
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

rep.to_csv(
    OUTDIR /
    "GSE155925_RSV_replicated_pathway_validation.tsv",
    sep="\t",
    index=False
)

domain.to_csv(
    OUTDIR /
    "GSE155925_domain_validation_summary.tsv",
    sep="\t",
    index=False
)

interferon.to_csv(
    OUTDIR /
    "GSE155925_interferon_validation_detail.tsv",
    sep="\t",
    index=False
)

agree.to_csv(
    OUTDIR /
    "GSE155925_microarray_concordant_pathway_validation.tsv",
    sep="\t",
    index=False
)

print("\nAudit complete.")
