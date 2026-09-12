import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

FROZEN = Path(
    "results/DS001_GSE38900/integrated_architecture/tables/"
    "Figure7A_high_confidence_gene_architecture.tsv"
)

DE_DIR = Path(
    "results/rnaseq_validation/GSE155925_DE/model_sensitivity"
)

OUT = Path(
    "results/rnaseq_validation/validation/"
    "GSE155925_frozen170_model_sensitivity.tsv"
)

f = pd.read_csv(FROZEN, sep="\t")

models = {
    "M0_unadjusted":
        "GSE155925_M0_unadjusted_RSV_vs_negative.tsv",

    "M1_age":
        "GSE155925_M1_age_RSV_vs_negative.tsv",

    "M2_batch":
        "GSE155925_M2_batch_RSV_vs_negative.tsv",

    "M3_full":
        "GSE155925_M3_full_RSV_vs_negative.tsv",
}

rows = []

for model, filename in models.items():

    d = pd.read_csv(DE_DIR / filename, sep="\t")

    d["gene_symbol"] = (
        d["gene_id"]
        .astype(str)
        .str.split(":", n=1)
        .str[0]
    )

    d["_p"] = d["padj"].fillna(d["pvalue"])

    d = (
        d.sort_values(["gene_symbol", "_p"])
         .drop_duplicates("gene_symbol")
    )

    z = f.merge(
        d[
            [
                "gene_symbol",
                "log2FoldChange",
                "padj"
            ]
        ],
        on="gene_symbol",
        how="left"
    )

    for cls, sub in [
        ("ALL_170", z),
        *list(z.groupby("Figure7A_class"))
    ]:

        e = sub.dropna(subset=["log2FoldChange"]).copy()

        concord = (
            np.sign(e["log2FoldChange"]) ==
            np.sign(e["rsv6884_logFC"])
        )

        pr = pearsonr(
            e["rsv6884_logFC"],
            e["log2FoldChange"]
        ) if len(e) >= 3 else (np.nan, np.nan)

        sr = spearmanr(
            e["rsv6884_logFC"],
            e["log2FoldChange"]
        ) if len(e) >= 3 else (np.nan, np.nan)

        rows.append({
            "model": model,
            "class": cls,
            "evaluable": len(e),
            "direction_concordant": int(concord.sum()),
            "direction_concordance_fraction":
                float(concord.mean()),
            "pearson_r": pr[0],
            "pearson_p": pr[1],
            "spearman_rho": sr[0],
            "spearman_p": sr[1],
            "FDR05": int((e["padj"] < 0.05).sum()),
            "median_RNAseq_log2FC":
                e["log2FoldChange"].median()
        })

out = pd.DataFrame(rows)

out.to_csv(OUT, sep="\t", index=False)

print(out.to_string(index=False))
print("\nWritten:", OUT)
