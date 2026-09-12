import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

FROZEN = Path(
    "results/DS001_GSE38900/integrated_architecture/tables/"
    "Figure7A_high_confidence_gene_architecture.tsv"
)

DE_DIR = Path(
    "results/rnaseq_validation/GSE155925_DE/age_sensitivity"
)

OUT = Path(
    "results/rnaseq_validation/validation/"
    "GSE155925_frozen170_age_sensitivity.tsv"
)

f = pd.read_csv(FROZEN, sep="\t")

models = {
    "A_age30_allsex":
        "A_age30_allsex_RSV_vs_negative.tsv",

    "B_age30_maleonly":
        "B_age30_maleonly_RSV_vs_negative.tsv",
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

    d["_rank"] = d["padj"].fillna(d["pvalue"])

    d = (
        d.sort_values(["gene_symbol", "_rank"])
         .drop_duplicates("gene_symbol", keep="first")
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

    groups = [("ALL_170", z)] + list(
        z.groupby("Figure7A_class", sort=False)
    )

    for label, sub in groups:

        e = sub.dropna(
            subset=["log2FoldChange"]
        ).copy()

        concord = (
            np.sign(e["log2FoldChange"]) ==
            np.sign(e["rsv6884_logFC"])
        )

        if len(e) >= 3:
            pr = pearsonr(
                e["rsv6884_logFC"],
                e["log2FoldChange"]
            )

            sr = spearmanr(
                e["rsv6884_logFC"],
                e["log2FoldChange"]
            )
        else:
            pr = (np.nan, np.nan)
            sr = (np.nan, np.nan)

        rows.append({
            "model": model,
            "class": label,
            "evaluable": len(e),
            "direction_concordant": int(concord.sum()),
            "direction_concordance_fraction":
                float(concord.mean()),
            "pearson_r": pr[0],
            "pearson_p": pr[1],
            "spearman_rho": sr[0],
            "spearman_p": sr[1],
            "FDR05":
                int((e["padj"] < 0.05).sum()),
            "median_RNAseq_log2FC":
                e["log2FoldChange"].median()
        })

out = pd.DataFrame(rows)

out.to_csv(
    OUT,
    sep="\t",
    index=False
)

print(out.to_string(index=False))
print("\nWritten:", OUT)
