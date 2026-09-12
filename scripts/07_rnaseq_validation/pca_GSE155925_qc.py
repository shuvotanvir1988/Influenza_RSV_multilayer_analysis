import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

COUNTFILE = Path(
    "data/rnaseq_validation/GSE155925/raw_counts/"
    "GSE155925_Raw_counts_matrix.txt.gz"
)

METAFILE = Path(
    "results/rnaseq_validation/tables/"
    "GSE155925_sample_eligibility_FROZEN_v1.0.tsv"
)

OUTDIR = Path("results/rnaseq_validation/qc/GSE155925")
OUTDIR.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv(COUNTFILE, sep="\t")
counts = counts.rename(columns={counts.columns[0]: "gene_id"})
counts = counts.set_index("gene_id")

meta = pd.read_csv(METAFILE, sep="\t")
meta = meta[meta["primary_eligible"] == True].copy()

samples = meta["count_matrix_label"].tolist()
x = counts[samples]

# QC-only expression filter:
# count >=10 in at least 5 eligible samples
keep = (x >= 10).sum(axis=1) >= 5
x = x.loc[keep]

print("Genes used for QC PCA:", x.shape[0])
print("Samples:", x.shape[1])

# Counts per million
lib = x.sum(axis=0)
cpm = x.div(lib, axis=1) * 1e6

# Log CPM
logcpm = np.log2(cpm + 1)

# PCA on samples.
# Center genes but do NOT variance-scale each gene.
X = logcpm.T
X_centered = X - X.mean(axis=0)

pca = PCA(n_components=10)
pcs = pca.fit_transform(X_centered)

pc = pd.DataFrame(
    pcs,
    index=X.index,
    columns=[f"PC{i+1}" for i in range(10)]
)

pc["count_matrix_label"] = pc.index

pc = pc.merge(
    meta[
        [
            "count_matrix_label",
            "geo_accession",
            "primary_group",
            "age_months",
            "sex",
            "hospital_batch",
            "enrollment_year_batch",
        ]
    ],
    on="count_matrix_label",
    how="left"
)

pc.to_csv(
    OUTDIR / "GSE155925_PCA_coordinates.tsv",
    sep="\t",
    index=False
)

variance = pd.DataFrame({
    "PC": [f"PC{i+1}" for i in range(10)],
    "variance_fraction": pca.explained_variance_ratio_
})

variance.to_csv(
    OUTDIR / "GSE155925_PCA_variance.tsv",
    sep="\t",
    index=False
)

print("\n=== PCA VARIANCE ===")
print(variance.to_string(index=False))

print("\n=== PCA COORDINATES ===")
print(
    pc[
        [
            "count_matrix_label",
            "primary_group",
            "hospital_batch",
            "enrollment_year_batch",
            "PC1",
            "PC2",
            "PC3",
        ]
    ].to_string(index=False)
)

# Plot phenotype
fig, ax = plt.subplots(figsize=(7, 6))

for group, sub in pc.groupby("primary_group"):
    ax.scatter(
        sub["PC1"],
        sub["PC2"],
        label=group,
        s=55,
        alpha=0.8
    )

for _, r in pc.iterrows():
    ax.annotate(
        r["count_matrix_label"].replace("Case ", ""),
        (r["PC1"], r["PC2"]),
        fontsize=6,
        alpha=0.7
    )

ax.set_xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
)
ax.set_ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
)
ax.set_title("GSE155925 — eligible-sample RNA-seq QC")
ax.legend(frameon=False)

fig.tight_layout()
fig.savefig(
    OUTDIR / "GSE155925_PCA_primary_group.pdf"
)
fig.savefig(
    OUTDIR / "GSE155925_PCA_primary_group.png",
    dpi=300
)

plt.close(fig)

print("\nPCA files written to:", OUTDIR)
