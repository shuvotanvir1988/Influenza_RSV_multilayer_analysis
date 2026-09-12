import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path.home() / "Influenza_RSV_Project"

INPUT_DIR = ROOT / "data/processed/DS001_GSE38900/analysis_ready"
FIG_DIR = ROOT / "results/DS001_GSE38900/gene_level/figures/postcollapse_qc"
TABLE_DIR = ROOT / "results/DS001_GSE38900/gene_level/tables"

FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

platforms = {
    "GPL10558": INPUT_DIR / "GPL10558_gene_expression.tsv.gz",
    "GPL6884": INPUT_DIR / "GPL6884_gene_expression.tsv.gz",
}

summary_rows = []

for platform, path in platforms.items():

    print("\n" + "=" * 80)
    print(f"Post-collapse QC: {platform}")
    print("=" * 80)

    df = pd.read_csv(path, sep="\t")

    id_cols = ["gene_symbol", "entrez_gene_id", "probe_id"]
    sample_cols = [c for c in df.columns if c not in id_cols]

    X = df[sample_cols].apply(pd.to_numeric, errors="coerce")

    print(f"Genes: {X.shape[0]:,}")
    print(f"Samples: {X.shape[1]:,}")
    print(f"Missing values: {int(X.isna().sum().sum()):,}")

    # ------------------------------------------------------
    # Sample summary statistics
    # ------------------------------------------------------
    sample_stats = pd.DataFrame({
        "sample": sample_cols,
        "mean": X.mean(axis=0).values,
        "median": X.median(axis=0).values,
        "sd": X.std(axis=0, ddof=1).values,
        "min": X.min(axis=0).values,
        "max": X.max(axis=0).values,
    })

    sample_stats["platform"] = platform

    sample_stats.to_csv(
        TABLE_DIR / f"{platform}_postcollapse_sample_statistics.tsv",
        sep="\t",
        index=False,
    )

    # ------------------------------------------------------
    # BOXPLOT
    # ------------------------------------------------------
    plt.figure(figsize=(max(12, len(sample_cols) * 0.12), 7))

    plt.boxplot(
        [X[c].dropna().values for c in sample_cols],
        showfliers=False
    )

    plt.xticks(
        range(1, len(sample_cols) + 1),
        sample_cols,
        rotation=90,
        fontsize=6
    )

    plt.ylabel("VSN-normalized expression")
    plt.title(f"{platform} Gene-Level Expression Boxplot")
    plt.tight_layout()

    plt.savefig(
        FIG_DIR / f"{platform}_gene_level_boxplot.png",
        dpi=300
    )
    plt.close()

    # ------------------------------------------------------
    # DENSITY PLOT
    # ------------------------------------------------------
    plt.figure(figsize=(10, 7))

    for sample in sample_cols:
        values = X[sample].dropna()
        values.plot(
            kind="density",
            linewidth=0.7,
            alpha=0.45
        )

    plt.xlabel("VSN-normalized expression")
    plt.ylabel("Density")
    plt.title(f"{platform} Gene-Level Expression Density")
    plt.tight_layout()

    plt.savefig(
        FIG_DIR / f"{platform}_gene_level_density.png",
        dpi=300
    )
    plt.close()

    # ------------------------------------------------------
    # PCA
    # samples × genes
    # ------------------------------------------------------
    pca_input = X.T.copy()

    # Remove zero-variance genes
    gene_sd = pca_input.std(axis=0)
    pca_input = pca_input.loc[:, gene_sd > 0]

    # Median-imputation safeguard if any missing values
    if pca_input.isna().any().any():
        pca_input = pca_input.fillna(
            pca_input.median(axis=0)
        )

    scaled = StandardScaler().fit_transform(pca_input)

    pca = PCA(n_components=min(5, scaled.shape[0], scaled.shape[1]))
    pcs = pca.fit_transform(scaled)

    pca_df = pd.DataFrame({
        "sample": sample_cols,
        "PC1": pcs[:, 0],
        "PC2": pcs[:, 1],
    })

    pca_df["platform"] = platform

    pca_df.to_csv(
        TABLE_DIR / f"{platform}_postcollapse_PCA_coordinates.tsv",
        sep="\t",
        index=False,
    )

    plt.figure(figsize=(8, 7))

    plt.scatter(
        pca_df["PC1"],
        pca_df["PC2"],
        s=35
    )

    plt.xlabel(
        f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% variance)"
    )
    plt.ylabel(
        f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% variance)"
    )
    plt.title(f"{platform} Gene-Level PCA")

    plt.tight_layout()

    plt.savefig(
        FIG_DIR / f"{platform}_gene_level_PCA.png",
        dpi=300
    )
    plt.close()

    # ------------------------------------------------------
    # SAMPLE CORRELATION HEATMAP
    # ------------------------------------------------------
    corr = X.corr(method="pearson")

    corr.to_csv(
        TABLE_DIR / f"{platform}_postcollapse_sample_correlation.tsv",
        sep="\t"
    )

    plt.figure(
        figsize=(
            max(8, len(sample_cols) * 0.08),
            max(8, len(sample_cols) * 0.08)
        )
    )

    plt.imshow(
        corr.values,
        aspect="auto",
        interpolation="nearest"
    )

    plt.colorbar(label="Pearson correlation")

    if len(sample_cols) <= 50:
        plt.xticks(
            range(len(sample_cols)),
            sample_cols,
            rotation=90,
            fontsize=6
        )
        plt.yticks(
            range(len(sample_cols)),
            sample_cols,
            fontsize=6
        )
    else:
        plt.xticks([])
        plt.yticks([])

    plt.title(f"{platform} Gene-Level Sample Correlation")
    plt.tight_layout()

    plt.savefig(
        FIG_DIR / f"{platform}_gene_level_correlation_heatmap.png",
        dpi=300
    )
    plt.close()

    # ------------------------------------------------------
    # Summary
    # ------------------------------------------------------
    corr_vals = corr.values[np.triu_indices_from(corr, k=1)]

    summary_rows.append({
        "platform": platform,
        "gene_count": X.shape[0],
        "sample_count": X.shape[1],
        "missing_values": int(X.isna().sum().sum()),
        "global_mean": float(X.values.mean()),
        "global_median": float(np.median(X.values)),
        "sample_mean_min": float(sample_stats["mean"].min()),
        "sample_mean_max": float(sample_stats["mean"].max()),
        "sample_median_min": float(sample_stats["median"].min()),
        "sample_median_max": float(sample_stats["median"].max()),
        "mean_pairwise_sample_correlation": float(np.mean(corr_vals)),
        "min_pairwise_sample_correlation": float(np.min(corr_vals)),
        "PC1_variance_percent": float(
            pca.explained_variance_ratio_[0] * 100
        ),
        "PC2_variance_percent": float(
            pca.explained_variance_ratio_[1] * 100
        ),
    })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    TABLE_DIR / "DS001_GSE38900_postcollapse_qc_summary.tsv",
    sep="\t",
    index=False
)

print("\n" + "=" * 80)
print("POST-COLLAPSE QC SUMMARY")
print("=" * 80)

print(summary.to_string(index=False))

print("\nPost-collapse QC completed successfully.")
