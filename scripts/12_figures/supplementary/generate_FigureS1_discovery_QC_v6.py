#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path.home() / "Influenza_RSV_Project"

GENE = ROOT / "results/DS001_GSE38900/gene_level/tables"
META = (
    ROOT / "results/DS001_GSE38900/differential_expression/tables/"
    "DS001_GSE38900_expression_metadata_sample_mapping.tsv"
)
OUT = (
    ROOT / "results/supplementary_figures/"
    "FigureS1_discovery_QC/figures"
)
TAB = (
    ROOT / "results/supplementary_figures/"
    "FigureS1_discovery_QC/tables"
)

OUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Session 42 — Supplementary Figure S1 v6
# Presentation-only refinement of frozen Figure S1 v4 / Session 42 v5.
#
# Scientific content is unchanged:
#   - same upstream PCA coordinates
#   - same sample-level expression statistics
#   - same group assignments
#   - same variance percentages
#
# v6 presentation changes:
#   - more horizontal and vertical space between panels
#   - panel labels aligned cleanly with panel-title baselines
#   - readability-first typography
#   - PDF + editable SVG + 600-dpi PNG
# ---------------------------------------------------------------------

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 18,
    "font.weight": "bold",
    "axes.titlesize": 23,
    "axes.titleweight": "bold",
    "axes.labelsize": 20,
    "axes.labelweight": "bold",
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "legend.fontsize": 17,
    "axes.linewidth": 1.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

COLORS = {
    "Healthy control": "#6E6E6E",
    "Influenza A acute": "#E36C0A",
    "RSV acute": "#149C7E",
    "RSV recovery": "#6F6AB1",
    "HRV acute": "#B8B8B8",
}

MARKERS = {
    "Healthy control": "o",
    "Influenza A acute": "^",
    "RSV acute": "s",
    "RSV recovery": "D",
    "HRV acute": "v",
}

meta = pd.read_csv(META, sep="\t")


def load_platform(platform):
    pca = pd.read_csv(
        GENE / f"{platform}_postcollapse_PCA_coordinates.tsv",
        sep="\t",
    )
    stats = pd.read_csv(
        GENE / f"{platform}_postcollapse_sample_statistics.tsv",
        sep="\t",
    )

    m = meta.loc[
        meta["platform_id"].eq(platform),
        ["expression_sample_id", "harmonized_group"],
    ].drop_duplicates()

    pca = pca.merge(
        m,
        left_on="sample",
        right_on="expression_sample_id",
        how="left",
        validate="one_to_one",
    )

    stats = stats.merge(
        m,
        left_on="sample",
        right_on="expression_sample_id",
        how="left",
        validate="one_to_one",
    )

    if pca["harmonized_group"].isna().any():
        missing = pca.loc[
            pca["harmonized_group"].isna(), "sample"
        ].tolist()
        raise RuntimeError(
            f"Missing PCA metadata for {platform}: {missing}"
        )

    if stats["harmonized_group"].isna().any():
        missing = stats.loc[
            stats["harmonized_group"].isna(), "sample"
        ].tolist()
        raise RuntimeError(
            f"Missing QC metadata for {platform}: {missing}"
        )

    pca.to_csv(
        TAB / f"FigureS1_{platform}_PCA_source.tsv",
        sep="\t",
        index=False,
    )
    stats.to_csv(
        TAB / f"FigureS1_{platform}_sample_statistics_source.tsv",
        sep="\t",
        index=False,
    )

    return pca, stats


g6884_pca, g6884_stats = load_platform("GPL6884")
g10558_pca, g10558_stats = load_platform("GPL10558")

summary = pd.read_csv(
    GENE / "DS001_GSE38900_postcollapse_qc_summary.tsv",
    sep="\t",
).set_index("platform")


# Larger canvas with deliberately increased inter-panel spacing.
fig = plt.figure(figsize=(19.5, 14.0))

gs = fig.add_gridspec(
    2,
    2,
    left=0.078,
    right=0.975,
    bottom=0.155,
    top=0.945,
    wspace=0.42,   # increased from v5
    hspace=0.50,   # increased from v5
)

axes = [
    fig.add_subplot(gs[0, 0]),
    fig.add_subplot(gs[0, 1]),
    fig.add_subplot(gs[1, 0]),
    fig.add_subplot(gs[1, 1]),
]


def style_axes(ax):
    ax.tick_params(
        axis="both",
        which="major",
        width=1.6,
        length=6,
        labelsize=17,
    )
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["bottom"].set_linewidth(1.8)


def plot_pca(ax, df, platform, title):
    groups = [
        g for g in COLORS
        if g in set(df["harmonized_group"])
    ]

    for group in groups:
        d = df[df["harmonized_group"] == group]
        ax.scatter(
            d["PC1"],
            d["PC2"],
            s=78,
            marker=MARKERS[group],
            c=COLORS[group],
            edgecolors="white",
            linewidths=0.9,
            alpha=0.90,
            zorder=3,
        )

    pc1 = summary.loc[platform, "PC1_variance_percent"]
    pc2 = summary.loc[platform, "PC2_variance_percent"]

    ax.set_xlabel(
        f"PC1 ({pc1:.1f}% variance)",
        labelpad=10,
    )
    ax.set_ylabel(
        f"PC2 ({pc2:.1f}% variance)",
        labelpad=10,
    )
    ax.set_title(
        title,
        loc="left",
        pad=16,
    )

    ax.axhline(
        0,
        lw=1.15,
        color="#9A9A9A",
        alpha=0.32,
        zorder=1,
    )
    ax.axvline(
        0,
        lw=1.15,
        color="#9A9A9A",
        alpha=0.32,
        zorder=1,
    )

    style_axes(ax)


def plot_medians(ax, df, title):
    groups = [
        g for g in COLORS
        if g in set(df["harmonized_group"])
    ]

    positions = {g: i for i, g in enumerate(groups)}

    for group in groups:
        d = df[df["harmonized_group"] == group].copy()
        x0 = positions[group]

        n = len(d)
        if n == 1:
            offsets = [0.0]
        else:
            offsets = [
                -0.20 + 0.40 * i / (n - 1)
                for i in range(n)
            ]

        ax.scatter(
            [x0 + o for o in offsets],
            d["median"],
            s=54,
            marker=MARKERS[group],
            c=COLORS[group],
            edgecolors="white",
            linewidths=0.75,
            alpha=0.88,
            zorder=3,
        )

        med = d["median"].median()
        ax.plot(
            [x0 - 0.27, x0 + 0.27],
            [med, med],
            color="black",
            lw=2.8,
            solid_capstyle="round",
            zorder=4,
        )

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(
        [g.replace(" ", "\n", 1) for g in groups],
        fontweight="bold",
        linespacing=1.05,
    )
    ax.set_ylabel(
        "Sample median expression",
        labelpad=10,
    )
    ax.set_title(
        title,
        loc="left",
        pad=16,
    )

    style_axes(ax)


plot_pca(
    axes[0],
    g6884_pca,
    "GPL6884",
    "GPL6884 post-collapse sample structure",
)

plot_medians(
    axes[1],
    g6884_stats,
    "GPL6884 sample-level expression QC",
)

plot_pca(
    axes[2],
    g10558_pca,
    "GPL10558",
    "GPL10558 post-collapse sample structure",
)

plot_medians(
    axes[3],
    g10558_stats,
    "GPL10558 sample-level expression QC",
)


# ---------------------------------------------------------------------
# Panel labels in figure coordinates.
# This avoids the "floating" label appearance seen in v5 and keeps
# A/B/C/D consistently aligned with the title baseline.
# ---------------------------------------------------------------------
fig.canvas.draw()

for ax, label in zip(axes, "ABCD"):
    bbox = ax.get_position()

    fig.text(
        bbox.x0 - 0.035,
        bbox.y1 + 0.014,
        label,
        fontsize=25,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


# Shared legend.
legend_groups = [
    "Healthy control",
    "Influenza A acute",
    "RSV acute",
    "RSV recovery",
    "HRV acute",
]

legend_handles = [
    Line2D(
        [0],
        [0],
        marker=MARKERS[g],
        linestyle="None",
        markerfacecolor=COLORS[g],
        markeredgecolor="white",
        markeredgewidth=0.9,
        markersize=11,
        label=g,
    )
    for g in legend_groups
]

legend = fig.legend(
    handles=legend_handles,
    labels=legend_groups,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.038),
    ncol=5,
    frameon=False,
    fontsize=17,
    handletextpad=0.55,
    columnspacing=1.70,
)

for text in legend.get_texts():
    text.set_fontweight("bold")


pdf = OUT / "FigureS1_v6_discovery_postcollapse_QC.pdf"
png = OUT / "FigureS1_v6_discovery_postcollapse_QC.png"
svg = OUT / "FigureS1_v6_discovery_postcollapse_QC.svg"

fig.savefig(
    pdf,
    bbox_inches="tight",
    facecolor="white",
)

fig.savefig(
    svg,
    bbox_inches="tight",
    facecolor="white",
)

fig.savefig(
    png,
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
)

plt.close(fig)

print("=== FIGURE S1 v6 GENERATION COMPLETE ===")
print(pdf)
print(svg)
print(png)

print("\n=== PANEL SAMPLE COUNTS ===")
for platform, pca in [
    ("GPL6884", g6884_pca),
    ("GPL10558", g10558_pca),
]:
    print(platform)
    print(
        pca["harmonized_group"]
        .value_counts()
        .to_string()
    )

print("\n=== QC SUMMARY ===")
print(
    summary[
        [
            "gene_count",
            "sample_count",
            "missing_values",
            "mean_pairwise_sample_correlation",
            "min_pairwise_sample_correlation",
            "PC1_variance_percent",
            "PC2_variance_percent",
        ]
    ].to_string()
)

print(
    "\nNo PCA coordinates, expression statistics, "
    "sample assignments, or upstream QC values were recomputed."
)
print(
    "Session 42 v6 is a presentation-only refinement "
    "with increased inter-panel spacing and aligned panel labels."
)
