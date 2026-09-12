#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ============================================================
# Paths
# ============================================================

ROOT = Path.home() / "Influenza_RSV_Project"

PBASE = (
    ROOT /
    "results/proteomics_validation/Influenza-P01"
)

OUT = (
    ROOT /
    "results/supplementary_figures/"
    "FigureS8_proteomic_validation"
)

FIGDIR = OUT / "figures"
TABDIR = OUT / "tables"

FIGDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)

PDF = FIGDIR / "FigureS8_v7_proteomic_validation.pdf"
PNG = FIGDIR / "FigureS8_v7_proteomic_validation.png"
SVG = FIGDIR / "FigureS8_v7_proteomic_validation.svg"


# ============================================================
# Authoritative frozen inputs
# ============================================================

COVERAGE_FILE = (
    PBASE /
    "frozen_mapping/"
    "Influenza-P01_FROZEN_MEASURABLE_COVERAGE_v1.0.tsv"
)

GENE_FILE = (
    PBASE /
    "validation/"
    "Influenza-P01_frozen84_gene_level_validation_v1.0.tsv"
)

DOMAIN_FILE = (
    PBASE /
    "pathway_validation/"
    "Influenza-P01_FROZEN_7_domain_validation_summary_v1.0.tsv"
)

CROSS_LAYER_FILE = (
    PBASE /
    "figure_source_data/"
    "Figure_PROTEOMICS_v3_cross_layer_summary.tsv"
)


# ============================================================
# Load
# ============================================================

coverage = pd.read_csv(
    COVERAGE_FILE,
    sep="\t",
)

genes = pd.read_csv(
    GENE_FILE,
    sep="\t",
)

domains = pd.read_csv(
    DOMAIN_FILE,
    sep="\t",
)

cross = pd.read_csv(
    CROSS_LAYER_FILE,
    sep="\t",
)


# ============================================================
# Hard QC
# ============================================================

if len(genes) != 84:
    raise RuntimeError(
        f"Expected 84 frozen measurable genes; found {len(genes)}."
    )

required_classes = {
    "shared_core",
    "influenza_amplified",
    "RSV_amplified",
}

observed_classes = set(
    coverage["architecture_class"].astype(str)
)

if not required_classes.issubset(observed_classes):
    raise RuntimeError(
        "Coverage table does not contain all three expected "
        "architecture classes."
    )

if len(domains) != 7:
    raise RuntimeError(
        f"Expected 7 preregistered domains; found {len(domains)}."
    )

if set(cross["validation_layer"]) != {
    "Gene",
    "Pathway",
    "Regulatory program",
}:
    raise RuntimeError(
        "Unexpected validation layers in frozen cross-layer summary."
    )


# ============================================================
# Figure-source tables
# ============================================================

class_order = [
    "shared_core",
    "influenza_amplified",
    "RSV_amplified",
]

A = (
    coverage
    .set_index("architecture_class")
    .loc[class_order]
    .reset_index()
    .copy()
)

A["n_unmeasured_genes"] = (
    A["n_frozen_genes"] -
    A["n_measurable_genes"]
)

A.to_csv(
    TABDIR /
    "FigureS8_panelA_frozen_architecture_coverage.tsv",
    sep="\t",
    index=False,
)


B = genes.copy()

B.to_csv(
    TABDIR /
    "FigureS8_panelB_frozen84_gene_validation.tsv",
    sep="\t",
    index=False,
)


domain_order = [
    "Interferon",
    "Innate_Cytokine",
    "Antigen_Presentation",
    "Neutrophil_Monocyte",
    "Adaptive_Immunity",
    "Metabolism",
    "Translation",
]

C = (
    domains
    .set_index("preregistered_domain")
    .loc[domain_order]
    .reset_index()
    .copy()
)

C.to_csv(
    TABDIR /
    "FigureS8_panelC_frozen7_domain_validation.tsv",
    sep="\t",
    index=False,
)


layer_order = [
    "Gene",
    "Pathway",
    "Regulatory program",
]

D = (
    cross
    .set_index("validation_layer")
    .loc[layer_order]
    .reset_index()
    .copy()
)

D.to_csv(
    TABDIR /
    "FigureS8_panelD_cross_layer_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# Style
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 18,
    "axes.titlesize": 23,
    "axes.titleweight": "bold",
    "axes.labelsize": 19,
    "axes.labelweight": "bold",
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def add_panel_header(ax, letter, title):
    """
    Explicit axis-coordinate positioning keeps all four
    panel labels and titles on the same baseline.
    """

    ax.text(
        -0.12,
        1.10,
        letter,
        transform=ax.transAxes,
        fontsize=26,
        fontweight=800,
        ha="left",
        va="bottom",
        clip_on=False,
    )

    ax.text(
        0.0,
        1.10,
        title,
        transform=ax.transAxes,
        fontsize=23,
        fontweight=800,
        ha="left",
        va="bottom",
        clip_on=False,
    )


# ============================================================
# Figure geometry
# ============================================================

fig = plt.figure(
    figsize=(22.0, 18.5)
)

gs = fig.add_gridspec(
    2,
    2,
    width_ratios=[0.92, 1.34],
    height_ratios=[1.0, 1.0],
    left=0.095,
    right=0.975,
    bottom=0.090,
    top=0.94,
    wspace=0.50,
    hspace=0.56,
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])


# ============================================================
# PANEL A — SomaScan coverage
# ============================================================

A_display = {
    "shared_core": "Shared core",
    "influenza_amplified": "Influenza-amplified",
    "RSV_amplified": "RSV-amplified",
}

yA = np.arange(len(A))

measured = A["n_measurable_genes"].to_numpy()
unmeasured = A["n_unmeasured_genes"].to_numpy()
total = A["n_frozen_genes"].to_numpy()

axA.barh(
    yA,
    measured,
    height=0.58,
    label="Measurable",
)

axA.barh(
    yA,
    unmeasured,
    left=measured,
    height=0.58,
    alpha=0.22,
    label="Not measurable",
)

for i in range(len(A)):
    axA.text(
        measured[i] / 2,
        i,
        f"{int(measured[i])}",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
    )

    axA.text(
        total[i] + 2,
        i,
        (
            f"{int(measured[i])}/{int(total[i])} "
            f"({100 * measured[i] / total[i]:.1f}%)"
        ),
        ha="left",
        va="center",
        fontsize=10.5,
    )

axA.set_yticks(yA)

axA.set_yticklabels(
    [
        A_display[x]
        for x in A["architecture_class"]
    ],
    fontweight="bold",
)

axA.invert_yaxis()

axA.set_xlabel(
    "Frozen genes"
)

axA.set_xlim(
    0,
    max(total) * 1.34,
)

overall_measured = int(
    A["n_measurable_genes"].sum()
)

overall_total = int(
    A["n_frozen_genes"].sum()
)

if overall_measured != 84 or overall_total != 170:
    raise RuntimeError(
        "Expected overall final frozen mapping 84/170; "
        f"found {overall_measured}/{overall_total}."
    )

axA.text(
    0.0,
    1.015,
    (
        f"Overall final mapping: "
        f"{overall_measured}/{overall_total} measurable "
        f"({100 * overall_measured / overall_total:.1f}%)"
    ),
    transform=axA.transAxes,
    fontsize=15.5,
    fontweight="bold",
    ha="left",
    va="bottom",
)

axA.legend(
    loc="lower left",
    bbox_to_anchor=(0.0, -0.22),
    ncol=2,
    frameon=False,
    borderaxespad=0.0,
    handletextpad=0.5,
    columnspacing=1.4,
)

add_panel_header(
    axA,
    "A",
    "Proteomic coverage of the frozen architecture",
)


# ============================================================
# PANEL B — gene-level transcript / protein concordance
# ============================================================

concordant = (
    B["gene_direction_concordant"]
    .astype(bool)
)

fdr_supported = (
    B["all_assays_FDR05"]
    .astype(bool)
)

xB = B["transcript_logFC"].to_numpy(float)

yB = (
    B["protein_median_standardized_effect"]
    .to_numpy(float)
)

axB.axhline(
    0,
    linewidth=1.0,
    linestyle="--",
    alpha=0.55,
)

axB.axvline(
    0,
    linewidth=1.0,
    linestyle="--",
    alpha=0.55,
)

# Discordant
axB.scatter(
    xB[~concordant],
    yB[~concordant],
    s=96,
    facecolor="white",
    edgecolor="0.45",
    linewidth=1.0,
    zorder=2,
)

# Concordant
axB.scatter(
    xB[concordant],
    yB[concordant],
    s=88,
    edgecolor="black",
    linewidth=0.45,
    zorder=3,
)

# FDR-supported assays receive a larger open ring.
axB.scatter(
    xB[fdr_supported],
    yB[fdr_supported],
    s=140,
    facecolor="none",
    edgecolor="black",
    linewidth=1.35,
    zorder=4,
)


# ------------------------------------------------------------
# Deterministic sparse labels:
# 6 strongest concordant-FDR genes + 2 strongest
# discordant-FDR genes by minimum assay FDR.
# ------------------------------------------------------------

conc_fdr = (
    B[
        concordant &
        fdr_supported
    ]
    .sort_values(
        ["min_assay_FDR", "gene_symbol"]
    )
    .head(6)
)

disc_fdr = (
    B[
        (~concordant) &
        fdr_supported
    ]
    .sort_values(
        ["min_assay_FDR", "gene_symbol"]
    )
    .head(2)
)

labelB = pd.concat(
    [conc_fdr, disc_fdr],
    ignore_index=True,
)

B_offsets = {
    "STAT1": (7, 8, "left", "bottom"),
    "ISG15": (7, -8, "left", "top"),
    "LAP3": (-8, 8, "right", "bottom"),
    "ZBP1": (-8, -8, "right", "top"),
    "MX1": (8, 8, "left", "bottom"),
    "FKBP5": (-8, 8, "right", "bottom"),
    "MATK": (-8, -10, "right", "top"),
}

for _, r in labelB.iterrows():

    # ATG7 is positioned separately below with a leader line.
    if r["gene_symbol"] == "ATG7":
        continue

    dx, dy, ha, va = B_offsets.get(
        r["gene_symbol"],
        (7, 7, "left", "bottom"),
    )

    axB.annotate(
        r["gene_symbol"],
        (
            r["transcript_logFC"],
            r["protein_median_standardized_effect"],
        ),
        xytext=(dx, dy),
        textcoords="offset points",
        fontsize=14,
        fontweight="bold",
        ha=ha,
        va=va,
    )


# ------------------------------------------------------------
# ATG7 — explicit placement with leader line
# ------------------------------------------------------------

atg7_row = B.loc[
    B["gene_symbol"] == "ATG7"
]

if len(atg7_row) != 1:
    raise RuntimeError(
        f"Expected exactly one ATG7 row; found {len(atg7_row)}."
    )

atg7_row = atg7_row.iloc[0]

axB.annotate(
    "ATG7",
    xy=(
        atg7_row["transcript_logFC"],
        atg7_row["protein_median_standardized_effect"],
    ),
    xytext=(42, 18),
    textcoords="offset points",
    fontsize=10.5,
    fontweight="bold",
    ha="left",
    va="bottom",
    arrowprops=dict(
        arrowstyle="-",
        linewidth=0.9,
        color="0.30",
        shrinkA=2,
        shrinkB=4,
    ),
    annotation_clip=False,
)


summary_all = (
    B.groupby(
        lambda _: True
    )
)

n_concordant = int(
    concordant.sum()
)

n_fdr_concordant = int(
    (concordant & fdr_supported).sum()
)

n_fdr_discordant = int(
    ((~concordant) & fdr_supported).sum()
)

axB.text(
    0.0,
    1.015,
    (
        f"{n_concordant}/84 direction concordant "
        f"({100*n_concordant/84:.1f}%)  ·  "
        f"FDR-supported: {n_fdr_concordant} concordant, "
        f"{n_fdr_discordant} discordant"
    ),
    transform=axB.transAxes,
    fontsize=15,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)

axB.set_xlabel(
    "Influenza transcript effect (logFC)"
)

axB.set_ylabel(
    "Protein standardized effect"
)

legend_B = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="",
        markersize=7,
        markerfacecolor="white",
        markeredgecolor="0.45",
        label="Direction discordant",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="",
        markersize=7,
        markerfacecolor="0.55",
        markeredgecolor="black",
        markeredgewidth=0.45,
        label="Direction concordant",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="",
        markersize=10,
        markerfacecolor="none",
        markeredgecolor="black",
        markeredgewidth=1.3,
        label="All assays FDR < 0.05",
    ),
]

axB.legend(
    handles=legend_B,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.25),
    ncol=3,
    frameon=False,
    borderaxespad=0.0,
    handletextpad=0.45,
    columnspacing=1.2,
)

add_panel_header(
    axB,
    "B",
    "Gene-level transcript–protein concordance",
)


# ============================================================
# PANEL C — seven biological domains
# ============================================================

domain_display = {
    "Interferon": "Interferon",
    "Innate_Cytokine": "Innate cytokine",
    "Antigen_Presentation": "Antigen presentation",
    "Neutrophil_Monocyte": "Neutrophil / monocyte",
    "Adaptive_Immunity": "Adaptive immunity",
    "Metabolism": "Metabolism",
    "Translation": "Translation",
}

yC = np.arange(len(C))

trans = C["median_transcript_NES"].to_numpy(float)
prot = C["median_protein_NES"].to_numpy(float)

for i in range(len(C)):
    axC.plot(
        [trans[i], prot[i]],
        [i, i],
        linewidth=1.5,
        alpha=0.55,
        zorder=1,
    )

axC.scatter(
    trans,
    yC,
    s=78,
    marker="o",
    edgecolor="black",
    linewidth=0.5,
    label="Transcript median NES",
    zorder=3,
)

axC.scatter(
    prot,
    yC,
    s=100,
    marker="s",
    edgecolor="black",
    linewidth=0.5,
    label="Protein median NES",
    zorder=3,
)

axC.axvline(
    0,
    linestyle="--",
    linewidth=1.0,
    alpha=0.55,
)

axC.set_yticks(yC)

axC.set_yticklabels(
    [
        domain_display[x]
        for x in C["preregistered_domain"]
    ],
    fontweight="bold",
)

axC.invert_yaxis()

axC.set_xlabel(
    "Median normalized enrichment score"
)

xmin = min(
    float(np.nanmin(trans)),
    float(np.nanmin(prot)),
)

xmax = max(
    float(np.nanmax(trans)),
    float(np.nanmax(prot)),
)

span = xmax - xmin

axC.set_xlim(
    xmin - 0.08 * span,
    xmax + 0.38 * span,
)

# Compact pathway-count annotation strip INSIDE Panel C.
# Using axes coordinates guarantees that annotations cannot
# intrude into Panel D.

axC.text(
    0.985,
    1.015,
    "Pathways",
    transform=axC.transAxes,
    fontsize=15,
    fontweight="bold",
    ha="right",
    va="bottom",
    clip_on=False,
)

for i, (_, r) in enumerate(C.iterrows()):

    y_frac = 1.0 - (i + 0.5) / len(C)

    axC.text(
        0.985,
        y_frac,
        (
            f"{int(r['n_direction_concordant'])}/"
            f"{int(r['n_unique_pathways'])}; "
            f"{int(r['n_concordant_FDR05'])} FDR"
        ),
        transform=axC.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="right",
        va="center",
        clip_on=False,
    )

axC.legend(
    loc="lower left",
    bbox_to_anchor=(0.0, -0.24),
    ncol=2,
    frameon=False,
    borderaxespad=0.0,
    handletextpad=0.45,
    columnspacing=1.2,
)

add_panel_header(
    axC,
    "C",
    "Proteomic support across biological domains",
)


# ============================================================
# PANEL D — cross-layer validation summary
# ============================================================

yD = np.arange(len(D))

fraction = (
    D["concordance_fraction"]
    .to_numpy(float)
)

axD.barh(
    yD,
    fraction,
    height=0.56,
)

axD.axvline(
    0.5,
    linestyle="--",
    linewidth=1.0,
    alpha=0.45,
)

axD.set_yticks(yD)

axD.set_yticklabels(
    D["validation_layer"],
    fontweight="bold",
)

axD.invert_yaxis()

axD.set_xlim(
    0,
    1.08,
)

axD.set_xlabel(
    "Direction concordance fraction"
)

axD.set_xticks(
    np.arange(0, 1.01, 0.2)
)

for i, (_, r) in enumerate(D.iterrows()):

    axD.text(
        r["concordance_fraction"] + 0.018,
        i - 0.06,
        (
            f"{int(r['n_direction_concordant'])}/"
            f"{int(r['n_testable'])} "
            f"({100*r['concordance_fraction']:.1f}%)"
        ),
        fontsize=15,
        fontweight="bold",
        va="center",
        ha="left",
    )

    axD.text(
        0.02,
        i + 0.18,
        (
            f"FDR: {int(r['FDR_concordant'])} concordant / "
            f"{int(r['FDR_discordant'])} discordant"
        ),
        fontsize=14,
        fontweight="bold",
        va="center",
        ha="left",
    )

    
add_panel_header(
    axD,
    "D",
    "Concordance across validation layers",
)

axD.text(
    0.0,
    1.015,
    (
        "Gene, pathway, and regulatory-program results are "
        "distinct prespecified validation layers."
    ),
    transform=axD.transAxes,
    fontsize=15,
    fontweight="bold",
    ha="left",
    va="bottom",
)


# ============================================================
# Save
# ============================================================

fig.savefig(
    PDF,
    bbox_inches="tight",
)

fig.savefig(
    PNG,
    dpi=600,
    bbox_inches="tight",
)

fig.savefig(
    SVG,
    bbox_inches="tight",
)

plt.close(fig)


# ============================================================
# Console QC
# ============================================================

print("=== FIGURE S8 V7 GENERATION COMPLETE ===")
print(PDF)
print(PNG)
print(SVG)

print()
print("=== PANEL A ===")
print(
    A[
        [
            "architecture_class",
            "n_frozen_genes",
            "n_measurable_genes",
            "coverage_percent",
        ]
    ].to_string(index=False)
)

print()
print(
    f"Overall final frozen mapping: "
    f"{overall_measured}/{overall_total}"
)

print()
print("=== PANEL B ===")
print("Genes:", len(B))
print("Direction concordant:", n_concordant)
print("FDR concordant:", n_fdr_concordant)
print("FDR discordant:", n_fdr_discordant)

print()
print("=== PANEL C ===")
print(
    C[
        [
            "preregistered_domain",
            "n_unique_pathways",
            "n_direction_concordant",
            "n_concordant_FDR05",
            "median_transcript_NES",
            "median_protein_NES",
        ]
    ].to_string(index=False)
)

print()
print("=== PANEL D ===")
print(
    D[
        [
            "validation_layer",
            "n_testable",
            "n_direction_concordant",
            "concordance_fraction",
            "FDR_concordant",
            "FDR_discordant",
        ]
    ].to_string(index=False)
)

print()
print(
    "No differential proteomics was rerun.\n"
    "No SomaScan mapping was rerun.\n"
    "No FGSEA pathway analysis was rerun.\n"
    "No regulatory-program analysis was rerun.\n"
    "No validation category was recalculated.\n"
    "Figure S8 v7 visualizes frozen proteomics outputs only."
)
