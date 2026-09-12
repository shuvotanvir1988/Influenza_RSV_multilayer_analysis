#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import (
    LinearSegmentedColormap,
    ListedColormap,
    BoundaryNorm,
)
from matplotlib.lines import Line2D


# ============================================================
# Paths
# ============================================================

ROOT = Path.home() / "Influenza_RSV_Project"

REG = ROOT / "results/regulatory_driver_analysis"

MASTER37 = (
    REG /
    "tables/REGULATORY_HIGH_CONFIDENCE_DRIVERS_v1.0.tsv"
)

RNASEQ37 = (
    REG /
    "rnaseq_validation/"
    "GSE155925_frozen37_driver_validation_v1.0.tsv"
)

CELL_SUPPORT = (
    REG /
    "tables/"
    "REGULATORY_DRIVER_37x11_CELLTYPE_SIGNED_SUPPORT_MATRIX_v1.0.tsv"
)

BASE = (
    ROOT /
    "results/supplementary_figures/"
    "FigureS6_regulatory_architecture"
)

FIGDIR = BASE / "figures"
TABDIR = BASE / "tables"

FIGDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Output filenames — v11
# ============================================================

PDF = (
    FIGDIR /
    "FigureS6_v11_extended_regulatory_architecture.pdf"
)

PNG = (
    FIGDIR /
    "FigureS6_v11_extended_regulatory_architecture.png"
)

SVG = (
    FIGDIR /
    "FigureS6_v11_extended_regulatory_architecture.svg"
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
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "legend.fontsize": 16,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ============================================================
# Load frozen inputs
# ============================================================

drivers = pd.read_csv(
    MASTER37,
    sep="\t",
)

rnaseq = pd.read_csv(
    RNASEQ37,
    sep="\t",
)

cell = pd.read_csv(
    CELL_SUPPORT,
    sep="\t",
)


# ============================================================
# Hard QC
# ============================================================

if len(drivers) != 37:
    raise RuntimeError(
        f"Expected 37 high-confidence drivers; found {len(drivers)}."
    )

if len(rnaseq) != 37:
    raise RuntimeError(
        f"Expected 37 GSE155925 validation rows; found {len(rnaseq)}."
    )

if len(cell) != 37:
    raise RuntimeError(
        f"Expected 37 rows in cell-support matrix; found {len(cell)}."
    )

if cell.shape[1] != 12:
    raise RuntimeError(
        "Expected TF + 11 cell-type columns in signed-support matrix; "
        f"found shape {cell.shape}."
    )

if drivers["TF"].duplicated().any():
    raise RuntimeError("Duplicate TFs in frozen 37-driver master table.")

if rnaseq["TF"].duplicated().any():
    raise RuntimeError("Duplicate TFs in GSE155925 validation table.")

if cell["TF"].duplicated().any():
    raise RuntimeError("Duplicate TFs in cell-support matrix.")

if set(drivers["TF"]) != set(rnaseq["TF"]):
    raise RuntimeError(
        "TF identities differ between frozen37 master and RNA-seq validation."
    )

if set(drivers["TF"]) != set(cell["TF"]):
    raise RuntimeError(
        "TF identities differ between frozen37 master and cell-support matrix."
    )


# ============================================================
# Stable driver ordering
#
# Prioritize the frozen regulatory-driver score, then
# influenza activity magnitude. This is visualization ordering
# only and does not create a new ranking statistic.
# ============================================================

order = (
    drivers
    .sort_values(
        ["regulatory_driver_score", "flu_delta", "TF"],
        ascending=[False, False, True],
    )
    ["TF"]
    .tolist()
)

drivers = (
    drivers
    .set_index("TF")
    .loc[order]
    .reset_index()
)

rnaseq = (
    rnaseq
    .set_index("TF")
    .loc[order]
    .reset_index()
)

cell = (
    cell
    .set_index("TF")
    .loc[order]
    .reset_index()
)


# ============================================================
# Panel A source
# Complete frozen regulatory activity architecture
# ============================================================

panelA = drivers[
    [
        "TF",
        "flu_delta",
        "rsv_delta",
        "flu_vs_rsv_delta",
        "high_confidence_pathogen_biased_driver",
        "regulatory_driver_score",
    ]
].copy()

panelA.to_csv(
    TABDIR /
    "FigureS6_panelA_frozen37_regulatory_activity.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# Panel B source
# Connectivity to frozen gene architecture
# ============================================================

panelB = drivers[
    [
        "TF",
        "shared_core_coverage_fraction",
        "shared_influenza_amplified_coverage_fraction",
        "frozen170_coverage_fraction",
        "frozen170_targets",
        "frozen170_connectivity_tier",
    ]
].copy()

panelB.to_csv(
    TABDIR /
    "FigureS6_panelB_frozen_architecture_connectivity.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# Panel C source
# Independent bulk RNA-seq regulatory portability
# ============================================================

panelC = rnaseq[
    [
        "TF",
        "flu_delta",
        "rnaseq_delta_activity",
        "rnaseq_P",
        "rnaseq_FDR",
        "evaluable_rnaseq",
        "direction_concordant",
        "nominal_support",
        "FDR_support",
        "high_confidence_pathogen_biased_driver",
    ]
].copy()

panelC.to_csv(
    TABDIR /
    "FigureS6_panelC_GSE155925_frozen37_validation.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# Panel D source
# Frozen cell-type signed-support matrix
# ============================================================

panelD = cell.copy()

panelD.to_csv(
    TABDIR /
    "FigureS6_panelD_37x11_signed_celltype_support.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# Figure layout
# ============================================================

fig = plt.figure(
    figsize=(21.5, 26.0)
)

gs = fig.add_gridspec(
    nrows=3,
    ncols=2,
    height_ratios=[1.84, 1.02, 2.05],
    left=0.125,
    right=0.975,
    bottom=0.070,
    top=0.965,
    wspace=0.40,
    hspace=0.64,
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])

# Panel C uses a dedicated plot + information region.
gsC = gs[1, :].subgridspec(
    1,
    2,
    width_ratios=[3.45, 1.45],
    wspace=0.10,
)

axC = fig.add_subplot(gsC[0, 0])
axCinfo = fig.add_subplot(gsC[0, 1])

axD = fig.add_subplot(gs[2, :])


# ============================================================
# Utility functions
# ============================================================

def add_panel_label(ax, label):
    ax.text(
        -0.125,
        1.075,
        label,
        transform=ax.transAxes,
        fontsize=26,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def add_panel_title(ax, text):
    ax.text(
        0.0,
        1.075,
        text,
        transform=ax.transAxes,
        fontsize=23,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )


# ============================================================
# Panel A
# ============================================================

A_cols = [
    "flu_delta",
    "rsv_delta",
    "flu_vs_rsv_delta",
]

A_labels = [
    "Influenza vs\ncontrol",
    "RSV vs\ncontrol",
    "Influenza vs\nRSV",
]

A_matrix = panelA[A_cols].to_numpy(float)

absmax_A = np.nanmax(np.abs(A_matrix))

activity_cmap = LinearSegmentedColormap.from_list(
    "activity_diverging",
    ["#3B6FB6", "#F7F7F7", "#C9493D"],
)

imA = axA.imshow(
    A_matrix,
    aspect="auto",
    cmap=activity_cmap,
    vmin=-absmax_A,
    vmax=absmax_A,
    interpolation="nearest",
)

axA.set_xticks(
    np.arange(len(A_labels))
)

axA.set_xticklabels(
    A_labels,
    fontweight="bold",
)

axA.set_yticks(
    np.arange(len(panelA))
)

A_tf_labels = [
    ("● " + tf) if bool(flag) else ("   " + tf)
    for tf, flag in zip(
        panelA["TF"],
        panelA["high_confidence_pathogen_biased_driver"],
    )
]

axA.set_yticklabels(
    A_tf_labels,
    fontsize=14.0,
    fontweight="bold",
)

axA.tick_params(
    axis="both",
    length=0,
)


cbarA = fig.colorbar(
    imA,
    ax=axA,
    fraction=0.040,
    pad=0.025,
)

cbarA.set_label(
    "Differential TF activity",
    fontsize=17.0,
    fontweight="bold",
)


add_panel_title(
    axA,
    "Extended regulatory-driver activity"
)

add_panel_label(
    axA,
    "A"
)

axA.text(
    0.0,
    -0.125,
    "● Focused pathogen-biased driver",
    transform=axA.transAxes,
    fontsize=14.5,
    fontweight="bold",
    ha="left",
    va="top",
    clip_on=False,
)



# ============================================================
# Panel B
# ============================================================

B_cols = [
    "shared_core_coverage_fraction",
    "shared_influenza_amplified_coverage_fraction",
    "frozen170_coverage_fraction",
]

B_labels = [
    "Shared core",
    "Influenza-\nenhanced\ncomponent",
    "Frozen\n170-gene\nresponse",
]

B_matrix = panelB[B_cols].to_numpy(float)

imB = axB.imshow(
    B_matrix,
    aspect="auto",
    cmap="Blues",
    vmin=0,
    vmax=max(
        0.20,
        np.nanmax(B_matrix),
    ),
    interpolation="nearest",
)

axB.set_xticks(
    np.arange(len(B_labels))
)

axB.set_xticklabels(
    B_labels,
    fontsize=15.5,
    fontweight="bold",
    linespacing=1.05,
)

axB.tick_params(
    axis="x",
    pad=10,
)

axB.set_yticks(
    np.arange(len(panelB))
)

axB.set_yticklabels(
    [""] * len(panelB)
)

axB.tick_params(
    axis="both",
    length=0,
)

cbarB = fig.colorbar(
    imB,
    ax=axB,
    fraction=0.040,
    pad=0.025,
)

cbarB.set_label(
    "Target coverage fraction",
    fontsize=17.0,
    fontweight="bold",
)

add_panel_title(
    axB,
    "Connectivity to frozen host-response architecture"
)

add_panel_label(
    axB,
    "B"
)

axB.text(
    0.0,
    -0.165,
    "TF order as in Panel A",
    transform=axB.transAxes,
    fontsize=14.5,
    fontweight="bold",
    ha="left",
    va="top",
)


# ============================================================
# Panel C
# ============================================================

C = panelC.loc[
    panelC["evaluable_rnaseq"].astype(bool)
].copy()

if len(C) != 37:
    raise RuntimeError(
        f"Expected all 37 drivers evaluable in RNA-seq; found {len(C)}."
    )

concordant = C[
    "direction_concordant"
].astype(bool)

xC = C["flu_delta"].to_numpy(float)
yC = C["rnaseq_delta_activity"].to_numpy(float)

# Identity/reference axes.
axC.axhline(
    0,
    color="0.55",
    linewidth=0.9,
    zorder=0,
)

axC.axvline(
    0,
    color="0.55",
    linewidth=0.9,
    zorder=0,
)

# Directionally discordant
axC.scatter(
    xC[~concordant],
    yC[~concordant],
    s=82,
    marker="o",
    facecolor="white",
    edgecolor="0.45",
    linewidth=1.0,
    label="Direction discordant",
    zorder=2,
)

# Directionally concordant
axC.scatter(
    xC[concordant],
    yC[concordant],
    s=90,
    marker="o",
    facecolor="#3B6FB6",
    edgecolor="black",
    linewidth=0.45,
    label="Direction concordant",
    zorder=3,
)

# FDR-supported points, if any.
fdr_supported = C[
    "FDR_support"
].astype(bool)

if fdr_supported.any():
    axC.scatter(
        C.loc[fdr_supported, "flu_delta"],
        C.loc[fdr_supported, "rnaseq_delta_activity"],
        s=135,
        facecolor="none",
        edgecolor="#C9493D",
        linewidth=1.7,
        label="RNA-seq FDR support",
        zorder=4,
    )

# Deliberately omit individual TF labels.
# Panel C communicates overall portability of the frozen
# regulatory architecture rather than emphasizing selected TFs.

axC.set_xlabel(
    "Discovery influenza differential TF activity"
)

axC.set_ylabel(
    "Independent RNA-seq differential TF activity"
)

# Deliberately compact plotting limits around the frozen
# validation observations.
x_pad = 0.08
y_pad = 0.04

axC.set_xlim(
    min(-0.05, float(C["flu_delta"].min()) - x_pad),
    float(C["flu_delta"].max()) + x_pad,
)

axC.set_ylim(
    float(C["rnaseq_delta_activity"].min()) - y_pad,
    float(C["rnaseq_delta_activity"].max()) + y_pad,
)

add_panel_title(
    axC,
    "Independent bulk RNA-seq portability of frozen regulators"
)

n_conc = int(
    C["direction_concordant"].astype(bool).sum()
)

n_nom = int(
    C["nominal_support"].astype(bool).sum()
)

n_fdr = int(
    C["FDR_support"].astype(bool).sum()
)

summary_text = (
    f"Evaluable: {len(C)}/37   |   "
    f"Direction concordant: {n_conc}/37   |   "
    f"Nominal support: {n_nom}/37   |   "
    f"FDR support: {n_fdr}/37"
)


add_panel_label(
    axC,
    "C"
)

# ----------------------------------------------------------
# Dedicated Panel C information/key area.
# ----------------------------------------------------------

axCinfo.axis("off")

axCinfo.text(
    0.02,
    0.94,
    "Frozen-regulator portability",
    transform=axCinfo.transAxes,
    fontsize=17.0,
    fontweight="bold",
    ha="left",
    va="top",
)

axCinfo.text(
    0.02,
    0.79,
    "37 / 37",
    transform=axCinfo.transAxes,
    fontsize=20,
    fontweight="bold",
    ha="left",
    va="center",
)

axCinfo.text(
    0.33,
    0.79,
    "evaluable",
    transform=axCinfo.transAxes,
    fontsize=14.0,
    ha="left",
    va="center",
)

axCinfo.text(
    0.02,
    0.63,
    f"{n_conc} / 37",
    transform=axCinfo.transAxes,
    fontsize=20,
    fontweight="bold",
    ha="left",
    va="center",
)

axCinfo.text(
    0.33,
    0.63,
    "direction concordant",
    transform=axCinfo.transAxes,
    fontsize=14.0,
    ha="left",
    va="center",
)

axCinfo.text(
    0.02,
    0.47,
    f"{n_nom} / 37",
    transform=axCinfo.transAxes,
    fontsize=20,
    fontweight="bold",
    ha="left",
    va="center",
)

axCinfo.text(
    0.33,
    0.47,
    "nominal support",
    transform=axCinfo.transAxes,
    fontsize=14.0,
    ha="left",
    va="center",
)

axCinfo.text(
    0.02,
    0.31,
    f"{n_fdr} / 37",
    transform=axCinfo.transAxes,
    fontsize=20,
    fontweight="bold",
    ha="left",
    va="center",
)

axCinfo.text(
    0.33,
    0.31,
    "FDR support",
    transform=axCinfo.transAxes,
    fontsize=14.0,
    ha="left",
    va="center",
)

panelC_handles = [
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
        markerfacecolor="#3B6FB6",
        markeredgecolor="black",
        markeredgewidth=0.45,
        label="Direction concordant",
    ),
]

axCinfo.legend(
    handles=panelC_handles,
    loc="lower left",
    bbox_to_anchor=(0.0, 0.03),
    frameon=False,
    fontsize=14.0,
    handletextpad=0.5,
    borderaxespad=0.0,
    labelspacing=0.8,
)







# ============================================================
# Panel D
# ============================================================

cell_types = [
    c for c in panelD.columns
    if c != "TF"
]

D_matrix = (
    panelD[cell_types]
    .to_numpy(float)
)

unique_D = set(
    np.unique(
        D_matrix[
            ~np.isnan(D_matrix)
        ]
    )
)

if not unique_D.issubset({-1.0, 0.0, 1.0}):
    raise RuntimeError(
        "Unexpected value in signed-support matrix: "
        f"{sorted(unique_D)}"
    )

support_cmap = ListedColormap(
    [
        "#3B6FB6",
        "#F2F2F2",
        "#C9493D",
    ]
)

support_norm = BoundaryNorm(
    [-1.5, -0.5, 0.5, 1.5],
    support_cmap.N,
)

imD = axD.imshow(
    D_matrix,
    aspect="auto",
    cmap=support_cmap,
    norm=support_norm,
    interpolation="nearest",
)

axD.set_xticks(
    np.arange(len(cell_types))
)

axD.set_xticklabels(
    cell_types,
    rotation=35,
    ha="right",
    rotation_mode="anchor",
    fontsize=20.0,
    fontweight="bold",
)

axD.tick_params(
    axis="x",
    pad=10,
)

axD.set_yticks(
    np.arange(len(panelD))
)

axD.set_yticklabels(
    panelD["TF"],
    fontsize=15.5,
    fontweight="bold",
)

axD.tick_params(
    axis="both",
    length=0,
)

# Thin grid to separate matrix cells.
axD.set_xticks(
    np.arange(-0.5, len(cell_types), 1),
    minor=True,
)

axD.set_yticks(
    np.arange(-0.5, len(panelD), 1),
    minor=True,
)

axD.grid(
    which="minor",
    linewidth=0.35,
    color="white",
)

axD.tick_params(
    which="minor",
    bottom=False,
    left=False,
)

add_panel_title(
    axD,
    "Cell-type support for the RSV-derived regulatory direction"
)

add_panel_label(
    axD,
    "D"
)

legend_D = [
    Line2D(
        [0],
        [0],
        marker="s",
        linestyle="",
        markersize=9,
        markerfacecolor="#C9493D",
        markeredgecolor="none",
        label="Supported, concordant direction",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        linestyle="",
        markersize=9,
        markerfacecolor="#F2F2F2",
        markeredgecolor="0.7",
        label="No FDR-supported activity",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        linestyle="",
        markersize=9,
        markerfacecolor="#3B6FB6",
        markeredgecolor="none",
        label="Supported, opposite direction",
    ),
]

axD.legend(
    handles=legend_D,
    loc="lower center",
    bbox_to_anchor=(0.50, 1.015),
    ncol=3,
    frameon=False,
    fontsize=15.0,
    columnspacing=1.35,
    handletextpad=0.4,
    borderaxespad=0.0,
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

print("=== FIGURE S6 V11 GENERATION COMPLETE ===")
print(PDF)
print(PNG)
print(SVG)

print()
print("=== PANEL A ===")
print("Drivers:", len(panelA))
print(
    "Focused pathogen-biased:",
    int(
        panelA[
            "high_confidence_pathogen_biased_driver"
        ].astype(bool).sum()
    )
)

print()
print("=== PANEL B ===")
print(
    "Frozen170 target range:",
    int(panelB["frozen170_targets"].min()),
    "to",
    int(panelB["frozen170_targets"].max()),
)

print()
print("=== PANEL C ===")
print("Evaluable:", len(C))
print("Direction concordant:", n_conc)
print("Nominal support:", n_nom)
print("FDR support:", n_fdr)

print()
print("=== PANEL D ===")
print("Drivers:", len(panelD))
print("Cell types:", len(cell_types))
print("Matrix values:", sorted(unique_D))

print()
print(
    "No regulatory inference, ULM analysis, limma model, "
    "RNA-seq validation model, or single-cell analysis was rerun."
)

print(
    "Figure S6 v11 is visualization-only and consumes frozen outputs."
)
