#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path.home() / "Influenza_RSV_Project"

SIG = (
    ROOT / "results/DS001_GSE38900/"
    "signature_analysis/tables"
)

FROZEN170 = (
    ROOT / "results/figure2_influenza_centered/"
    "frozen_submission_v5/tables/"
    "Figure2C_v8_replicated_170_gene_architecture.tsv"
)

BASE = (
    ROOT / "results/supplementary_figures/"
    "FigureS3_signature_robustness"
)

OUT = BASE / "figures"
TAB = BASE / "tables"

OUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

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
    "legend.fontsize": 16,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

COLORS = {
    "Shared core": "#7A7A7A",
    "Influenza-amplified": "#D95F02",
    "RSV-amplified": "#1B9E77",
}

# ---------------------------------------------------------
# Load frozen sources
# ---------------------------------------------------------

x170 = pd.read_csv(FROZEN170, sep="\t")

required_170 = {
    "gene_symbol",
    "flu_logFC",
    "rsv6884_logFC",
    "architecture_class",
}

missing = required_170 - set(x170.columns)

if missing:
    raise RuntimeError(
        f"Frozen 170-gene table missing columns: {missing}"
    )

if len(x170) != 170:
    raise RuntimeError(
        f"Expected frozen 170-gene set, found {len(x170)} rows."
    )

expected_counts = {
    "Shared core": 130,
    "Influenza-amplified": 37,
    "RSV-amplified": 3,
}

observed_counts = (
    x170["architecture_class"]
    .value_counts()
    .to_dict()
)

if observed_counts != expected_counts:
    raise RuntimeError(
        f"Unexpected frozen class counts: {observed_counts}"
    )

robust = pd.read_csv(
    SIG / "DS001_GSE38900_shared_signature_robustness_summary.tsv",
    sep="\t",
)

robust_map = dict(
    zip(
        robust["metric"],
        robust["value"],
    )
)

needed_robust = {
    "shared_GPL6884_total",
    "shared_GPL6884_crossplatform_evaluable",
    "shared_GPL6884_RSV_replicated",
    "shared_replicated_direct_not_significant",
    "shared_replicated_direct_significant",
    "shared_replicated_influenza_biased",
    "shared_replicated_RSV_biased",
}

missing_robust = needed_robust - set(robust_map)

if missing_robust:
    raise RuntimeError(
        f"Missing robustness metrics: {missing_robust}"
    )


# Save exact figure source tables.
x170.to_csv(
    TAB / "FigureS3_frozen170_source.tsv",
    sep="\t",
    index=False,
)

robust.to_csv(
    TAB / "FigureS3_derivation_source.tsv",
    sep="\t",
    index=False,
)


# ---------------------------------------------------------
# Figure canvas
# ---------------------------------------------------------

fig = plt.figure(figsize=(20.8, 15.0))

gs = fig.add_gridspec(
    2, 2,
    left=0.070,
    right=0.975,
    bottom=0.125,
    top=0.945,
    wspace=0.44,
    hspace=0.50,
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

def panel_label(ax, lab):
    ax.text(
        -0.15, 1.08, lab,
        transform=ax.transAxes,
        fontsize=25,
        fontweight="bold",
        ha="left",
        va="top",
    )

# ---------------------------------------------------------
# Panel A — filtering cascade
# ---------------------------------------------------------

axA.axis("off")

axA.set_title(
    "Derivation of the replicated shared architecture",
    loc="left",
    pad=16,
)

def flow_box(
    ax,
    x,
    y,
    w,
    h,
    number,
    label,
    fontsize_number=24,
    fontsize_label=16,
):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        transform=ax.transAxes,
        boxstyle="round,pad=0.018,rounding_size=0.018",
        linewidth=1.8,
        facecolor="white",
        edgecolor="black",
    )

    ax.add_patch(box)

    ax.text(
        x + w/2,
        y + h*0.64,
        f"{number:,}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fontsize_number,
        fontweight="bold",
    )

    ax.text(
        x + w/2,
        y + h*0.27,
        label,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fontsize_label,
        fontweight="bold",
    )


# Three-step vertical filtering cascade.
x_main = 0.14
w_main = 0.72
h_main = 0.145

y1 = 0.77
y2 = 0.56
y3 = 0.35

flow_box(
    axA,
    x_main,
    y1,
    w_main,
    h_main,
    1735,
    "Shared on GPL6884",
)

flow_box(
    axA,
    x_main,
    y2,
    w_main,
    h_main,
    1544,
    "Cross-platform evaluable",
)

flow_box(
    axA,
    x_main,
    y3,
    w_main,
    h_main,
    170,
    "Cross-platform replicated shared genes",
)

# Vertical arrows between stages.
for y_top, y_bottom in [
    (y1, y2),
    (y2, y3),
]:
    axA.annotate(
        "",
        xy=(0.50, y_bottom + h_main + 0.012),
        xytext=(0.50, y_top - 0.012),
        xycoords=axA.transAxes,
        textcoords=axA.transAxes,
        arrowprops=dict(
            arrowstyle="->",
            lw=2.0,
        ),
    )


# Final split of the frozen 170 genes.
left_x = 0.06
right_x = 0.55
bottom_y = 0.055
bottom_w = 0.39
bottom_h = 0.17

flow_box(
    axA,
    left_x,
    bottom_y,
    bottom_w,
    bottom_h,
    130,
    "Shared core",
    fontsize_number=22,
    fontsize_label=16,
)

flow_box(
    axA,
    right_x,
    bottom_y,
    bottom_w,
    bottom_h,
    40,
    "Directly divergent",
    fontsize_number=22,
    fontsize_label=16,
)

# Branch arrows from the 170-gene box.
axA.annotate(
    "",
    xy=(left_x + bottom_w/2, bottom_y + bottom_h + 0.01),
    xytext=(0.47, y3 - 0.012),
    xycoords=axA.transAxes,
    textcoords=axA.transAxes,
    arrowprops=dict(
        arrowstyle="->",
        lw=1.9,
    ),
)

axA.annotate(
    "",
    xy=(right_x + bottom_w/2, bottom_y + bottom_h + 0.01),
    xytext=(0.53, y3 - 0.012),
    xycoords=axA.transAxes,
    textcoords=axA.transAxes,
    arrowprops=dict(
        arrowstyle="->",
        lw=1.9,
    ),
)

# Explain the 40-gene divergent subset without another box.

# ---------------------------------------------------------
# Panel B — class composition
# ---------------------------------------------------------

axB.set_title(
    "Composition of the frozen 170-gene set",
    loc="left",
    pad=16,
)

classes = [
    "Shared core",
    "Influenza-amplified",
    "RSV-amplified",
]

counts = [
    expected_counts[c]
    for c in classes
]

bars = axB.barh(
    np.arange(len(classes)),
    counts,
    color=[COLORS[c] for c in classes],
)

axB.set_yticks(np.arange(len(classes)))
axB.set_yticklabels(
    classes,
    fontweight="bold",
)

axB.invert_yaxis()

axB.set_xlabel("Genes", labelpad=10)

for bar, n in zip(bars, counts):
    axB.text(
        bar.get_width() + 2,
        bar.get_y() + bar.get_height()/2,
        f"{n}",
        va="center",
        fontsize=17,
        fontweight="bold",
    )

axB.set_xlim(0, 145)
axB.grid(axis="x", alpha=0.15)

# ---------------------------------------------------------
# Panel C — frozen 170 effect-size structure
# ---------------------------------------------------------

axC.set_title(
    "Effect-size structure of the frozen 170 genes",
    loc="left",
    pad=16,
)

for cls in classes:
    d = x170[x170["architecture_class"] == cls]

    size = 85 if cls != "RSV-amplified" else 115

    axC.scatter(
        d["rsv6884_logFC"],
        d["flu_logFC"],
        s=size,
        c=COLORS[cls],
        alpha=0.82,
        edgecolors="white",
        linewidths=0.8,
        label=f"{cls} (n={len(d)})",
    )

vals = np.concatenate([
    x170["rsv6884_logFC"].to_numpy(float),
    x170["flu_logFC"].to_numpy(float),
])

lo = vals.min()
hi = vals.max()
pad = (hi - lo) * 0.07

axC.plot(
    [lo-pad, hi+pad],
    [lo-pad, hi+pad],
    linestyle="--",
    color="black",
    linewidth=1.5,
    alpha=0.6,
)

axC.axhline(
    0,
    color="black",
    linewidth=1.35,
    alpha=0.15,
)

axC.axvline(
    0,
    color="black",
    linewidth=1.35,
    alpha=0.15,
)

axC.set_xlim(lo-pad, hi+pad)
axC.set_ylim(lo-pad, hi+pad)

axC.set_xlabel("RSV vs control logFC", labelpad=10)
axC.set_ylabel("Influenza vs control logFC", labelpad=10)

axC.legend(
    frameon=False,
    loc="upper left",
)

# Label the three RSV-amplified genes because n=3.
# Use explicit separated offsets to prevent label collisions.
rsv3 = x170[
    x170["architecture_class"] == "RSV-amplified"
]

panelC_offsets = {
    "CCR2": (12, 12),
    "LOC650518": (12, -15),
    "TIMD4": (12, 2),
}

for _, row in rsv3.iterrows():

    gene = row["gene_symbol"]

    dx, dy = panelC_offsets.get(
        gene,
        (12, 0),
    )

    axC.annotate(
        gene,
        (
            row["rsv6884_logFC"],
            row["flu_logFC"],
        ),
        xytext=(dx, dy),
        textcoords="offset points",
        ha="left",
        va="center",
        fontsize=14.5,
        fontweight="bold",
        arrowprops=dict(
            arrowstyle="-",
            lw=1.0,
            color="black",
            shrinkA=2,
            shrinkB=3,
        ),
    )

# ---------------------------------------------------------
# Panel D — influenza amplification within frozen 37 genes
# ---------------------------------------------------------

axD.set_title(
    "Magnitude of influenza amplification",
    loc="left",
    pad=16,
)

flu37 = (
    x170[
        x170["architecture_class"] == "Influenza-amplified"
    ]
    .copy()
)

if len(flu37) != 37:
    raise RuntimeError(
        f"Expected 37 influenza-amplified genes, found {len(flu37)}"
    )

# This is descriptive arithmetic from frozen effects only.
# No statistical model or classification is recomputed.
flu37["amplification_logFC"] = (
    flu37["flu_logFC"] -
    flu37["rsv6884_logFC"]
)

flu37 = flu37.sort_values(
    "amplification_logFC",
    ascending=True,
).reset_index(drop=True)

flu37.to_csv(
    TAB / "FigureS3_panelD_frozen37_influenza_amplification.tsv",
    sep="\t",
    index=False,
)

y = np.arange(len(flu37))

axD.hlines(
    y,
    0,
    flu37["amplification_logFC"],
    linewidth=1.35,
    alpha=0.45,
    color=COLORS["Influenza-amplified"],
)

axD.scatter(
    flu37["amplification_logFC"],
    y,
    s=58,
    color=COLORS["Influenza-amplified"],
    edgecolors="white",
    linewidths=0.7,
    zorder=3,
)

axD.axvline(
    0,
    color="black",
    linewidth=1.35,
    alpha=0.45,
)

axD.set_xlabel(
    "Comparative influenza amplification (ΔlogFC)",
    labelpad=10,
)

axD.set_ylabel("")

# Label only the strongest amplification signals to keep
# the panel readable while plotting all 37 frozen genes.
top_n = 12

top_idx = set(
    flu37.nlargest(
        top_n,
        "amplification_logFC"
    ).index
)

# Suppress the 37-row y-axis label/tick structure.
axD.set_yticks([])

# Put the strongest 12 signals into a dedicated label column.
# This preserves every point while preventing text collisions.
top12 = (
    flu37.loc[list(top_idx)]
    .sort_values(
        "amplification_logFC",
        ascending=False,
    )
    .copy()
)

label_x = (
    flu37["amplification_logFC"].max()
    + 0.18
)

# Evenly spaced fixed label positions spanning a much taller
# vertical range than v9. This prevents the top-12 gene labels
# from stacking on one another while preserving the exact ranking.
label_y_positions = np.linspace(
    len(flu37) - 0.5,
    len(flu37) - 17.0,
    len(top12),
)

for label_y, (_, row) in zip(
    label_y_positions,
    top12.iterrows(),
):
    data_y = row.name

    axD.annotate(
        row["gene_symbol"],
        xy=(
            row["amplification_logFC"],
            data_y,
        ),
        xytext=(
            label_x,
            label_y,
        ),
        textcoords="data",
        ha="left",
        va="center",
        fontsize=13.6,
        fontweight="bold",
        arrowprops=dict(
            arrowstyle="-",
            lw=0.85,
            color="black",
            alpha=0.55,
            shrinkA=2,
            shrinkB=3,
            connectionstyle="arc3,rad=0",
        ),
        clip_on=False,
    )

axD.grid(
    axis="x",
    alpha=0.15,
)

# Reserve explicit horizontal space for the label column.
dmax = flu37["amplification_logFC"].max()

axD.set_xlim(
    min(-0.05, flu37["amplification_logFC"].min() - 0.08),
    dmax + 1.15,
)

axD.text(
    0.02,
    -0.19,
    "All 37 frozen influenza-amplified genes shown; "
    "top 12 signals labeled",
    transform=axD.transAxes,
    ha="left",
    va="top",
    fontsize=14.0,
    clip_on=False,
)

axD.text(
    0.02,
    -0.255,
    "ΔlogFC = influenza-vs-control logFC − RSV-vs-control logFC",
    transform=axD.transAxes,
    ha="left",
    va="top",
    fontsize=13.5,
    clip_on=False,
)

fig.canvas.draw()
for ax, lab in zip(
    [axA, axB, axC, axD],
    "ABCD",
):
    bbox = ax.get_position()
    fig.text(
        bbox.x0 - 0.034,
        bbox.y1 + 0.013,
        lab,
        fontsize=25,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

pdf = (
    OUT /
    "FigureS3_v10_influenza_amplified_architecture.pdf"
)

png = (
    OUT /
    "FigureS3_v10_influenza_amplified_architecture.png"
)

svg = (
    OUT /
    "FigureS3_v10_influenza_amplified_architecture.svg"
)

fig.savefig(
    pdf,
    bbox_inches="tight",
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

print("=== FIGURE S3 v10 GENERATION COMPLETE ===")
print(pdf)
print(svg)
print(png)

print("\n=== FROZEN 170-GENE CLASS COUNTS ===")
print(
    x170["architecture_class"]
    .value_counts()
    .to_string()
)

print("\n=== DERIVATION METRICS ===")
for key in [
    "shared_GPL6884_total",
    "shared_GPL6884_crossplatform_evaluable",
    "shared_GPL6884_RSV_replicated",
    "shared_replicated_direct_not_significant",
    "shared_replicated_direct_significant",
    "shared_replicated_influenza_biased",
    "shared_replicated_RSV_biased",
]:
    print(f"{key}: {robust_map[key]}")


print(
    "\nNo signature classification or differential-expression "
    "analysis was recomputed."
)
