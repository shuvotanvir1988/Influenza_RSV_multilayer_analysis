#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.home() / "Influenza_RSV_Project"

SENS = (
    ROOT / "results/DS001_GSE38900/"
    "cell_composition/sensitivity"
)

BASE = (
    ROOT / "results/supplementary_figures/"
    "FigureS5_composition_robustness"
)

FIG = BASE / "figures"
TAB = BASE / "tables"

FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------
# Global style — consistent with manuscript figure standard
# ---------------------------------------------------------

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 18,
    "axes.titlesize": 23,
    "axes.titleweight": "bold",
    "axes.labelsize": 19,
    "axes.labelweight": "bold",
    "xtick.labelsize": 16.5,
    "ytick.labelsize": 16.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ---------------------------------------------------------
# Load frozen existing outputs
# ---------------------------------------------------------

pcvar = pd.read_csv(
    SENS / "DS001_GSE38900_MCPcounter_composition_PC_variance.tsv",
    sep="\t",
)

pathway = pd.read_csv(
    SENS / "DS001_GSE38900_composition_adjusted_pathway_effects.tsv",
    sep="\t",
)

genes = pd.read_csv(
    SENS / "DS001_GSE38900_interferon_focus_gene_composition_effects.tsv",
    sep="\t",
)

signature = pd.read_csv(
    SENS / "DS001_GSE38900_signature_composition_adjusted_effects.tsv",
    sep="\t",
)

# Influenza-centered manuscript-facing subsets only.
A = pcvar.loc[
    pcvar["platform"].eq("GPL6884")
].copy()

B = pathway.loc[
    pathway["contrast"].eq("GPL6884_Influenza_vs_control")
].copy()

C = genes.loc[
    genes["contrast"].eq("GPL6884_Influenza_vs_control")
].copy()

D = signature.loc[
    signature["contrast"].eq("GPL6884_Influenza_vs_control")
].copy()

# Explicit display orders.
domain_order = [
    "Interferon",
    "Antigen_Presentation",
    "Neutrophil_Monocyte",
    "Adaptive_Immunity",
    "Metabolism",
    "Translation",
]

gene_order = [
    "ISG15",
    "IRF7",
    "IFIH1",
    "STAT1",
    "EIF2AK2",
    "OAS3",
    "OASL",
    "TRIM25",
]

signature_order = [
    "Shared_170",
    "Shared_core_130",
    "Interferon_focus_8",
]

B["preregistered_domain"] = pd.Categorical(
    B["preregistered_domain"],
    categories=domain_order,
    ordered=True,
)
B = B.sort_values("preregistered_domain")

C["gene_symbol"] = pd.Categorical(
    C["gene_symbol"],
    categories=gene_order,
    ordered=True,
)
C = C.sort_values("gene_symbol")

D["signature"] = pd.Categorical(
    D["signature"],
    categories=signature_order,
    ordered=True,
)
D = D.sort_values("signature")

# Save exact figure source tables.
A.to_csv(
    TAB / "FigureS5_panelA_GPL6884_PC_variance.tsv",
    sep="\t", index=False
)
B.to_csv(
    TAB / "FigureS5_panelB_influenza_pathway_adjustment.tsv",
    sep="\t", index=False
)
C.to_csv(
    TAB / "FigureS5_panelC_influenza_gene_adjustment.tsv",
    sep="\t", index=False
)
D.to_csv(
    TAB / "FigureS5_panelD_influenza_signature_adjustment.tsv",
    sep="\t", index=False
)

# ---------------------------------------------------------
# Figure
# ---------------------------------------------------------

fig = plt.figure(figsize=(21.0, 16.5))

gs = fig.add_gridspec(
    2, 2,
    left=0.085,
    right=0.975,
    bottom=0.085,
    top=0.955,
    wspace=0.42,
    hspace=0.48,
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

# ---------------------------------------------------------
# Panel A
# ---------------------------------------------------------

x = np.arange(1, len(A) + 1)

axA.bar(
    x,
    A["variance_explained"] * 100,
    alpha=0.78,
)

axA.plot(
    x,
    A["cumulative_variance"] * 100,
    marker="o",
    linewidth=2.2,
)

axA.axvline(
    5.5,
    linestyle="--",
    linewidth=1.4,
    alpha=0.65,
)

axA.text(
    5.35,
    96,
    "PC1–PC5 used\nfor adjustment",
    ha="right",
    va="top",
    fontsize=15.5,
    fontweight="bold",
)

axA.set_xticks(x)
axA.set_xticklabels(A["PC"])
axA.set_ylim(0, 105)

axA.set_ylabel("Variance explained (%)")
axA.set_xlabel("Composition principal component")
axA.text(
    0.0,
    1.145,
    "Estimated leukocyte-composition structure",
    transform=axA.transAxes,
    fontsize=23,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)

# ---------------------------------------------------------
# Utility for identity plots
# ---------------------------------------------------------

def identity_plot(
    ax,
    df,
    xcol,
    ycol,
    labels,
    title,
    label_offsets=None,
):
    vals = np.concatenate([
        df[xcol].to_numpy(float),
        df[ycol].to_numpy(float),
    ])

    lo = np.nanmin(vals)
    hi = np.nanmax(vals)
    pad = max((hi - lo) * 0.18, 0.12)

    lo -= pad
    hi += pad

    ax.plot(
        [lo, hi],
        [lo, hi],
        linestyle="--",
        linewidth=1.5,
        alpha=0.65,
    )

    ax.scatter(
        df[xcol],
        df[ycol],
        s=100,
        alpha=0.86,
        edgecolor="black",
        linewidth=0.5,
    )

    if label_offsets is None:
        label_offsets = {}

    for _, r in df.iterrows():
        label = str(r[labels])

        dx, dy, ha, va = label_offsets.get(
            label,
            (6, 5, "left", "bottom"),
        )

        arrow = None

        if label in {"OAS3", "ISG15", "OASL"}:
            arrow = dict(
                arrowstyle="-",
                lw=0.8,
                color="black",
                alpha=0.7,
                shrinkA=2,
                shrinkB=3,
            )

        ax.annotate(
            label,
            (r[xcol], r[ycol]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=14.0,
            fontweight="bold",
            ha=ha,
            va=va,
            arrowprops=arrow,
        )

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)

    ax.axhline(0, linewidth=0.8, alpha=0.35)
    ax.axvline(0, linewidth=0.8, alpha=0.35)

    ax.set_xlabel("Unadjusted effect (β)")
    ax.set_ylabel("Composition-adjusted effect (β)")
    ax.set_title(title, loc="left")


# ---------------------------------------------------------
# Panel B
# ---------------------------------------------------------

B["display"] = (
    B["preregistered_domain"]
    .astype(str)
    .str.replace("_", " ", regex=False)
)

B["display"] = B["display"].replace({
    "Neutrophil Monocyte": "Neutrophil / monocyte",
})

B_offsets = {
    "Interferon": (7, 7, "left", "bottom"),
    "Antigen Presentation": (7, 5, "left", "bottom"),
    "Neutrophil / monocyte": (8, -2, "left", "center"),
    "Adaptive Immunity": (-8, 8, "right", "bottom"),
    "Metabolism": (7, 5, "left", "bottom"),
    "Translation": (7, 4, "left", "bottom"),
}

identity_plot(
    axB,
    B,
    "beta_base",
    "beta_adjusted",
    "display",
    "Pathway effects after composition adjustment",
    label_offsets=B_offsets,
)

# Replace automatic Panel B title with explicitly positioned text.
axB.set_title("", loc="left")

axB.text(
    0.0,
    1.145,
    "Pathway effects after composition adjustment",
    transform=axB.transAxes,
    fontsize=23,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)

# ---------------------------------------------------------
# Panel C
# ---------------------------------------------------------

C_offsets = {
    "ISG15": (24, -18, "left", "top"),
    "IRF7": (8, -12, "left", "top"),
    "IFIH1": (8, -2, "left", "center"),
    "STAT1": (-8, 10, "right", "bottom"),
    "EIF2AK2": (8, 10, "left", "bottom"),
    "OAS3": (16, 26, "left", "bottom"),
    "OASL": (-20, -14, "right", "top"),
    "TRIM25": (8, 6, "left", "bottom"),
}

identity_plot(
    axC,
    C,
    "beta_base",
    "beta_adjusted",
    "gene_symbol",
    "Antiviral-gene effects remain robust",
    label_offsets=C_offsets,
)

# Explicit title positioning ensures C and D align exactly.
axC.set_title(
    "Antiviral-gene effects remain robust",
    loc="left",
    y=1.075,
    pad=0,
)

# ---------------------------------------------------------
# Panel D
# ---------------------------------------------------------

display_map = {
    "Shared_170": "Replicated 170-gene\nresponse",
    "Shared_core_130": "Shared core\n(130 genes)",
    "Interferon_focus_8": "Interferon focus\n(8 genes)",
}

D["display"] = D["signature"].astype(str).map(display_map)

y = np.arange(len(D))

axD.scatter(
    D["beta_base"],
    y + 0.11,
    s=110,
    label="Unadjusted",
)

axD.scatter(
    D["beta_adjusted"],
    y - 0.11,
    s=110,
    marker="s",
    label="Composition adjusted",
)

for i, (_, r) in enumerate(D.iterrows()):
    axD.plot(
        [r["beta_base"], r["beta_adjusted"]],
        [i + 0.11, i - 0.11],
        linewidth=1.5,
        alpha=0.55,
    )

axD.axvline(
    0,
    linestyle="--",
    linewidth=1.1,
    alpha=0.55,
)

axD.set_yticks(y)
axD.set_yticklabels(
    D["display"],
    fontweight="bold",
)

axD.invert_yaxis()

axD.set_xlabel("Influenza vs control effect (β)")
axD.set_title(
    "Influenza signatures remain composition robust",
    loc="left",
    y=1.075,
    pad=0,
)

from matplotlib.lines import Line2D

panelD_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="",
        markersize=8,
        markerfacecolor="C0",
        markeredgecolor="C0",
        label="Unadjusted",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        linestyle="",
        markersize=8,
        markerfacecolor="C1",
        markeredgecolor="C1",
        label="Composition adjusted",
    ),
]

axD.legend(
    handles=panelD_handles,
    loc="upper left",
    bbox_to_anchor=(0.035, 0.965),
    frameon=True,
    framealpha=0.95,
    edgecolor="0.80",
    fontsize=14.5,
    borderpad=0.6,
    labelspacing=0.6,
    handletextpad=0.6,
)


# ---------------------------------------------------------
# Panel labels
# ---------------------------------------------------------

for label, ax in zip(
    ["A", "B", "C", "D"],
    [axA, axB, axC, axD],
):
    ax.text(
        -0.13,
        1.075,
        label,
        transform=ax.transAxes,
        fontsize=26,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

# ---------------------------------------------------------
# Output
# ---------------------------------------------------------

pdf = FIG / "FigureS5_v14_influenza_composition_robustness.pdf"
png = FIG / "FigureS5_v14_influenza_composition_robustness.png"
svg = FIG / "FigureS5_v14_influenza_composition_robustness.svg"

fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=600, bbox_inches="tight")
fig.savefig(svg, bbox_inches="tight")
plt.close(fig)

print("=== FIGURE S5 GENERATION COMPLETE ===")
print(pdf)
print(png)
print(svg)

print()
print("=== PANEL A PC VARIANCE ===")
print(
    A[
        ["PC", "variance_explained", "cumulative_variance"]
    ].to_string(index=False)
)

print()
print("=== PANEL B INFLUENZA PATHWAY EFFECTS ===")
print(
    B[
        [
            "preregistered_domain",
            "beta_base",
            "beta_adjusted",
            "FDR_base",
            "FDR_adjusted",
            "sign_preserved",
        ]
    ].to_string(index=False)
)

print()
print("=== PANEL C INFLUENZA ANTIVIRAL GENES ===")
print(
    C[
        [
            "gene_symbol",
            "beta_base",
            "beta_adjusted",
            "FDR_base",
            "FDR_adjusted",
            "sign_preserved",
        ]
    ].to_string(index=False)
)

print()
print("=== PANEL D INFLUENZA SIGNATURES ===")
print(
    D[
        [
            "signature",
            "beta_base",
            "beta_adjusted",
            "FDR_base",
            "FDR_adjusted",
            "sign_preserved",
        ]
    ].to_string(index=False)
)

print()
print("No MCP-counter estimation, pathway model, gene model,")
print("signature model, or differential-expression analysis was rerun.")
print("Figure S5 uses existing prespecified sensitivity outputs only.")
