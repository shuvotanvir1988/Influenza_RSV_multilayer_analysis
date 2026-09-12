#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.home() / "Influenza_RSV_Project"

PATH = (
    ROOT /
    "results/DS001_GSE38900/pathway_analysis/tables"
)

BASE = (
    ROOT /
    "results/supplementary_figures/"
    "FigureS4_influenza_pathways"
)

OUT = BASE / "figures"
TAB = BASE / "tables"

OUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 18,
    "axes.titlesize": 23,
    "axes.titleweight": "bold",
    "axes.labelsize": 19,
    "axes.labelweight": "bold",
    "xtick.labelsize": 16.5,
    "ytick.labelsize": 16.5,
    "legend.fontsize": 15.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

ORANGE = "#D95F02"
GREY = "#7A7A7A"

THEME_COLORS = {
    "Interferon": "#D95F02",
    "Antigen_Presentation": "#7570B3",
    "Translation": "#666666",
    "Neutrophil_Myeloid": "#1B9E77",
    "Metabolism": "#A6761D",
}

# ============================================================
# Load frozen/pre-existing pathway outputs
# ============================================================

hall = pd.read_csv(
    PATH /
    "fgsea_hallmark/"
    "GPL6884_InfluenzaA_vs_control_Hallmark_fgsea.tsv",
    sep="\t",
)

react = pd.read_csv(
    PATH /
    "fgsea_reactome/"
    "GPL6884_InfluenzaA_vs_control_Reactome_fgsea.tsv",
    sep="\t",
)

le = pd.read_csv(
    PATH /
    "robustness_audit/"
    "DS001_GSE38900_corrected_leading_edge_gene_recurrence.tsv",
    sep="\t",
)

theme_summary = pd.read_csv(
    PATH /
    "robustness_audit/"
    "DS001_GSE38900_corrected_leading_edge_theme_summary.tsv",
    sep="\t",
)

robust = pd.read_csv(
    PATH /
    "robustness_audit/"
    "DS001_GSE38900_major_theme_robustness_summary.tsv",
    sep="\t",
)

# ============================================================
# Panel A source
# ============================================================

hall_sig = (
    hall[hall["padj"] < 0.05]
    .copy()
    .sort_values("NES", ascending=True)
)

if len(hall_sig) != 25:
    raise RuntimeError(
        f"Expected 25 significant Hallmark pathways; found {len(hall_sig)}."
    )

def clean_hallmark(x):
    return (
        x.replace("HALLMARK_", "")
         .replace("_", " ")
         .title()
         .replace("Tnfa", "TNFα")
         .replace("Il6", "IL6")
         .replace("Il2", "IL2")
         .replace("Jak", "JAK")
         .replace("Stat3", "STAT3")
         .replace("Stat5", "STAT5")
         .replace("G2M", "G2M")
         .replace("E2F", "E2F")
         .replace("Mtorc1", "mTORC1")
         .replace("Myc", "MYC")
         .replace("P53", "p53")
    )

hall_sig["display"] = hall_sig["pathway"].map(clean_hallmark)

hall_sig["display"] = (
    hall_sig["display"]
    .str.replace("Uv Response Up", "UV Response Up", regex=False)
    .str.replace(
        "TNFα Signaling Via Nfkb",
        "TNFα Signaling via NFκB",
        regex=False,
    )
)


hall_sig.to_csv(
    TAB / "FigureS4_panelA_significant_Hallmark.tsv",
    sep="\t",
    index=False,
)

# ============================================================
# Panel B source — representative extended Reactome programs
# ============================================================

representative = {
    "REACTOME_INTERFERON_GAMMA_SIGNALING":
        "IFN-γ signaling",

    "REACTOME_ANTIMICROBIAL_MECHANISM_OF_IFN_STIMULATED_GENES":
        "Antimicrobial ISG mechanisms",

    "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES":
        "Toll-like receptor cascades",

    "REACTOME_INFLAMMASOMES":
        "Inflammasomes",

    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION":
        "Antigen cross-presentation",

    "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC":
        "MHC-I peptide loading",

    "REACTOME_NEUTROPHIL_DEGRANULATION":
        "Neutrophil degranulation",

    "REACTOME_CELL_CYCLE_MITOTIC":
        "Mitotic cell cycle",

    "REACTOME_TRANSLATION":
        "Translation",

    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION":
        "Translation elongation",
}

react_B = react[
    react["pathway"].isin(representative)
].copy()

missing_B = set(representative) - set(react_B["pathway"])

if missing_B:
    raise RuntimeError(
        f"Missing expected Reactome pathways: {sorted(missing_B)}"
    )

if not (react_B["padj"] < 0.05).all():
    bad = react_B.loc[
        react_B["padj"] >= 0.05,
        ["pathway", "padj"]
    ]
    raise RuntimeError(
        "One or more representative Reactome pathways "
        f"are not significant:\n{bad}"
    )

react_B["display"] = react_B["pathway"].map(representative)

react_B = react_B.sort_values(
    "NES",
    ascending=True,
)

react_B.to_csv(
    TAB / "FigureS4_panelB_representative_Reactome.tsv",
    sep="\t",
    index=False,
)

# ============================================================
# Panel C source — corrected recurrent leading-edge genes
# ============================================================

le_flu = le[
    le["contrast"].eq("Influenza_GPL6884")
].copy()

themes_C = [
    "Interferon",
    "Antigen_Presentation",
    "Translation",
]

selected_C = []

for theme in themes_C:

    z = le_flu[
        le_flu["theme"].eq(theme)
    ].copy()

    z["abs_median_NES"] = z[
        "median_pathway_NES"
    ].abs()

    z = z.sort_values(
        [
            "pathways_containing_gene",
            "abs_median_NES",
            "leading_edge_gene",
        ],
        ascending=[False, False, True],
    ).head(5)

    selected_C.append(z)

panelC = pd.concat(
    selected_C,
    ignore_index=True,
)

if len(panelC) != 15:
    raise RuntimeError(
        f"Expected 15 recurrent genes for Panel C; found {len(panelC)}."
    )

panelC.to_csv(
    TAB / "FigureS4_panelC_recurrent_leading_edge_genes.tsv",
    sep="\t",
    index=False,
)

# ============================================================
# Panel D source — influenza-only theme robustness
# ============================================================

rob_flu = robust[
    robust["contrast"].eq("Influenza_GPL6884")
].copy()

theme_flu = theme_summary[
    theme_summary["contrast"].eq("Influenza_GPL6884")
].copy()

panelD = rob_flu.merge(
    theme_flu[
        [
            "theme",
            "maximum_pathway_recurrence",
        ]
    ],
    on="theme",
    how="left",
    validate="one_to_one",
)

theme_order = [
    "Interferon",
    "Antigen_Presentation",
    "Neutrophil_Myeloid",
    "Translation",
    "Metabolism",
]

panelD["theme"] = pd.Categorical(
    panelD["theme"],
    categories=theme_order,
    ordered=True,
)

panelD = panelD.sort_values("theme")

panelD["significant_fraction"] = (
    panelD["significant_FDR05"] /
    panelD["pathways"]
)

panelD.to_csv(
    TAB / "FigureS4_panelD_theme_robustness.tsv",
    sep="\t",
    index=False,
)

# ============================================================
# Figure
# ============================================================

fig = plt.figure(
    figsize=(21.5, 17.0)
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.105,
    right=0.975,
    bottom=0.085,
    top=0.955,
    wspace=0.46,
    hspace=0.46,
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

def panel_label(ax, lab):
    ax.text(
        -0.16,
        1.075,
        lab,
        transform=ax.transAxes,
        fontsize=26,
        fontweight="bold",
        ha="left",
        va="top",
    )

# ============================================================
# A — Hallmark landscape
# ============================================================

colors_A = [
    ORANGE if v > 0 else GREY
    for v in hall_sig["NES"]
]

yA = np.arange(len(hall_sig))

axA.hlines(
    yA,
    0,
    hall_sig["NES"],
    color=colors_A,
    linewidth=2.2,
    alpha=0.65,
)

axA.scatter(
    hall_sig["NES"],
    yA,
    c=colors_A,
    s=72,
    edgecolors="white",
    linewidths=0.45,
    zorder=3,
)

axA.axvline(
    0,
    color="black",
    linewidth=0.9,
    alpha=0.5,
)

axA.set_yticks(yA)

axA.set_yticklabels(
    hall_sig["display"],
    fontsize=13.2,
    fontweight="bold",
)


axA.set_xlabel(
    "Normalized enrichment score (NES)"
)

axA.set_title(
    "Significant Hallmark programs",
    loc="left",
    pad=12,
)

axA.grid(
    axis="x",
    alpha=0.15,
)

# ============================================================
# B — representative Reactome programs
# ============================================================

yB = np.arange(len(react_B))

colors_B = [
    ORANGE if v > 0 else GREY
    for v in react_B["NES"]
]

axB.barh(
    yB,
    react_B["NES"],
    color=colors_B,
    alpha=0.9,
)

axB.axvline(
    0,
    color="black",
    linewidth=0.9,
    alpha=0.5,
)

axB.set_yticks(yB)

axB.set_yticklabels(
    react_B["display"],
    fontsize=15.0,
    fontweight="bold",
)

axB.set_xlabel(
    "Normalized enrichment score (NES)"
)

axB.set_title(
    "Extended Reactome pathway architecture",
    loc="left",
    pad=12,
)

axB.grid(
    axis="x",
    alpha=0.15,
)

for i, row in react_B.reset_index(drop=True).iterrows():

    if row["NES"] >= 0:
        x = row["NES"] - 0.08
        ha = "right"
    else:
        x = row["NES"] + 0.08
        ha = "left"

    axB.text(
        x,
        i,
        f'{row["NES"]:.2f}',
        va="center",
        ha=ha,
        fontsize=13.5,
        fontweight="bold",
        color="white",
    )

# ============================================================
# C — recurrent leading-edge genes
# ============================================================

theme_display = {
    "Interferon": "Interferon",
    "Antigen_Presentation": "Antigen presentation",
    "Translation": "Translation",
}

panelC_plot = panelC.copy()

panelC_plot["display"] = (
    panelC_plot["leading_edge_gene"]
)

panelC_plot = panelC_plot.sort_values(
    [
        "theme",
        "pathways_containing_gene",
        "leading_edge_gene",
    ],
    ascending=[True, True, True],
)

yC = np.arange(len(panelC_plot))

axC.hlines(
    yC,
    0,
    panelC_plot["pathways_containing_gene"],
    color=[
        THEME_COLORS[t]
        for t in panelC_plot["theme"]
    ],
    linewidth=2.0,
    alpha=0.5,
)

axC.scatter(
    panelC_plot["pathways_containing_gene"],
    yC,
    s=92,
    c=[
        THEME_COLORS[t]
        for t in panelC_plot["theme"]
    ],
    edgecolors="white",
    linewidths=0.5,
)

axC.set_yticks(yC)

axC.set_yticklabels(
    panelC_plot["display"],
    fontsize=15.5,
    fontweight="bold",
)

axC.set_xlabel(
    "Pathways containing gene"
)

axC.set_title(
    "Recurrent influenza leading-edge genes",
    loc="left",
    pad=12,
)

axC.set_xlim(
    0,
    panelC_plot["pathways_containing_gene"].max() + 1.2,
)

axC.grid(
    axis="x",
    alpha=0.15,
)

from matplotlib.lines import Line2D

theme_handles = [
    Line2D(
        [0], [0],
        marker="o",
        linestyle="",
        markersize=7,
        markerfacecolor=THEME_COLORS["Antigen_Presentation"],
        markeredgecolor="white",
        label="Antigen presentation",
    ),
    Line2D(
        [0], [0],
        marker="o",
        linestyle="",
        markersize=7,
        markerfacecolor=THEME_COLORS["Interferon"],
        markeredgecolor="white",
        label="Interferon",
    ),
    Line2D(
        [0], [0],
        marker="o",
        linestyle="",
        markersize=7,
        markerfacecolor=THEME_COLORS["Translation"],
        markeredgecolor="white",
        label="Translation",
    ),
]

axC.legend(
    handles=theme_handles,
    loc="lower right",
    frameon=False,
    fontsize=14.0,
)

# ============================================================
# D — theme robustness
# ============================================================

display_D = {
    "Interferon": "Interferon",
    "Antigen_Presentation": "Antigen presentation",
    "Neutrophil_Myeloid": "Neutrophil / myeloid",
    "Translation": "Translation",
    "Metabolism": "Metabolism",
}

panelD_plot = panelD.copy()

panelD_plot["display"] = (
    panelD_plot["theme"]
    .astype(str)
    .map(display_D)
)

panelD_plot = panelD_plot.iloc[::-1].reset_index(drop=True)

yD = np.arange(len(panelD_plot))

bars = axD.barh(
    yD,
    panelD_plot["significant_fraction"],
    color=[
        THEME_COLORS[str(t)]
        for t in panelD_plot["theme"]
    ],
    alpha=0.88,
)

axD.set_yticks(yD)

axD.set_yticklabels(
    panelD_plot["display"],
    fontsize=15.5,
    fontweight="bold",
)

axD.set_xlim(0, 1.0)

axD.set_xlabel(
    "Fraction of pathways significant (FDR < 0.05)"
)

axD.set_title(
    "Preregistered-theme robustness",
    loc="left",
    pad=12,
)

axD.grid(
    axis="x",
    alpha=0.15,
)

for bar, (_, row) in zip(
    bars,
    panelD_plot.iterrows(),
):

    ymid = (
        bar.get_y()
        + bar.get_height()/2
    )

    axD.text(
        min(
            row["significant_fraction"] + 0.025,
            0.94,
        ),
        ymid + 0.10,
        (
            f'{int(row["significant_FDR05"])}'
            f'/{int(row["pathways"])}'
        ),
        va="center",
        ha="left",
        fontsize=14.5,
        fontweight="bold",
    )

    axD.text(
        min(
            row["significant_fraction"] + 0.025,
            0.94,
        ),
        ymid - 0.13,
        (
            f'max recurrence '
            f'{int(row["maximum_pathway_recurrence"])}×'
        ),
        va="center",
        ha="left",
        fontsize=12.8,
    )

for ax, lab in zip(
    [axA, axB, axC, axD],
    "ABCD",
):
    panel_label(ax, lab)

pdf = (
    OUT /
    "FigureS4_v4_extended_influenza_pathway_architecture.pdf"
)

png = (
    OUT /
    "FigureS4_v4_extended_influenza_pathway_architecture.png"
)

svg = (
    OUT /
    "FigureS4_v4_extended_influenza_pathway_architecture.svg"
)

fig.savefig(
    pdf,
    bbox_inches="tight",
)

fig.savefig(
    png,
    dpi=600,
    bbox_inches="tight",
)

fig.savefig(
    svg,
    bbox_inches="tight",
)

plt.close(fig)

# ============================================================
# Console QC
# ============================================================

print("=== FIGURE S4 GENERATION COMPLETE ===")

print(pdf)
print(png)
print(svg)

print("\n=== PANEL A ===")
print(
    f"Significant Hallmark pathways: {len(hall_sig)}"
)
print(
    f"Positive NES: {(hall_sig['NES'] > 0).sum()}"
)
print(
    f"Negative NES: {(hall_sig['NES'] < 0).sum()}"
)

print("\n=== PANEL B ===")
print(
    react_B[
        ["display", "NES", "padj"]
    ].to_string(index=False)
)

print("\n=== PANEL C ===")
print(
    panelC[
        [
            "theme",
            "leading_edge_gene",
            "pathways_containing_gene",
            "median_pathway_NES",
            "minimum_pathway_FDR",
        ]
    ].to_string(index=False)
)

print("\n=== PANEL D ===")
print(
    panelD[
        [
            "theme",
            "pathways",
            "significant_FDR05",
            "maximum_pathway_recurrence",
        ]
    ].to_string(index=False)
)

print(
    "\nNo FGSEA, leading-edge, or pathway classification "
    "analysis was rerun."
)

print(
    "Corrected leading-edge recurrence table was used."
)

print(
    "Figure S4 contains influenza-only pathway results."
)
