#!/usr/bin/env python3

from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


# ==========================================================
# PATHS
# ==========================================================

ROOT = Path.home() / "Influenza_RSV_Project"

OUT = (
    ROOT / "results/supplementary_figures/"
    "FigureS7_crispr_functional_convergence"
)

FIGDIR = OUT / "figures"
TABDIR = OUT / "tables"

FIGDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)

PDF = FIGDIR / "FigureS7_v9_crispr_functional_convergence.pdf"
PNG = FIGDIR / "FigureS7_v9_crispr_functional_convergence.png"
SVG = FIGDIR / "FigureS7_v9_crispr_functional_convergence.svg"


# ==========================================================
# INPUTS — EXISTING FROZEN ANALYSIS OUTPUTS ONLY
# ==========================================================

A = pd.read_csv(
    ROOT / "results/crispr_validation/analysis/"
    "CRISPR_RANK_ENRICHMENT_FAMILY_FDR_v1.0.tsv",
    sep="\t",
)

B = pd.read_csv(
    ROOT / "results/crispr_validation/pathways/"
    "CRISPR_FROZEN58_PATHWAY_RANK_ENRICHMENT_v1.0.tsv",
    sep="\t",
)

C = pd.read_csv(
    ROOT / "results/crispr_validation/pathways/"
    "CRISPR_FROZEN7_DOMAIN_DIRECT_GENESET_ENRICHMENT_v1.0.tsv",
    sep="\t",
)

D = pd.read_csv(
    ROOT / "results/crispr_validation/analysis/"
    "CRISPR_FROZEN170_OBJECTIVE_FUNCTIONAL_CANDIDATES_v1.0.tsv",
    sep="\t",
)


# ==========================================================
# FILTER EXISTING RESULTS
# ==========================================================

A = A[A["screen"] == "IAV_P01"].copy()

A_order = [
    "shared_170",
    "shared_core_130",
    "influenza_amplified_37",
]

A["gene_set"] = pd.Categorical(
    A["gene_set"],
    categories=A_order,
    ordered=True,
)
A = A.sort_values("gene_set")


B = B[
    (B["screen"] == "IAV_P01")
    & (B["FDR_within_screen"] < 0.05)
].copy()

# Display a representative set rather than overcrowding the
# figure with all 32 significant pathways.
#
# Selection is deterministic:
# top two most significant/highest-enrichment pathways within
# each preregistered domain, then highest enrichment overall
# until a maximum of 14 displayed pathways.
B = B.sort_values(
    ["preregistered_domain",
     "FDR_within_screen",
     "rank_enrichment_ratio"],
    ascending=[True, True, False],
)

selected = (
    B.groupby("preregistered_domain", observed=True)
    .head(2)
)

if len(selected) < 14:
    remaining = B.loc[~B.index.isin(selected.index)].sort_values(
        ["FDR_within_screen", "rank_enrichment_ratio"],
        ascending=[True, False],
    )
    selected = pd.concat(
        [selected, remaining.head(14 - len(selected))],
        axis=0,
    )

Bplot = (
    selected
    .drop_duplicates()
    .sort_values("rank_enrichment_ratio", ascending=True)
    .copy()
)


C = C[C["screen"] == "IAV_P01"].copy()

domain_order = [
    "Interferon",
    "Innate_Cytokine",
    "Antigen_Presentation",
    "Neutrophil_Monocyte",
    "Adaptive_Immunity",
    "Metabolism",
    "Translation",
]

C["preregistered_domain"] = pd.Categorical(
    C["preregistered_domain"],
    categories=domain_order,
    ordered=True,
)
C = C.sort_values("preregistered_domain")


D = D[
    D["functional_evidence_layers_n"] >= 2
].copy()

D = D.sort_values(
    ["functional_evidence_layers_n", "gene_symbol"],
    ascending=[False, True],
)


# ==========================================================
# WRITE EXACT FIGURE-SOURCE TABLES
# ==========================================================

A.to_csv(
    TABDIR / "FigureS7_panelA_IAV_signature_enrichment.tsv",
    sep="\t",
    index=False,
)

Bplot.to_csv(
    TABDIR / "FigureS7_panelB_IAV_significant_pathways_displayed.tsv",
    sep="\t",
    index=False,
)

C.to_csv(
    TABDIR / "FigureS7_panelC_IAV_domain_enrichment.tsv",
    sep="\t",
    index=False,
)

D.to_csv(
    TABDIR / "FigureS7_panelD_multilayer_gene_evidence.tsv",
    sep="\t",
    index=False,
)


# ==========================================================
# STYLE
# ==========================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 18,
    "axes.titlesize": 23,
    "axes.titleweight": "bold",
    "axes.labelsize": 19,
    "axes.labelweight": "bold",
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def panel_label(ax, letter):
    ax.text(
        -0.12, 1.08,
        letter,
        transform=ax.transAxes,
        fontsize=26,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def panel_title(ax, title):
    ax.text(
        0.0, 1.08,
        title,
        transform=ax.transAxes,
        fontsize=23,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def clean_pathway(x):
    x = x.replace("REACTOME_", "")
    x = x.replace("_", " ")
    x = x.title()

    replacements = {
        "Nf Kappab": "NF-κB",
        "Nf Kb": "NF-κB",
        "Rrna": "rRNA",
        "Nmd": "NMD",
        "SrP": "SRP",
        "Bcr": "BCR",
        "Fceri": "FcεRI",
        "Tnfr2": "TNFR2",
    }

    for a, b in replacements.items():
        x = x.replace(a, b)

    return "\n".join(
        textwrap.wrap(x, width=42)
    )


domain_display = {
    "Interferon": "Interferon",
    "Innate_Cytokine": "Innate cytokine",
    "Antigen_Presentation": "Antigen presentation",
    "Neutrophil_Monocyte": "Neutrophil / monocyte",
    "Adaptive_Immunity": "Adaptive immunity",
    "Metabolism": "Metabolism",
    "Translation": "Translation",
}


# ==========================================================
# FIGURE GEOMETRY
# ==========================================================

fig = plt.figure(figsize=(26.0, 21.5))

gs = fig.add_gridspec(
    2, 2,
    width_ratios=[0.88, 1.72],
    height_ratios=[1.10, 1.05],
    left=0.095,
    right=0.980,
    bottom=0.090,
    top=0.945,
    wspace=0.82,
    hspace=0.70,
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])


# ==========================================================
# PANEL A — FROZEN SIGNATURE ENRICHMENT
# ==========================================================

A_labels = {
    "shared_170": "Shared 170",
    "shared_core_130": "Shared core 130",
    "influenza_amplified_37": "Influenza-amplified 37",
}

y = np.arange(len(A))

axA.axvline(
    1.0,
    linestyle="--",
    linewidth=1.3,
    color="0.45",
)

axA.scatter(
    A["rank_enrichment_ratio"],
    y,
    s=120,
    edgecolor="black",
    linewidth=0.6,
    zorder=3,
)

for i, (_, r) in enumerate(A.iterrows()):
    axA.text(
        r["rank_enrichment_ratio"] + 0.015,
        i,
        f"FDR={r['family_FDR']:.2f}",
        fontsize=17,
        fontweight="bold",
        va="center",
        ha="left",
    )

axA.set_yticks(y)
axA.set_yticklabels(
    [A_labels[str(x)] for x in A["gene_set"]],
    fontsize=18,
    fontweight="bold",
)

axA.invert_yaxis()
axA.set_xlim(0.75, 1.30)

axA.set_xlabel(
    "CRISPR rank-enrichment ratio\n(>1 = stronger dependency enrichment)"
)

panel_label(axA, "A")
panel_title(
    axA,
    "Frozen transcriptional signatures",
)

axA.text(
    0.0,
    1.015,
    "No prespecified signature reached family-level FDR < 0.05",
    transform=axA.transAxes,
    fontsize=16.5,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)


# ==========================================================
# PANEL B — PATHWAY-LEVEL CONVERGENCE
# ==========================================================

yb = np.arange(len(Bplot))

fdr_score = -np.log10(
    Bplot["FDR_within_screen"].clip(lower=1e-6)
)

sc = axB.scatter(
    Bplot["rank_enrichment_ratio"],
    yb,
    s=75 + 18 * fdr_score,
    c=fdr_score,
    cmap="viridis",
    edgecolor="black",
    linewidth=0.45,
    zorder=3,
)

axB.axvline(
    1.0,
    linestyle="--",
    linewidth=1.2,
    color="0.45",
)

axB.set_yticks(yb)
axB.set_yticklabels(
    [clean_pathway(x) for x in Bplot["pathway"]],
    fontsize=16.0,
    fontweight="bold",
)

axB.set_xlabel(
    "IAV CRISPR rank-enrichment ratio"
)

panel_label(axB, "B")
panel_title(
    axB,
    "Pathway-level functional convergence",
)

axB.text(
    0.0,
    1.015,
    "32 / 56 pathways FDR < 0.05",
    transform=axB.transAxes,
    fontsize=16.5,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)

cb = fig.colorbar(
    sc,
    ax=axB,
    fraction=0.032,
    pad=0.020,
)

cb.set_label(
    "−log10(FDR)",
    fontsize=17.0,
    fontweight="bold",
)

cb.ax.tick_params(labelsize=15.5)


# ==========================================================
# PANEL C — SEVEN PRESPECIFIED DOMAINS
# ==========================================================

yc = np.arange(len(C))

sig = C["FDR_within_screen"] < 0.05

axC.axvline(
    1.0,
    linestyle="--",
    linewidth=1.3,
    color="0.45",
)

axC.scatter(
    C.loc[~sig, "rank_enrichment_ratio"],
    yc[~sig],
    s=110,
    facecolor="white",
    edgecolor="0.45",
    linewidth=1.2,
    zorder=3,
)

axC.scatter(
    C.loc[sig, "rank_enrichment_ratio"],
    yc[sig],
    s=120,
    edgecolor="black",
    linewidth=0.55,
    zorder=3,
)

for i, (_, r) in enumerate(C.iterrows()):
    fdr = r["FDR_within_screen"]

    if fdr < 0.001:
        txt = "FDR<0.001"
    else:
        txt = f"FDR={fdr:.3f}"

    axC.text(
        r["rank_enrichment_ratio"] + 0.018,
        i,
        txt,
        fontsize=16.0,
        fontweight="bold",
        va="center",
        ha="left",
    )

axC.set_yticks(yc)
axC.set_yticklabels(
    [
        domain_display[str(x)]
        for x in C["preregistered_domain"]
    ],
    fontsize=18,
    fontweight="bold",
)

axC.invert_yaxis()
axC.set_xlim(
    0.95,
    max(1.62, C["rank_enrichment_ratio"].max() + 0.12)
)

axC.set_xlabel(
    "Direct domain gene-set rank-enrichment ratio"
)

panel_label(axC, "C")
panel_title(
    axC,
    "Prespecified biological domains",
)

axC.text(
    0.0,
    1.015,
    "6 / 7 domains FDR < 0.05",
    transform=axC.transAxes,
    fontsize=16.5,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)


# ==========================================================
# PANEL D — MULTILAYER GENE EVIDENCE
# ==========================================================

genes = D["gene_symbol"].tolist()

evidence_cols = [
    "iav_functional_support",
    "rsv_viability_positive_support",
    "rsv_early_positive_support",
    "rsv_viability_negative_support",
    "rsv_early_negative_support",
]

evidence_labels = [
    "IAV\nCRISPR",
    "RSV viability\npositive",
    "RSV early\npositive",
    "RSV viability\nnegative",
    "RSV early\nnegative",
]


def evidence_level(value):
    v = str(value)

    if (
        "FDR_lt_0.05" in v
        or "top1pct" in v
    ):
        return 3

    if "top5pct" in v:
        return 2

    if "top10pct" in v:
        return 1

    return 0


M = np.array([
    [
        evidence_level(row[c])
        for c in evidence_cols
    ]
    for _, row in D.iterrows()
])

# Use neutral sequential values; categories are additionally
# encoded with text/symbols so interpretation is not dependent
# on color alone.
im = axD.imshow(
    M,
    aspect="auto",
    interpolation="nearest",
    vmin=0,
    vmax=3,
    cmap="Blues",
)

axD.set_xticks(np.arange(len(evidence_labels)))
axD.set_xticklabels(
    evidence_labels,
    rotation=0,
    ha="center",
    fontsize=14.0,
    fontweight="bold",
    linespacing=1.25,
)

axD.set_yticks(np.arange(len(genes)))
axD.set_yticklabels(
    genes,
    fontsize=18,
    fontweight="bold",
)

for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        val = M[i, j]

        symbol = {
            0: "—",
            1: "10%",
            2: "5%",
            3: "HC",
        }[val]

        axD.text(
            j,
            i,
            symbol,
            ha="center",
            va="center",
            fontsize=15.0,
            fontweight="bold" if val > 0 else "normal",
        )

# Architecture class and influenza transcriptional effect
# shown as a sixth text-only column immediately adjacent to
# the five CRISPR evidence columns.

annotation_x = len(evidence_cols) + 0.18

# Compact evidence key and separate annotation-column header.

axD.text(
    0.0,
    1.015,
    "HC = top 1% / targeted FDR < 0.05   ·   "
    "5% = top 5%   ·   — = no displayed support",
    transform=axD.transAxes,
    fontsize=13.0,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)

axD.text(
    1.04,
    1.015,
    "Architecture / influenza logFC",
    transform=axD.transAxes,
    fontsize=14.0,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)

for i, (_, r) in enumerate(D.iterrows()):
    cls = str(r["proteomics_architecture_class"]).replace(
        "_", " "
    )

    axD.text(
        annotation_x,
        i,
        f"{cls} / {r['flu_logFC']:+.2f}",
        fontsize=16.0,
        fontweight="bold",
        va="center",
        ha="left",
        clip_on=False,
    )

axD.set_xlim(-0.5, len(evidence_cols) + 2.75)

axD.text(
    -0.12,
    1.08,
    "D",
    transform=axD.transAxes,
    fontsize=26,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)
axD.text(
    0.0,
    1.08,
    "Genes supported across multiple CRISPR evidence layers",
    transform=axD.transAxes,
    fontsize=23,
    fontweight="bold",
    ha="left",
    va="bottom",
    clip_on=False,
)



# ==========================================================
# SAVE
# ==========================================================

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


# ==========================================================
# QC REPORT
# ==========================================================

print("=== FIGURE S7 V9 GENERATION COMPLETE ===")
print(PDF)
print(PNG)
print(SVG)

print("\n=== PANEL A ===")
print(
    A[
        [
            "gene_set",
            "covered_n",
            "rank_enrichment_ratio",
            "family_FDR",
        ]
    ].to_string(index=False)
)

print("\n=== PANEL B ===")
print("IAV pathways tested: 56")
print("IAV pathways FDR < 0.05: 32")
print("Displayed:", len(Bplot))

print("\n=== PANEL C ===")
print(
    C[
        [
            "preregistered_domain",
            "evaluable_genes_n",
            "rank_enrichment_ratio",
            "FDR_within_screen",
        ]
    ].to_string(index=False)
)

print("\n=== PANEL D ===")
print(
    D[
        [
            "gene_symbol",
            "proteomics_architecture_class",
            "functional_evidence_layers_n",
            "iav_functional_support",
            "rsv_viability_positive_support",
            "rsv_early_positive_support",
            "rsv_viability_negative_support",
            "rsv_early_negative_support",
        ]
    ].to_string(index=False)
)

print(
    "\nNo CRISPR screen was rerun.\n"
    "No permutation analysis was rerun.\n"
    "No enrichment analysis was rerun.\n"
    "Figure S7 uses existing frozen CRISPR outputs only."
)
