#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.home() / "Influenza_RSV_Project"

BASE = ROOT / "results/DS001_GSE38900/cell_composition"
OUTDIR = ROOT / "results/figure4_influenza_centered/figures"
PROVDIR = ROOT / "results/figure4_influenza_centered/provenance"

OUTDIR.mkdir(parents=True, exist_ok=True)
PROVDIR.mkdir(parents=True, exist_ok=True)

DISEASE = BASE / "tables/DS001_GSE38900_MCPcounter_disease_comparisons_final.tsv"
PATHWAY = BASE / "sensitivity/DS001_GSE38900_pathway_composition_classification.tsv"
GENES = BASE / "sensitivity/DS001_GSE38900_interferon_focus_gene_composition_effects.tsv"

for f in [DISEASE, PATHWAY, GENES]:
    if not f.exists():
        raise FileNotFoundError(f)

disease = pd.read_csv(DISEASE, sep="\t")
pathway = pd.read_csv(PATHWAY, sep="\t")
genes = pd.read_csv(GENES, sep="\t")

print("=== INPUT SCHEMA CHECK ===")
print("Disease columns:", list(disease.columns))
print("Pathway columns:", list(pathway.columns))
print("Gene columns:", list(genes.columns))
print()

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 16.0,
    "axes.titlesize": 22.0,
    "axes.labelsize": 19.0,
    "axes.labelweight": "bold",
    "xtick.labelsize": 16.0,
    "ytick.labelsize": 16.0,
    "legend.fontsize": 15.0,
    "axes.linewidth": 2.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

fig = plt.figure(figsize=(19.2, 13.4))

gs = fig.add_gridspec(
    2, 2,
    width_ratios=[1, 1],
    height_ratios=[1, 1],
    hspace=0.54,
    wspace=0.38,
    left=0.095,
    right=0.975,
    top=0.925,
    bottom=0.090
)

cell_order = [
    "Monocytic lineage",
    "Neutrophils",
    "B lineage",
    "T cells",
    "CD8 T cells",
    "Cytotoxic lymphocytes",
    "NK cells",
    "Myeloid dendritic cells",
    "Endothelial cells",
    "Fibroblasts",
]


def get_fdr_column(df):
    """
    Use the existing frozen multiple-testing column if present.
    No p-values or FDR values are recalculated.
    """
    candidates = [
        "model_FDR_all10",
        "model_FDR",
        "model_FDR_primary8",
        "FDR",
        "fdr",
        "adj_p",
        "padj",
        "p_adj",
        "qvalue",
        "q_value",
        "BH_FDR",
        "fdr_bh",
    ]

    for col in candidates:
        if col in df.columns:
            return col

    return None


def forest_panel(ax, data, contrast, xlabel, title, show_legend=False):
    x = data[data["contrast"].eq(contrast)].copy()

    x["cell_type"] = pd.Categorical(
        x["cell_type"],
        categories=cell_order,
        ordered=True
    )

    x = x.sort_values("cell_type", ascending=False).reset_index(drop=True)

    y = np.arange(len(x))

    for yi, row in zip(y, x.itertuples()):
        ax.plot(
            [row.hedges_g_ci_low, row.hedges_g_ci_high],
            [yi, yi],
            color="0.48",
            lw=2.8,
            zorder=1
        )

    primary = (
        x["primary_immune_population"]
        .astype(str)
        .str.upper()
        .eq("TRUE")
    )

    fdr_col = get_fdr_column(x)

    if fdr_col is not None:
        significant = pd.to_numeric(
            x[fdr_col],
            errors="coerce"
        ).lt(0.05)
    else:
        significant = pd.Series(False, index=x.index)

    # Primary immune populations
    for sig_state, face, label_suffix in [
        (True, None, "FDR < 0.05"),
        (False, "white", "FDR ≥ 0.05"),
    ]:
        mask = primary & significant.eq(sig_state)

        if mask.any():
            ax.scatter(
                x.loc[mask, "hedges_g"],
                y[mask.to_numpy()],
                s=115,
                marker="o",
                facecolors=face if face is not None else None,
                edgecolors="black",
                linewidth=0.85,
                zorder=3,
                label=(
                    f"Primary immune, {label_suffix}"
                    if show_legend else None
                )
            )

    # Non-primary estimates
    for sig_state, face, label_suffix in [
        (True, None, "FDR < 0.05"),
        (False, "white", "FDR ≥ 0.05"),
    ]:
        mask = (~primary) & significant.eq(sig_state)

        if mask.any():
            ax.scatter(
                x.loc[mask, "hedges_g"],
                y[mask.to_numpy()],
                s=105,
                marker="s",
                facecolors=face if face is not None else None,
                edgecolors="black",
                linewidth=0.85,
                zorder=3,
                label=(
                    f"Non-primary, {label_suffix}"
                    if show_legend else None
                )
            )

    # If the table has no recognized FDR column, retain the original
    # primary/non-primary encoding instead of inventing significance.
    if fdr_col is None:
        ax.collections.clear()

        ax.scatter(
            x.loc[primary, "hedges_g"],
            y[primary.to_numpy()],
            s=115,
            marker="o",
            edgecolor="black",
            linewidth=0.85,
            zorder=3,
            label="Primary immune population" if show_legend else None
        )

        ax.scatter(
            x.loc[~primary, "hedges_g"],
            y[(~primary).to_numpy()],
            s=105,
            marker="s",
            edgecolor="black",
            linewidth=0.85,
            zorder=3,
            label="Non-primary estimate" if show_legend else None
        )

    ax.axvline(0, color="0.18", lw=2.0)

    ax.set_yticks(y)
    ax.set_yticklabels(x["cell_type"], fontweight="bold")

    ax.set_xlabel(xlabel, fontweight="bold")
    # Panel title added later at figure level for exact alignment.

    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", width=1.8, length=6.5)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")

    if show_legend:
        leg = ax.legend(
            frameon=True,
            loc="upper left",
            bbox_to_anchor=(0.015, 0.985),
            fontsize=14.0,
            handletextpad=0.5,
            labelspacing=0.42,
            borderpad=0.45,
            facecolor="white",
            edgecolor="0.80",
            framealpha=0.94,
        )
        for t in leg.get_texts():
            t.set_fontweight("bold")

    return x, fdr_col


# ============================================================
# PANEL A
# ============================================================

axA = fig.add_subplot(gs[0, 0])

a, fdr_A = forest_panel(
    axA,
    disease,
    "GPL6884_Influenza_vs_control",
    "Hedges' g, influenza A vs control",
    "Influenza-associated leukocyte-state shifts",
    show_legend=False
)


# ============================================================
# PANEL B
# ============================================================

axB = fig.add_subplot(gs[0, 1])

b, fdr_B = forest_panel(
    axB,
    disease,
    "GPL6884_Influenza_vs_RSV",
    "Hedges' g, influenza A vs RSV",
    "Cellular differentiation of influenza from RSV",
    show_legend=False
)


# ============================================================
# PANEL C
# ============================================================

axC = fig.add_subplot(gs[1, 0])

c = pathway[
    pathway["contrast"].eq("GPL6884_Influenza_vs_control")
].copy()

domain_order = [
    "Interferon",
    "Antigen_Presentation",
    "Neutrophil_Monocyte",
    "Adaptive_Immunity",
    "Metabolism",
    "Translation",
]

domain_labels = {
    "Interferon": "Interferon",
    "Antigen_Presentation": "Antigen presentation",
    "Neutrophil_Monocyte": "Neutrophil / monocyte",
    "Adaptive_Immunity": "Adaptive immunity",
    "Metabolism": "Metabolism",
    "Translation": "Translation",
}

class_labels = {
    "Composition_robust_amplified": "Robust",
    "Composition_sensitive": "Composition-sensitive",
    "Unclassified_base_not_significant": "Base NS",
}

c["preregistered_domain"] = pd.Categorical(
    c["preregistered_domain"],
    categories=domain_order,
    ordered=True
)

c = c.sort_values(
    "preregistered_domain",
    ascending=False
).reset_index(drop=True)

yC = np.arange(len(c))

for yi, row in zip(yC, c.itertuples()):
    axC.plot(
        [row.beta_base, row.beta_adjusted],
        [yi, yi],
        color="0.50",
        lw=2.8,
        zorder=1
    )

axC.scatter(
    c["beta_base"],
    yC,
    s=115,
    marker="o",
    edgecolor="black",
    linewidth=0.85,
    zorder=3,
    label="Base model"
)

axC.scatter(
    c["beta_adjusted"],
    yC,
    s=120,
    marker="D",
    edgecolor="black",
    linewidth=0.85,
    zorder=3,
    label="Composition-adjusted"
)

axC.axvline(0, color="0.18", lw=2.0)

axC.set_yticks(yC)
axC.set_yticklabels(
    [domain_labels[str(x)] for x in c["preregistered_domain"]],
    fontweight="bold"
)

axC.set_xlabel("Disease-associated pathway effect (β)", fontweight="bold")
# Panel title added later at figure level for exact alignment.

# Short, readable classification labels
for yi, row in zip(yC, c.itertuples()):
    cls = str(row.composition_class)

    if cls == "Composition_robust_amplified":
        short = "Robust"
    elif cls == "Composition_sensitive":
        short = "Composition-sensitive"
    else:
        # Do not create a new biological classification for
        # base-nonsignificant domains during figure redesign.
        short = None

    if short == "Composition-sensitive":
        xmax = max(row.beta_base, row.beta_adjusted)

        axC.text(
            xmax + 0.055,
            yi,
            short,
            va="center",
            ha="left",
            fontsize=14.0,
            fontweight="bold",
            color="0.28",
            bbox=dict(
                boxstyle="round,pad=0.10",
                facecolor="white",
                edgecolor="none",
                alpha=0.92
            )
        )

xmin = min(
    -0.50,
    float(c[["beta_base", "beta_adjusted"]].min().min()) - 0.06
)

xmax = max(
    0.62,
    float(c[["beta_base", "beta_adjusted"]].max().max()) + 0.18
)

axC.set_xlim(xmin, xmax)

legC = axC.legend(
    frameon=True,
    loc="upper left",
    bbox_to_anchor=(0.015, 0.985),
    fontsize=14.5,
    facecolor="white",
    edgecolor="0.80",
    framealpha=0.94,
    borderpad=0.45,
    labelspacing=0.42,
)
for t in legC.get_texts():
    t.set_fontweight("bold")

axC.spines[["top", "right"]].set_visible(False)
axC.tick_params(axis="both", width=1.55, length=5.8)
for tick in axC.get_xticklabels() + axC.get_yticklabels():
    tick.set_fontweight("bold")


# ============================================================
# PANEL D
# ============================================================

axD = fig.add_subplot(gs[1, 1])

g = genes[
    genes["contrast"].eq("GPL6884_Influenza_vs_control")
].copy()

gene_order = [
    "ISG15",
    "OAS3",
    "OASL",
    "IRF7",
    "EIF2AK2",
    "IFIH1",
    "STAT1",
    "TRIM25",
]

g["gene_symbol"] = pd.Categorical(
    g["gene_symbol"],
    categories=gene_order,
    ordered=True
)

g = (
    g[g["gene_symbol"].notna()]
    .sort_values("gene_symbol", ascending=False)
    .reset_index(drop=True)
)

yD = np.arange(len(g))

for yi, row in zip(yD, g.itertuples()):
    axD.plot(
        [row.beta_base, row.beta_adjusted],
        [yi, yi],
        color="0.50",
        lw=2.8,
        zorder=1
    )

axD.scatter(
    g["beta_base"],
    yD,
    s=115,
    marker="o",
    edgecolor="black",
    linewidth=0.85,
    zorder=3,
    label="Base model"
)

axD.scatter(
    g["beta_adjusted"],
    yD,
    s=120,
    marker="D",
    edgecolor="black",
    linewidth=0.85,
    zorder=3,
    label="Composition-adjusted"
)

axD.axvline(0, color="0.18", lw=2.0)

axD.set_yticks(yD)
axD.set_yticklabels(
    g["gene_symbol"],
    fontweight="bold"
)

axD.set_xlabel("Gene-level disease effect (β)", fontweight="bold")
# Panel title added later at figure level for exact alignment.

# No significance stars in v2.
# All displayed genes remain FDR significant after adjustment;
# exact FDR values remain in the source table and manuscript text.

legD = axD.legend(
    frameon=False,
    loc="lower right",
    bbox_to_anchor=(0.98, 0.02),
    fontsize=14.5
)
for t in legD.get_texts():
    t.set_fontweight("bold")

axD.spines[["top", "right"]].set_visible(False)
axD.tick_params(axis="both", width=1.55, length=5.8)
for tick in axD.get_xticklabels() + axD.get_yticklabels():
    tick.set_fontweight("bold")


# ============================================================
# PANEL HEADERS — LOCKED FIGURE 2/3 SUBMISSION STYLE
# ============================================================

PANEL_LABEL_SIZE = 25.0
PANEL_TITLE_SIZE = 23.0

panel_titles = {
    "A": "Influenza-associated leukocyte-state shifts",
    "B": "Cellular differentiation of influenza from RSV",
    "C": "Composition robustness of influenza\nbiological programs",
    "D": "Interferon genes remain composition-robust",
}

# Use figure coordinates so each panel letter and title share an exact baseline.
fig.canvas.draw()
for label, axis in zip(["A", "B", "C", "D"], [axA, axB, axC, axD]):
    box = axis.get_position()
    y = box.y1 + (0.030 if label == "C" else 0.024)
    fig.text(
        box.x0 - 0.034,
        y,
        label,
        fontsize=PANEL_LABEL_SIZE,
        fontweight="bold",
        ha="left",
        va="bottom"
    )
    fig.text(
        box.x0,
        y,
        panel_titles[label],
        fontsize=PANEL_TITLE_SIZE,
        fontweight="bold",
        ha="left",
        va="bottom"
    )

# No overall Figure X title and no footer: manuscript legend carries detailed
# marker/model definitions. Panel A retains the compact A/B marker encoding key.

# ============================================================
# SAVE
# ============================================================

pdf = OUTDIR / "Figure4_SUBMISSION_v8_influenza_cellular_context.pdf"
png = OUTDIR / "Figure4_SUBMISSION_v8_influenza_cellular_context.png"
svg = OUTDIR / "Figure4_SUBMISSION_v8_influenza_cellular_context.svg"

fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=600, bbox_inches="tight")
fig.savefig(svg, format="svg", bbox_inches="tight")
plt.close(fig)


# ============================================================
# PROVENANCE
# ============================================================

manifest = PROVDIR / "Figure4_SUBMISSION_v8_input_manifest.txt"

with open(manifest, "w") as fh:
    fh.write(
        "FIGURE 4 SUBMISSION v8\n"
        "Cellular context and composition robustness of the influenza A host response\n\n"
        "STATUS: submission-style rendering; no scientific analysis rerun\n"
        "No scientific analysis rerun.\n\n"
        "Frozen inputs:\n"
        f"{DISEASE}\n"
        f"{PATHWAY}\n"
        f"{GENES}\n\n"
        "Primary contrasts:\n"
        "GPL6884_Influenza_vs_control\n"
        "GPL6884_Influenza_vs_RSV\n\n"
        "Submission v8 visual changes:\n"
        "1. Removed the internal Panel A legend.\n"
        "2. Marker encoding retained in the figure data but explained in the manuscript legend rather than inside Panel A.\n"
        "3. Removed repetitive Robust annotations from Panel C while retaining the frozen Composition-sensitive classification.\n"
        "4. Standardized Panel C title and panel-heading alignment with Panels A/B/D.\n"
        "5. Refined Panel D legend placement only; gene-level data unchanged.\n"
        "6. No values, contrasts, domains, genes, FDR calculations, or models changed.\n\n"
        f"Panel A FDR column detected: {fdr_A}\n"
        f"Panel B FDR column detected: {fdr_B}\n"
    )

print("Written:")
print(pdf)
print(png)
print(svg)
print(manifest)

print()
print("=== PANEL A FDR COLUMN ===")
print(fdr_A)

print()
print("=== PANEL B FDR COLUMN ===")
print(fdr_B)

print()
print("=== PANEL C ===")
print(
    c[
        [
            "preregistered_domain",
            "beta_base",
            "beta_adjusted",
            "FDR_adjusted",
            "composition_class"
        ]
    ].to_string(index=False)
)

print()
print("=== PANEL D ===")
print(
    g[
        [
            "gene_symbol",
            "beta_base",
            "beta_adjusted",
            "attenuation_percent",
            "FDR_adjusted"
        ]
    ].to_string(index=False)
)
