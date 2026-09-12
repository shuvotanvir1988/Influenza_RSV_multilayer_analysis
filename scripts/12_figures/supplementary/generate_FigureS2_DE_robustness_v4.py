#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.home() / "Influenza_RSV_Project"

DE = (
    ROOT / "results/DS001_GSE38900/"
    "differential_expression/tables"
)

SENS = (
    ROOT / "results/DS001_GSE38900/"
    "differential_expression/sensitivity"
)

BASE = (
    ROOT / "results/supplementary_figures/"
    "FigureS2_DE_robustness"
)

OUT = BASE / "figures"
TAB = BASE / "tables"

OUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Session 42 — Supplementary Figure S2 v4
# Surgical Panel-D refinement of Session 42 v3.
#
# Scientific content unchanged:
#   - no DE models rerun
#   - same input DE/sensitivity tables
#   - same 25,440-gene merges in A-C
#   - same Pearson/Spearman calculations
#
# v4 changes:
#   - A-C retained from v3
#   - wider Panel D plotting region
#   - greater vertical separation among Panel D sensitivity rows
#   - cleaner value-label placement
#   - legend moved away from data
#   - more bottom breathing room
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
    "legend.fontsize": 16,
    "axes.linewidth": 1.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def read_table(path):
    x = pd.read_csv(path, sep="\t")
    required = {"gene_symbol", "logFC"}
    missing = required - set(x.columns)

    if missing:
        raise RuntimeError(
            f"{path.name}: missing required columns {missing}"
        )

    if x["gene_symbol"].duplicated().any():
        dup = x.loc[
            x["gene_symbol"].duplicated(False),
            "gene_symbol"
        ].head(10).tolist()

        raise RuntimeError(
            f"{path.name}: gene_symbol is not unique. "
            f"Examples: {dup}"
        )

    return x


def merge_effects(primary_path, sensitivity_path,
                  primary_name, sensitivity_name):

    a = read_table(primary_path)[
        ["gene_symbol", "logFC"]
    ].rename(columns={"logFC": primary_name})

    b = read_table(sensitivity_path)[
        ["gene_symbol", "logFC"]
    ].rename(columns={"logFC": sensitivity_name})

    m = a.merge(
        b,
        on="gene_symbol",
        how="inner",
        validate="one_to_one",
    )

    if len(m) == 0:
        raise RuntimeError(
            f"No genes merged between {primary_path.name} "
            f"and {sensitivity_path.name}"
        )

    return m


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


def scatter_concordance(
    ax,
    df,
    xcol,
    ycol,
    xlabel,
    ylabel,
    title,
):
    x = df[xcol].to_numpy(float)
    y = df[ycol].to_numpy(float)

    pearson = np.corrcoef(x, y)[0, 1]
    spearman = pd.Series(x).corr(
        pd.Series(y),
        method="spearman",
    )

    lim = max(
        abs(np.nanmin(x)),
        abs(np.nanmax(x)),
        abs(np.nanmin(y)),
        abs(np.nanmax(y)),
    ) * 1.05

    ax.scatter(
        x,
        y,
        s=15,
        alpha=0.18,
        linewidths=0,
        rasterized=True,
        zorder=2,
    )

    ax.plot(
        [-lim, lim],
        [-lim, lim],
        linestyle="--",
        linewidth=1.6,
        color="black",
        alpha=0.55,
        zorder=3,
    )

    ax.axhline(
        0,
        linewidth=1.0,
        color="#8A8A8A",
        alpha=0.28,
        zorder=1,
    )

    ax.axvline(
        0,
        linewidth=1.0,
        color="#8A8A8A",
        alpha=0.28,
        zorder=1,
    )

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel, labelpad=10)
    ax.set_ylabel(ylabel, labelpad=10)
    ax.set_title(title, loc="left", pad=16)

    ax.text(
        0.045,
        0.945,
        (
            f"Pearson r = {pearson:.3f}\n"
            f"Spearman ρ = {spearman:.3f}\n"
            f"Genes = {len(df):,}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
        linespacing=1.18,
    )

    style_axes(ax)
    return pearson, spearman


# ---------------------------------------------------------
# Panels A-C source data
# ---------------------------------------------------------

A = merge_effects(
    DE / "GPL6884_InfluenzaA_vs_control_full_DE.tsv",
    SENS / "GPL6884_unadjusted_Influenza_vs_control.tsv",
    "adjusted_logFC",
    "unadjusted_logFC",
)

B = merge_effects(
    DE / "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv",
    SENS / "GPL6884_unadjusted_Influenza_vs_RSV.tsv",
    "adjusted_logFC",
    "unadjusted_logFC",
)

C = merge_effects(
    DE / "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv",
    SENS / "GPL6884_age_restricted_Influenza_vs_RSV.tsv",
    "primary_logFC",
    "age_restricted_logFC",
)

A.to_csv(
    TAB / "FigureS2_panelA_Influenza_adjusted_vs_unadjusted.tsv",
    sep="\t",
    index=False,
)

B.to_csv(
    TAB / "FigureS2_panelB_Influenza_vs_RSV_adjusted_vs_unadjusted.tsv",
    sep="\t",
    index=False,
)

C.to_csv(
    TAB / "FigureS2_panelC_Influenza_vs_RSV_age_sensitivity.tsv",
    sep="\t",
    index=False,
)


# ---------------------------------------------------------
# Panel D source data
# ---------------------------------------------------------

D = pd.read_csv(
    SENS / "DS001_DE_sensitivity_summary.tsv",
    sep="\t",
)

required_D = {
    "analysis",
    "pearson_logFC",
    "spearman_logFC",
}

if not required_D.issubset(D.columns):
    raise RuntimeError(
        "Unexpected sensitivity summary schema: "
        f"{D.columns.tolist()}"
    )

D.to_csv(
    TAB / "FigureS2_panelD_global_sensitivity_summary.tsv",
    sep="\t",
    index=False,
)


# ---------------------------------------------------------
# Figure layout
# ---------------------------------------------------------

fig = plt.figure(figsize=(21.5, 15.2))

# Slightly asymmetric width allocation gives Panel D more room,
# while keeping Panels A-C visually unchanged in scale/style.
gs = fig.add_gridspec(
    2,
    2,
    left=0.070,
    right=0.982,
    bottom=0.090,
    top=0.945,
    wspace=0.48,
    hspace=0.50,
    width_ratios=[1.00, 1.13],
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])


scatter_concordance(
    axA,
    A,
    "adjusted_logFC",
    "unadjusted_logFC",
    "Adjusted logFC",
    "Unadjusted logFC",
    "Influenza vs control: covariate sensitivity",
)

scatter_concordance(
    axB,
    B,
    "adjusted_logFC",
    "unadjusted_logFC",
    "Adjusted logFC",
    "Unadjusted logFC",
    "Influenza vs RSV: covariate sensitivity",
)

scatter_concordance(
    axC,
    C,
    "primary_logFC",
    "age_restricted_logFC",
    "Primary-model logFC",
    "Age-restricted logFC",
    "Influenza vs RSV: age sensitivity",
)


# ---------------------------------------------------------
# Panel D — refined layout
# ---------------------------------------------------------

label_map = {
    "GPL10558 primary vs sex-adjusted":
        "GPL10558 RSV:\nsex adjustment",

    "GPL6884 RSV adjusted vs unadjusted":
        "GPL6884 RSV:\ncovariate adjustment",

    "GPL6884 Influenza adjusted vs unadjusted":
        "Influenza:\ncovariate adjustment",

    "GPL6884 Influenza-vs-RSV adjusted vs unadjusted":
        "Influenza vs RSV:\ncovariate adjustment",

    "GPL6884 Influenza-vs-RSV primary vs age-restricted":
        "Influenza vs RSV:\nage restriction",
}

Dplot = D.copy()
Dplot["display"] = Dplot["analysis"].map(label_map)

if Dplot["display"].isna().any():
    missing = Dplot.loc[
        Dplot["display"].isna(),
        "analysis"
    ].tolist()
    raise RuntimeError(
        f"Unhandled sensitivity labels: {missing}"
    )

Dplot = Dplot.iloc[::-1].reset_index(drop=True)

# More vertical separation than v3.
y = np.arange(len(Dplot)) * 1.24

pearson_y = y + 0.16
spearman_y = y - 0.16

axD.scatter(
    Dplot["pearson_logFC"],
    pearson_y,
    s=120,
    marker="o",
    label="Pearson r",
    zorder=3,
)

axD.scatter(
    Dplot["spearman_logFC"],
    spearman_y,
    s=120,
    marker="s",
    label="Spearman ρ",
    zorder=3,
)

for i, row in Dplot.iterrows():
    lo = min(
        row["pearson_logFC"],
        row["spearman_logFC"],
    )
    hi = max(
        row["pearson_logFC"],
        row["spearman_logFC"],
    )

    axD.plot(
        [lo, hi],
        [y[i], y[i]],
        linewidth=1.6,
        color="black",
        alpha=0.28,
        zorder=1,
    )

axD.set_yticks(y)
axD.set_yticklabels(
    Dplot["display"],
    fontweight="bold",
    fontsize=16,
    linespacing=1.10,
)

xmin = min(
    Dplot["pearson_logFC"].min(),
    Dplot["spearman_logFC"].min(),
)

axD.set_xlim(
    max(0.95, xmin - 0.010),
    1.0022,
)

# Explicit vertical margins so annotations never approach axes.
axD.set_ylim(
    y.min() - 0.58,
    y.max() + 0.72,
)

axD.axvline(
    1.0,
    linestyle="--",
    linewidth=1.4,
    color="black",
    alpha=0.35,
)

axD.set_xlabel(
    "Effect-size correlation",
    labelpad=10,
)

axD.set_title(
    "Prespecified differential-expression sensitivities",
    loc="left",
    pad=16,
)

# Legend moved above/rightward and kept clear of data.
legend = axD.legend(
    frameon=False,
    loc="upper right",
    bbox_to_anchor=(0.985, 0.985),
    fontsize=15.5,
    handletextpad=0.55,
    borderaxespad=0.1,
)

for text in legend.get_texts():
    text.set_fontweight("bold")

axD.grid(
    axis="x",
    alpha=0.14,
    linewidth=1.0,
)

# Correlation values positioned consistently:
# - Pearson above-left of the circle
# - Spearman below-left of the square
# This avoids the crowding observed in v3.
for i, row in Dplot.iterrows():
    axD.annotate(
        f'{row["pearson_logFC"]:.3f}',
        (row["pearson_logFC"], pearson_y[i]),
        xytext=(-10, 11),
        textcoords="offset points",
        ha="right",
        va="bottom",
        fontsize=14.5,
        fontweight="bold",
        clip_on=False,
    )

    axD.annotate(
        f'{row["spearman_logFC"]:.3f}',
        (row["spearman_logFC"], spearman_y[i]),
        xytext=(-10, -11),
        textcoords="offset points",
        ha="right",
        va="top",
        fontsize=14.5,
        fontweight="bold",
        clip_on=False,
    )

style_axes(axD)


# ---------------------------------------------------------
# Panel labels
# ---------------------------------------------------------

fig.canvas.draw()

for ax, label in zip(
    [axA, axB, axC, axD],
    "ABCD",
):
    bbox = ax.get_position()
    fig.text(
        bbox.x0 - 0.034,
        bbox.y1 + 0.013,
        label,
        fontsize=25,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


# ---------------------------------------------------------
# Exports
# ---------------------------------------------------------

pdf = (
    OUT /
    "FigureS2_v4_differential_expression_robustness.pdf"
)

png = (
    OUT /
    "FigureS2_v4_differential_expression_robustness.png"
)

svg = (
    OUT /
    "FigureS2_v4_differential_expression_robustness.svg"
)

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


# ---------------------------------------------------------
# Console QC
# ---------------------------------------------------------

print("=== FIGURE S2 v4 GENERATION COMPLETE ===")
print(pdf)
print(svg)
print(png)

print("\n=== PANEL MERGE COUNTS ===")
print(f"Panel A genes: {len(A):,}")
print(f"Panel B genes: {len(B):,}")
print(f"Panel C genes: {len(C):,}")

print("\n=== FROZEN SENSITIVITY SUMMARY ===")
print(D.to_string(index=False))

print("\nNo differential-expression models were rerun.")
print(
    "Figure S2 v4 uses existing primary and "
    "prespecified sensitivity outputs only."
)
print(
    "Session 42 v4 is a presentation-only refinement "
    "of v3, restricted to Panel D layout/readability "
    "plus final export."
)
