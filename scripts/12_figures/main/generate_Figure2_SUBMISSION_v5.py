#!/usr/bin/env python3

from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT = Path.home() / "Influenza_RSV_Project"

DE_FILE = (
    ROOT / "results/DS001_GSE38900/differential_expression/tables/"
    "GPL6884_InfluenzaA_vs_control_full_DE.tsv"
)

CLASS_FILE = (
    ROOT / "results/DS001_GSE38900/signature_analysis/tables/"
    "DS001_GSE38900_gene_evidence_master_final_classification.tsv"
)

SUMMARY_FILE = (
    ROOT / "results/DS001_GSE38900/integrated_architecture/tables/"
    "Figure7A_frozen_signature_summary.tsv"
)

OUTDIR = ROOT / "results/figure2_influenza_centered"
FIGDIR = OUTDIR / "figures"
TABDIR = OUTDIR / "tables"

FIGDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)

de = pd.read_csv(DE_FILE, sep="\t")
cl = pd.read_csv(CLASS_FILE, sep="\t")
summary = pd.read_csv(SUMMARY_FILE, sep="\t")

# ============================================================
# QC
# ============================================================

sig = de["adj.P.Val"] < 0.05
n_fdr = int(sig.sum())
n_fc05 = int((sig & (de["logFC"].abs() >= 0.5)).sum())
n_fc10 = int((sig & (de["logFC"].abs() >= 1.0)).sum())

boolcol = lambda x: cl[x].fillna(False).astype(bool)

n_shared = int(boolcol("shared_replicated").sum())
n_core = int(boolcol("shared_replicated_core").sum())
n_fluamp = int(boolcol("shared_replicated_influenza_amplified").sum())
n_rsvamp = int(boolcol("shared_replicated_RSV_amplified").sum())

assert n_fdr == 3501
assert n_fc05 == 1934
assert n_fc10 == 323
assert n_shared == 170
assert n_core == 130
assert n_fluamp == 37
assert n_rsvamp == 3
assert n_core + n_fluamp + n_rsvamp == n_shared

print("=== QC COUNTS ===")
print("FDR < 0.05:", n_fdr)
print("FDR + |logFC| >= 0.5:", n_fc05)
print("FDR + |logFC| >= 1.0:", n_fc10)
print("Shared replicated:", n_shared)
print("Core:", n_core)
print("Influenza-amplified:", n_fluamp)
print("RSV-amplified:", n_rsvamp)

# ============================================================
# BUILD 170-GENE TABLE
# ============================================================

rep = cl.loc[
    boolcol("shared_replicated"),
    [
        "gene_symbol",
        "flu_logFC",
        "flu_FDR",
        "rsv6884_logFC",
        "rsv6884_FDR",
        "shared_replicated_core",
        "shared_replicated_influenza_amplified",
        "shared_replicated_RSV_amplified",
    ],
].copy()

def architecture_class(row):
    if bool(row["shared_replicated_influenza_amplified"]):
        return "Influenza-amplified"
    if bool(row["shared_replicated_RSV_amplified"]):
        return "RSV-amplified"
    if bool(row["shared_replicated_core"]):
        return "Shared core"
    return "Unclassified"

rep["architecture_class"] = rep.apply(architecture_class, axis=1)

assert len(rep) == 170
assert (rep["architecture_class"] != "Unclassified").all()

rep.to_csv(
    TABDIR / "Figure2C_v8_replicated_170_gene_architecture.tsv",
    sep="\t",
    index=False
)

# ============================================================
# STYLE
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 13.0,
    "axes.titlesize": 15.0,
    "axes.labelsize": 13.8,
    "axes.labelweight": "bold",
    "xtick.labelsize": 11.8,
    "ytick.labelsize": 11.8,
    "legend.fontsize": 12.0,
    "axes.linewidth": 1.70,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

C_GRAY = "#B8BDC5"
C_DARK = "#333333"
C_FLU = "#C44E52"
C_CORE = "#4C72B0"
C_RSV = "#55A868"
C_LIGHT = "#E7E8EA"

# ============================================================
# LAYOUT: A-B / C-D
# ============================================================

fig = plt.figure(figsize=(18.6, 12.6))

gs = GridSpec(
    2, 2,
    figure=fig,
    width_ratios=[1.0, 1.22],
    height_ratios=[1.0, 1.06],
    wspace=0.30,
    hspace=0.50
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

# ============================================================
# PANEL A
# ============================================================

labels = [
    "FDR < 0.05",
    "FDR < 0.05\n|log$_2$FC| ≥ 0.5",
    "FDR < 0.05\n|log$_2$FC| ≥ 1.0",
]

counts = [n_fdr, n_fc05, n_fc10]
y = np.arange(3)

axA.barh(
    y,
    counts,
    color=[C_LIGHT, "#D7A3A5", C_FLU],
    height=0.58,
    edgecolor="none"
)

for yi, val in zip(y, counts):
    axA.text(
        val + 70,
        yi,
        f"{val:,}",
        va="center",
        fontsize=13.5,
        fontweight="bold"
    )

axA.set_yticks(y)
axA.set_yticklabels(labels)
axA.invert_yaxis()
axA.set_xlim(0, 3900)
axA.set_xlabel("Number of genes")

axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)
axA.spines["left"].set_visible(False)
axA.tick_params(axis="y", length=0)
axA.grid(axis="x", alpha=0.30, linewidth=1.0)

# ============================================================
# PANEL B
# ============================================================

plot_de = de.copy()

tiny = np.nextafter(0, 1)
plot_de["minuslog10FDR"] = -np.log10(
    plot_de["adj.P.Val"].clip(lower=tiny)
)

significant = plot_de["adj.P.Val"] < 0.05
strong = significant & (plot_de["logFC"].abs() >= 0.5)

axB.scatter(
    plot_de.loc[~significant, "logFC"],
    plot_de.loc[~significant, "minuslog10FDR"],
    s=6.5,
    alpha=0.18,
    color=C_GRAY,
    linewidths=0,
    rasterized=True
)

axB.scatter(
    plot_de.loc[significant & ~strong, "logFC"],
    plot_de.loc[significant & ~strong, "minuslog10FDR"],
    s=8,
    alpha=0.27,
    color="#929292",
    linewidths=0,
    rasterized=True
)

axB.scatter(
    plot_de.loc[strong, "logFC"],
    plot_de.loc[strong, "minuslog10FDR"],
    s=11,
    alpha=0.52,
    color=C_FLU,
    linewidths=0,
    rasterized=True
)

axB.axvline(0, color=C_DARK, lw=1.55)
axB.axvline(0.5, color=C_DARK, lw=1.30, ls="--", alpha=0.80)
axB.axvline(-0.5, color=C_DARK, lw=1.30, ls="--", alpha=0.80)
axB.axhline(
    -np.log10(0.05),
    color=C_DARK,
    lw=1.30,
    ls="--",
    alpha=0.80
)

# Reduced, deliberately positioned annotation set
label_offsets_B = {
    "IFI27": (8, 5),
    "OAS3": (-10, 8),
    "ISG15": (8, -1),
    "IRF7": (8, 8),
    "IFIH1": (-8, -10),
    "STAT1": (8, -10),
}

for gene, offset in label_offsets_B.items():
    sub = (
        plot_de.loc[plot_de["gene_symbol"] == gene]
        .sort_values("adj.P.Val")
    )

    if len(sub) == 0:
        continue

    row = sub.iloc[0]

    axB.annotate(
        gene,
        (row["logFC"], row["minuslog10FDR"]),
        xytext=offset,
        textcoords="offset points",
        fontsize=12.2,
        fontweight="bold",
        ha="left" if offset[0] >= 0 else "right",
        va="bottom",
        arrowprops=dict(
            arrowstyle="-",
            lw=1.05,
            color="#666666"
        )
    )

axB.set_xlabel("Influenza A vs control log$_2$ fold change")
axB.set_ylabel("−log$_{10}$ FDR")

axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

# ============================================================
# PANEL C
# ============================================================

class_style = {
    "Shared core": (C_CORE, 28, 0.67),
    "Influenza-amplified": (C_FLU, 45, 0.90),
    "RSV-amplified": (C_RSV, 52, 0.97),
}

for cls_name in [
    "Shared core",
    "Influenza-amplified",
    "RSV-amplified",
]:
    sub = rep[rep["architecture_class"] == cls_name]
    color, size, alpha = class_style[cls_name]

    axC.scatter(
        sub["rsv6884_logFC"],
        sub["flu_logFC"],
        s=size,
        color=color,
        alpha=alpha,
        edgecolor="white",
        linewidth=0.35,
        label=f"{cls_name} (n={len(sub)})",
        zorder=3
    )

vals = np.concatenate([
    rep["rsv6884_logFC"].dropna().values,
    rep["flu_logFC"].dropna().values
])

lo = min(vals.min(), 0) - 0.25
hi = max(vals.max(), 0) + 0.25

axC.plot(
    [lo, hi],
    [lo, hi],
    ls="--",
    lw=1.60,
    color=C_DARK,
    alpha=0.85
)

axC.axhline(0, color=C_DARK, lw=1.30, alpha=0.85)
axC.axvline(0, color=C_DARK, lw=1.30, alpha=0.85)

axC.set_xlim(lo, hi)
axC.set_ylim(lo, hi)

label_offsets_C = {
    "OAS3": (7, 10),
    "ISG15": (10, 2),
    "IRF7": (8, 8),
    "IFIH1": (8, -10),
}

for gene, offset in label_offsets_C.items():
    sub = rep[
        (rep["gene_symbol"] == gene) &
        (rep["architecture_class"] == "Influenza-amplified")
    ]

    if len(sub) == 0:
        continue

    row = sub.iloc[0]

    axC.annotate(
        gene,
        (row["rsv6884_logFC"], row["flu_logFC"]),
        xytext=offset,
        textcoords="offset points",
        fontsize=12.2,
        fontweight="bold",
        ha="left" if offset[0] >= 0 else "right",
        arrowprops=dict(
            arrowstyle="-",
            lw=1.0,
            color="#666666"
        )
    )

axC.set_xlabel("RSV vs control log$_2$ fold change")
axC.set_ylabel("Influenza A vs control log$_2$ fold change")


axC.legend(
    frameon=False,
    loc="lower right",
    borderaxespad=0.5,
    handletextpad=0.6,
    labelspacing=0.55
)

axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)

# ============================================================
# PANEL D
# ============================================================

axD.axis("off")
axD.set_xlim(0, 1)
axD.set_ylim(0, 1)


# Source evidence
axD.text(
    0.50,
    0.865,
    f"Influenza A response\n{n_fdr:,} FDR-significant genes",
    ha="center",
    va="center",
    fontsize=14.4,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.72",
        facecolor="#F5F5F5",
        edgecolor="#777777",
        linewidth=1.60
    )
)

axD.annotate(
    "",
    xy=(0.50, 0.685),
    xytext=(0.50, 0.775),
    arrowprops=dict(
        arrowstyle="-|>",
        color="#777777",
        lw=2.00
    )
)

axD.text(
    0.50,
    0.615,
    f"Replicated shared antiviral backbone\nn = {n_shared}",
    ha="center",
    va="center",
    fontsize=14.8,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.72",
        facecolor="#EEF2F7",
        edgecolor=C_CORE,
        linewidth=2.20
    )
)

positions = [0.13, 0.50, 0.87]

boxes = [
    ("Shared core", n_core, C_CORE),
    ("Influenza-\namplified", n_fluamp, C_FLU),
    ("RSV-\namplified", n_rsvamp, C_RSV),
]

for xpos, (name, number, color) in zip(positions, boxes):

    axD.annotate(
        "",
        xy=(xpos, 0.305),
        xytext=(0.50, 0.515),
        arrowprops=dict(
            arrowstyle="-|>",
            color="#999999",
            lw=1.95
        )
    )

    axD.text(
        xpos,
        0.165,
        f"{name}\n{number} genes",
        ha="center",
        va="center",
        fontsize=14.4,
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.68",
            facecolor="white",
            edgecolor=color,
            linewidth=2.20
        )
    )

# ============================================================
# PANEL HEADINGS — LOCKED ALIGNMENT + BOLD TITLES
# ============================================================

PANEL_LABEL_SIZE = 20.0
PANEL_TITLE_SIZE = 16.0
PANEL_HEADING_Y = 1.075
PANEL_LABEL_X = -0.095
PANEL_TITLE_X = 0.00

def add_panel_heading(ax, letter, title):
    # Letter and title share exactly the same baseline/vertical anchor.
    ax.text(
        PANEL_LABEL_X, PANEL_HEADING_Y, letter,
        transform=ax.transAxes,
        fontsize=PANEL_LABEL_SIZE, fontweight="bold",
        ha="left", va="bottom", clip_on=False
    )
    ax.text(
        PANEL_TITLE_X, PANEL_HEADING_Y, title,
        transform=ax.transAxes,
        fontsize=PANEL_TITLE_SIZE, fontweight="bold",
        ha="left", va="bottom", clip_on=False,
        linespacing=1.05
    )

add_panel_heading(axA, "A", "Influenza A transcriptional response")
add_panel_heading(axB, "B", "Differential expression in acute influenza A")
add_panel_heading(axC, "C", "Influenza amplification within a\nreplicated antiviral backbone")
add_panel_heading(axD, "D", "Influenza-centered response architecture")

# No figure-level title in manuscript-facing version.


fig.subplots_adjust(top=0.94, bottom=0.075, left=0.075, right=0.985)


# ============================================================
# TYPOGRAPHY LOCK — AXES / TICKS / LEGENDS
# ============================================================

for ax in [axA, axB, axC]:
    ax.xaxis.label.set_fontweight("bold")
    ax.yaxis.label.set_fontweight("bold")
    ax.tick_params(axis="both", width=1.45, length=5.5)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")

# Panel A categorical labels are manuscript text, not numerical ticks.
for tick in axA.get_yticklabels():
    tick.set_fontweight("bold")

leg = axC.get_legend()
if leg is not None:
    for txt in leg.get_texts():
        txt.set_fontweight("bold")

# ============================================================
# SAVE
# ============================================================

pdf = FIGDIR / "Figure2_SUBMISSION_v5_influenza_centered_gene_architecture.pdf"
png = FIGDIR / "Figure2_SUBMISSION_v5_influenza_centered_gene_architecture.png"

fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=400, bbox_inches="tight")
plt.close(fig)

# ============================================================
# PROVENANCE
# ============================================================

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

prov = pd.DataFrame([
    {
        "file": str(DE_FILE.relative_to(ROOT)),
        "sha256": sha256(DE_FILE),
        "role": "Influenza differential expression"
    },
    {
        "file": str(CLASS_FILE.relative_to(ROOT)),
        "sha256": sha256(CLASS_FILE),
        "role": "Frozen gene classification"
    },
    {
        "file": str(SUMMARY_FILE.relative_to(ROOT)),
        "sha256": sha256(SUMMARY_FILE),
        "role": "Frozen signature summary"
    },
    {
        "file": str(pdf.relative_to(ROOT)),
        "sha256": sha256(pdf),
        "role": "Figure 2 submission PDF"
    },
    {
        "file": str(png.relative_to(ROOT)),
        "sha256": sha256(png),
        "role": "Figure 2 submission PNG"
    },
])

prov.to_csv(
    TABDIR / "Figure2_SUBMISSION_v5_provenance_sha256.tsv",
    sep="\t",
    index=False
)

print("\n=== WRITTEN ===")
print(pdf)
print(png)
print(TABDIR / "Figure2C_v8_replicated_170_gene_architecture.tsv")
print(TABDIR / "Figure2_SUBMISSION_v5_provenance_sha256.tsv")
print("\nSTATUS: FIGURE 2 SUBMISSION v5 COMPLETE")
