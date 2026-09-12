#!/usr/bin/env python3
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path.home() / "Influenza_RSV_Project"
S38B = ROOT / "results/session38_additional_analysis/frozen/Session38B/tables"

OUT = ROOT / "results/session40_script_organization/manuscript_supplementary_figures/regenerated/S10"
FIG = OUT / "figures"
PROV = OUT / "provenance"
for d in [FIG, PROV]:
    d.mkdir(parents=True, exist_ok=True)

MASTER = S38B / "Session38B_MOLECULAR_FUNCTIONAL_CONVERGENCE_MASTER_v1.0.tsv"
STATS = S38B / "Session38B_STATISTICAL_TESTS_v1.0.tsv"
BRIDGE = S38B / "Session38B_CRISPR_BRIDGE_MOLECULAR_CONTEXT_v1.0.tsv"
CLASSES = S38B / "Session38B_CONVERGENCE_CLASS_SUMMARY_v1.0.tsv"

for p in [MASTER, STATS, BRIDGE, CLASSES]:
    if not p.exists():
        raise FileNotFoundError(p)

PDF = FIG / "Supplementary_Figure_S10_molecular_functional_convergence_robustness.pdf"
PNG = FIG / "Supplementary_Figure_S10_molecular_functional_convergence_robustness.png"

master = pd.read_csv(MASTER, sep="\t")
stats = pd.read_csv(STATS, sep="\t")
bridge = pd.read_csv(BRIDGE, sep="\t")
classes = pd.read_csv(CLASSES, sep="\t")

assert len(master) == 170 and master["gene_symbol"].nunique() == 170
assert len(stats) == 9
assert len(bridge) == 7
assert int(classes["n_genes"].sum()) == 170

def stat(analysis, test=None):
    x = stats.loc[stats["analysis"].eq(analysis)]
    if test is not None:
        x = x.loc[x["test"].eq(test)]
    if len(x) != 1:
        raise RuntimeError(f"Expected one row for analysis={analysis}, test={test}; found {len(x)}")
    return x.iloc[0]

disc_pear = stat("Discovery_RNA_vs_protein", "Pearson")
disc_spear = stat("Discovery_RNA_vs_protein", "Spearman")
disc_sign = stat("Discovery_RNA_vs_protein_sign", "Exact_binomial_vs_0.5")

ind_pear = stat("Independent_RNA_vs_protein", "Pearson")
ind_spear = stat("Independent_RNA_vs_protein", "Spearman")
ind_sign = stat("Independent_RNA_vs_protein_sign", "Exact_binomial_vs_0.5")

crispr_rows = stats.loc[stats["analysis"].isin([
    "Abs_discovery_RNA_vs_CRISPR_dependency",
    "Abs_independent_RNA_vs_CRISPR_dependency",
    "Abs_protein_vs_CRISPR_dependency",
])].copy()
assert len(crispr_rows) == 3

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 18,
    "font.weight": "bold",
    "axes.labelsize": 21,
    "axes.labelweight": "bold",
    "xtick.labelsize": 16.5,
    "ytick.labelsize": 16.5,
    "legend.fontsize": 14.5,
    "axes.linewidth": 2.4,
    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

C_RNA = "#2A78C6"
C_PROT = "#F28C28"
C_CRISPR = "#D94B4B"
C_SHARED = "#76A9CF"
C_FLU = "#F47A3C"
C_HC = "#7A3FA1"
C_GRID = "#E5E5E5"
C_ZERO = "#666666"
C_TEXT = "#222222"
C_NA = "#E4E4E4"

# Wider canvas + extra lower-row separation
fig = plt.figure(figsize=(21.2, 17.6))
gs = fig.add_gridspec(
    2, 2,
    left=0.09, right=0.965, top=0.94, bottom=0.105,
    wspace=0.48, hspace=0.46,
    width_ratios=[1.0, 1.07]
)
axA = fig.add_subplot(gs[0,0])
axB = fig.add_subplot(gs[0,1])
axC = fig.add_subplot(gs[1,0])
axD = fig.add_subplot(gs[1,1])

# ---------------- A ----------------
a = master.loc[
    master["proteomics_evaluable"].fillna(False)
    & master["discovery_iav_logFC"].notna()
    & master["proteomics_median_logFC"].notna()
].copy()

arch_col = "proteomics_architecture_class" if "proteomics_architecture_class" in a.columns else None
if arch_col:
    colorsA = [C_FLU if str(v) == "influenza_amplified" else C_SHARED for v in a[arch_col]]
else:
    colorsA = C_SHARED

axA.scatter(
    a["discovery_iav_logFC"],
    a["proteomics_median_logFC"],
    s=72, c=colorsA, alpha=0.86,
    edgecolor="white", linewidth=0.7
)
axA.axhline(0, color=C_ZERO, lw=1.5, ls="--")
axA.axvline(0, color=C_ZERO, lw=1.5, ls="--")
axA.set_xlabel("Discovery influenza RNA logFC")
axA.set_ylabel("Protein median logFC")
axA.grid(color=C_GRID, lw=1.0)
axA.set_axisbelow(True)
axA.spines[["top","right"]].set_visible(False)

axA.text(
    0.03, 0.97,
    f'Pearson r = {disc_pear["effect"]:.3f}, P = {disc_pear["p_value"]:.4f}\n'
    f'Spearman ρ = {disc_spear["effect"]:.3f}, P = {disc_spear["p_value"]:.4f}\n'
    f'Same sign = {100*disc_sign["effect"]:.1f}% ({int(disc_sign["n"])} genes), '
    f'P = {disc_sign["p_value"]:.4f}',
    transform=axA.transAxes,
    ha="left", va="top",
    fontsize=16.2, fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.42", facecolor="white",
              edgecolor="#CCCCCC", alpha=0.95)
)

# ---------------- B ----------------
b = master.loc[
    master["proteomics_evaluable"].fillna(False)
    & master["rnaseq_mean_logFC"].notna()
    & master["proteomics_median_logFC"].notna()
].copy()

if arch_col:
    colorsB = [C_FLU if str(v) == "influenza_amplified" else C_SHARED for v in b[arch_col]]
else:
    colorsB = C_SHARED

axB.scatter(
    b["rnaseq_mean_logFC"],
    b["proteomics_median_logFC"],
    s=72, c=colorsB, alpha=0.86,
    edgecolor="white", linewidth=0.7
)
axB.axhline(0, color=C_ZERO, lw=1.5, ls="--")
axB.axvline(0, color=C_ZERO, lw=1.5, ls="--")
axB.set_xlabel("Independent RNA-seq mean logFC")
axB.set_ylabel("Protein median logFC")
axB.grid(color=C_GRID, lw=1.0)
axB.set_axisbelow(True)
axB.spines[["top","right"]].set_visible(False)

axB.text(
    0.03, 0.97,
    f'Pearson r = {ind_pear["effect"]:.3f}, P = {ind_pear["p_value"]:.4f}\n'
    f'Spearman ρ = {ind_spear["effect"]:.3f}, P = {ind_spear["p_value"]:.4f}\n'
    f'Same sign = {100*ind_sign["effect"]:.1f}% ({int(ind_sign["n"])} genes), '
    f'P = {ind_sign["p_value"]:.4f}',
    transform=axB.transAxes,
    ha="left", va="top",
    fontsize=16.2, fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.42", facecolor="white",
              edgecolor="#CCCCCC", alpha=0.95)
)

handles = [
    Line2D([0],[0], marker="o", linestyle="None", markersize=9,
           markerfacecolor=C_SHARED, markeredgecolor="white", label="Shared core / other"),
    Line2D([0],[0], marker="o", linestyle="None", markersize=9,
           markerfacecolor=C_FLU, markeredgecolor="white", label="Influenza-amplified"),
]
leg = axB.legend(handles=handles, frameon=False, loc="lower right")
for t in leg.get_texts():
    t.set_fontweight("bold")

# ---------------- C ----------------
order = [
    "Abs_discovery_RNA_vs_CRISPR_dependency",
    "Abs_independent_RNA_vs_CRISPR_dependency",
    "Abs_protein_vs_CRISPR_dependency",
]
labels = [
    "Discovery RNA\nmagnitude",
    "Independent RNA\nmagnitude",
    "Protein\nmagnitude",
]
cc = crispr_rows.set_index("analysis").loc[order].reset_index()
x = np.arange(3)
effects = cc["effect"].astype(float).to_numpy()

axC.bar(
    x, effects,
    color=[C_RNA, "#55A868", C_PROT],
    width=0.62,
    edgecolor="white", linewidth=1.0
)
axC.axhline(0, color=C_ZERO, lw=1.8)
axC.set_xticks(x)
axC.set_xticklabels(labels)
axC.set_ylabel("Spearman ρ with CRISPR\ndependency strength", labelpad=14)
axC.set_ylim(-0.16, 0.075)
axC.grid(axis="y", color=C_GRID, lw=1.0)
axC.set_axisbelow(True)
axC.spines[["top","right"]].set_visible(False)

for xi, (_, r) in enumerate(cc.iterrows()):
    y = float(r["effect"])
    fdr = float(r["FDR"])
    p = float(r["p_value"])
    n = int(r["n"])
    va = "top" if y < 0 else "bottom"
    offset = -0.010 if y < 0 else 0.010
    axC.text(
        xi, y + offset,
        f'ρ={y:.3f}\nP={p:.3f}\nFDR={fdr:.3f}\nn={n}',
        ha="center", va=va,
        fontsize=15.0, fontweight="bold"
    )

# Compact annotation, moved inside upper-left plot area
axC.text(
    0.03, 0.95,
    "No test survived FDR correction",
    transform=axC.transAxes,
    ha="left", va="top",
    fontsize=15.0, fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.34", facecolor="white",
              edgecolor="#CCCCCC", alpha=0.95)
)

# ---------------- D ----------------
bridge_order = ["KPNB1","CHMP5","TOP2A","FCGR1B","HIST2H2AC","OTOF","HERC5"]
d = bridge.set_index("gene_symbol").loc[bridge_order].reset_index()
y = np.arange(len(d))
cols = ["Discovery\nRNA", "Independent\nRNA", "Protein", "CRISPR"]
xx = np.arange(4)

for yi, (_, r) in enumerate(d.iterrows()):
    axD.plot([0,3], [yi, yi], color="#EFEFEF", lw=1.2, zorder=0)

    axD.scatter(0, yi, s=150, facecolor=C_RNA, edgecolor="white", linewidth=1.0, zorder=3)
    axD.text(0, yi+0.28, f'{r["discovery_iav_logFC"]:.2f}',
             ha="center", va="top", fontsize=13.3, color=C_TEXT)

    axD.scatter(1, yi, s=150, facecolor="#55A868", edgecolor="white", linewidth=1.0, zorder=3)
    axD.text(1, yi+0.28, f'{r["rnaseq_mean_logFC"]:.2f}',
             ha="center", va="top", fontsize=13.3, color=C_TEXT)

    if not bool(r["proteomics_evaluable"]):
        axD.scatter(2, yi, s=150, facecolor=C_NA, edgecolor="#AAAAAA", linewidth=1.4, zorder=3)
        axD.text(2, yi+0.28, "NA", ha="center", va="top",
                 fontsize=13.3, color="#777777")
    else:
        cat = str(r["proteomics_validation_category"])
        concordant = cat.startswith("PROTEIN_CONCORDANT")
        fc = C_PROT if concordant else "white"
        axD.scatter(2, yi, s=150, facecolor=fc, edgecolor=C_PROT, linewidth=2.0, zorder=3)
        axD.text(2, yi+0.28, f'{r["proteomics_median_logFC"]:.2f}',
                 ha="center", va="top", fontsize=13.3, color=C_TEXT)

    axD.scatter(3, yi, s=150, facecolor=C_CRISPR, edgecolor="white", linewidth=1.0, zorder=3)
    if bool(r["crispr_iav_dependency_high_confidence"]):
        axD.scatter(3, yi, s=220, facecolor="none", edgecolor=C_HC,
                    linewidth=2.6, zorder=4)
    axD.text(3, yi+0.28, f'{r["crispr_dependency_strength"]:.3f}',
             ha="center", va="top", fontsize=13.3, color=C_TEXT)

axD.set_yticks(y)
axD.set_yticklabels(d["gene_symbol"])
axD.set_xticks(xx)
axD.set_xticklabels(cols)
axD.set_xlim(-0.6, 3.6)
axD.set_ylim(len(d)-0.40, -0.65)
axD.tick_params(axis="both", length=0)
axD.tick_params(axis="x", pad=9)
axD.spines[["top","right","bottom","left"]].set_visible(False)

legendD = [
    Line2D([0],[0], marker="o", linestyle="None", markersize=9,
           markerfacecolor=C_PROT, markeredgecolor=C_PROT,
           label="Protein concordant"),
    Line2D([0],[0], marker="o", linestyle="None", markersize=9,
           markerfacecolor="white", markeredgecolor=C_PROT,
           markeredgewidth=2, label="Protein discordant"),
    Line2D([0],[0], marker="o", linestyle="None", markersize=9,
           markerfacecolor=C_NA, markeredgecolor="#AAAAAA",
           label="Protein not evaluable"),
    Line2D([0],[0], marker="o", linestyle="None", markersize=11,
           markerfacecolor="none", markeredgecolor=C_HC,
           markeredgewidth=2, label="High-confidence CRISPR"),
]
legD = axD.legend(
    handles=legendD,
    frameon=False,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.18),
    ncol=2,
    fontsize=13.8,
    columnspacing=1.6,
    handletextpad=0.6,
    borderaxespad=0.0
)
for t in legD.get_texts():
    t.set_fontweight("bold")

# ---------------- panel headers ----------------
titles = {
    "A": "Discovery RNA–protein concordance",
    "B": "Independent RNA–protein concordance",
    "C": "Molecular magnitude versus CRISPR dependency",
    "D": "Molecular context of response–dependency bridges",
}
axes = {"A": axA, "B": axB, "C": axC, "D": axD}

for lab, ax in axes.items():
    box = ax.get_position()
    fig.text(box.x0-0.028, box.y1+0.026, lab,
             ha="left", va="bottom", fontsize=26, fontweight="bold")
    fig.text(box.x0+0.012, box.y1+0.026, titles[lab],
             ha="left", va="bottom", fontsize=22.0, fontweight="bold")

fig.savefig(PDF, bbox_inches="tight")
fig.savefig(PNG, dpi=600, bbox_inches="tight")
plt.close(fig)

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

manifest = PROV / "Supplementary_Figure_S10_SHA256.txt"
with open(manifest, "w") as f:
    for p in [MASTER, STATS, BRIDGE, CLASSES, PDF, PNG]:
        f.write(f"{sha256(p)}  {p}\n")

print("=== MANUSCRIPT SUPPLEMENTARY FIGURE S10 GENERATED ===")
print(PDF)
print(PNG)
print(manifest)
