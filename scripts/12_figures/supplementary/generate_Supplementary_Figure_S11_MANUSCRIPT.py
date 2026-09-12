#!/usr/bin/env python3
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path.home() / "Influenza_RSV_Project"
S38C = ROOT / "results/session38_additional_analysis/frozen/Session38C/tables"

PRIMARY = S38C / "Session38C_ARCHITECTURE_LOCALIZATION_v1.0.tsv"
ROBUST = S38C / "Session38C_ARCHITECTURE_DETECTION_ROBUSTNESS_v1.0.tsv"
PRIORITY = S38C / "Session38C_PRIORITY5_CELLTYPE_EFFECTS_v1.0.tsv"
BRIDGE = S38C / "Session38C_BRIDGE5_CELLTYPE_EFFECTS_v1.0.tsv"

for p in [PRIMARY, ROBUST, PRIORITY, BRIDGE]:
    if not p.exists():
        raise FileNotFoundError(p)

OUT = ROOT / "results/session40_script_organization/manuscript_supplementary_figures/regenerated/S11"
FIG = OUT / "figures"
PROV = OUT / "provenance"
for d in [FIG, PROV]:
    d.mkdir(parents=True, exist_ok=True)

PDF = FIG / "Supplementary_Figure_S11_single_cell_localization_robustness.pdf"
PNG = FIG / "Supplementary_Figure_S11_single_cell_localization_robustness.png"

primary = pd.read_csv(PRIMARY, sep="\t")
robust = pd.read_csv(ROBUST, sep="\t")
priority = pd.read_csv(PRIORITY, sep="\t")
bridge = pd.read_csv(BRIDGE, sep="\t")

assert len(primary) == 5
assert len(robust) == 5
assert len(priority) == 25
assert len(bridge) == 25

CELL_ORDER = [
    "natural killer cell",
    "CD8-positive, alpha-beta T cell",
    "effector CD8-positive, alpha-beta T cell",
    "effector CD4-positive, alpha-beta T cell",
    "classical monocyte",
]
CELL_SHORT = {
    "natural killer cell": "NK",
    "CD8-positive, alpha-beta T cell": "CD8 T",
    "effector CD8-positive, alpha-beta T cell": "Effector CD8 T",
    "effector CD4-positive, alpha-beta T cell": "Effector CD4 T",
    "classical monocyte": "Classical monocyte",
}
CELL_HEAT = {
    "natural killer cell": "NK",
    "CD8-positive, alpha-beta T cell": "CD8 T",
    "effector CD8-positive, alpha-beta T cell": "Eff. CD8 T",
    "effector CD4-positive, alpha-beta T cell": "Eff. CD4 T",
    "classical monocyte": "Classical\nmonocyte",
}

PRIORITY_ORDER = ["IFIH1","IFIT3","ISG15","STAT1","KPNB1"]
BRIDGE_ORDER = ["CHMP5","HERC5","KPNB1","OTOF","TOP2A"]

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

C_PRIMARY = "#2468B4"
C_ROBUST = "#F28C28"
C_DIR = "#6A3D9A"
C_NEG = "#D95F5F"
C_GRID = "#E6E6E6"
C_ZERO = "#666666"
C_TEXT = "#222222"

# Wider canvas, more column separation, larger bottom margin
fig = plt.figure(figsize=(21.5, 17.8))
gs = fig.add_gridspec(
    2, 2,
    left=0.10, right=0.965, top=0.94, bottom=0.13,
    wspace=0.50, hspace=0.46,
    width_ratios=[1.02, 1.05]
)
axA = fig.add_subplot(gs[0,0])
axB = fig.add_subplot(gs[0,1])
axC = fig.add_subplot(gs[1,0])
axD = fig.add_subplot(gs[1,1])

# ------------------------------------------------------------
# A: primary architecture localization
# ------------------------------------------------------------
p = primary.set_index("cell_type").loc[CELL_ORDER].reset_index()
y = np.arange(len(p))

axA.barh(
    y, p["spearman_rho"].astype(float),
    color=[C_PRIMARY if v >= 0 else C_NEG for v in p["spearman_rho"]],
    height=0.62, edgecolor="white", linewidth=1.0
)
axA.axvline(0, color=C_ZERO, lw=1.7)
axA.set_yticks(y)
axA.set_yticklabels([CELL_SHORT[x] for x in p["cell_type"]])
axA.invert_yaxis()
axA.set_xlabel("Spearman ρ: single-cell vs bulk effect")
axA.set_xlim(-0.28, 0.39)
axA.grid(axis="x", color=C_GRID, lw=1.0)
axA.set_axisbelow(True)
axA.spines[["top","right"]].set_visible(False)

# safer annotation placement: positives outside bar; negative value displayed inside/right of bar
for yi, (_, r) in enumerate(p.iterrows()):
    rho = float(r["spearman_rho"])
    fdr = float(r["spearman_FDR"])
    if rho >= 0:
        xtext = rho + 0.014
        ha = "left"
    else:
        # place annotation on the right side of the negative bar, away from y labels
        xtext = -0.145
        ha = "left"
    axA.text(
        xtext, yi,
        f"{rho:.3f}\nFDR={fdr:.3f}",
        ha=ha, va="center",
        fontsize=14.6, fontweight="bold"
    )

# ------------------------------------------------------------
# B: primary vs detection-filtered robustness
# ------------------------------------------------------------
r = robust.set_index("cell_type").loc[CELL_ORDER].reset_index()
y = np.arange(len(r))

for yi, (_, row) in enumerate(r.iterrows()):
    x1 = float(row["primary_spearman_rho"])
    x2 = float(row["robust_spearman_rho"])
    axB.plot([x1, x2], [yi, yi], color="#AFAFAF", lw=3.0, zorder=1)
    axB.scatter(x1, yi, s=115, color=C_PRIMARY, edgecolor="white", linewidth=1.2, zorder=3)
    axB.scatter(x2, yi, s=115, color=C_ROBUST, edgecolor="white", linewidth=1.2, zorder=3)

axB.axvline(0, color=C_ZERO, lw=1.7)
axB.set_yticks(y)
axB.set_yticklabels([CELL_SHORT[x] for x in r["cell_type"]])
axB.invert_yaxis()
axB.set_xlabel("Spearman ρ")
axB.set_xlim(-0.25, 0.41)
axB.grid(axis="x", color=C_GRID, lw=1.0)
axB.set_axisbelow(True)
axB.spines[["top","right"]].set_visible(False)

for yi, (_, row) in enumerate(r.iterrows()):
    axB.text(
        0.395, yi,
        f'{100*float(row["robust_direction_concordance"]):.1f}% dir.',
        ha="right", va="center",
        fontsize=14.2, color=C_DIR
    )

# legend moved above lower-left plotting area, avoiding data
axB.scatter([], [], s=115, color=C_PRIMARY, label="Primary")
axB.scatter([], [], s=115, color=C_ROBUST, label="Detection-filtered")
leg = axB.legend(
    frameon=False,
    loc="lower left",
    bbox_to_anchor=(0.01, 0.01),
    borderaxespad=0.0,
    handletextpad=0.5,
    labelspacing=0.45
)
for t in leg.get_texts():
    t.set_fontweight("bold")

# ------------------------------------------------------------
# heatmap helper
# ------------------------------------------------------------
def effect_heatmap(ax, df, genes, fdr_col):
    pivot = (
        df.pivot(index="gene_symbol", columns="cell_type", values="sc_effect_logCPM")
          .reindex(index=genes, columns=CELL_ORDER)
    )
    same = (
        df.pivot(index="gene_symbol", columns="cell_type", values="same_direction")
          .reindex(index=genes, columns=CELL_ORDER)
    )
    fdr = (
        df.pivot(index="gene_symbol", columns="cell_type", values=fdr_col)
          .reindex(index=genes, columns=CELL_ORDER)
    )

    vals = pivot.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(vals))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.imshow(vals, cmap="RdBu_r", norm=norm, aspect="auto")

    ax.set_xticks(np.arange(len(CELL_ORDER)))
    ax.set_xticklabels([CELL_HEAT[x] for x in CELL_ORDER], rotation=22, ha="right")
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes)
    ax.tick_params(axis="x", pad=8)

    for i in range(len(genes)):
        for j in range(len(CELL_ORDER)):
            val = vals[i, j]
            if np.isnan(val):
                continue
            is_same = bool(same.iloc[i, j])
            symbol = "✓" if is_same else "×"
            txt_color = "white" if abs(val) > 0.56*vmax else C_TEXT
            ax.text(
                j, i-0.08, f"{val:.2f}",
                ha="center", va="center",
                fontsize=14.0, color=txt_color, fontweight="bold"
            )
            ax.text(
                j, i+0.24, symbol,
                ha="center", va="center",
                fontsize=15.2, color=txt_color, fontweight="bold"
            )

    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_linewidth(1.5)

    return im

# ------------------------------------------------------------
# C / D heatmaps
# ------------------------------------------------------------
imC = effect_heatmap(axC, priority, PRIORITY_ORDER, "priority_FDR_25")
cbC = fig.colorbar(imC, ax=axC, fraction=0.038, pad=0.028)
cbC.set_label("Single-cell effect (logCPM)", fontsize=17.5, fontweight="bold", labelpad=9)
cbC.ax.tick_params(labelsize=14.0)

imD = effect_heatmap(axD, bridge, BRIDGE_ORDER, "bridge_FDR_25")
cbD = fig.colorbar(imD, ax=axD, fraction=0.038, pad=0.028)
cbD.set_label("Single-cell effect (logCPM)", fontsize=17.5, fontweight="bold", labelpad=9)
cbD.ax.tick_params(labelsize=14.0)

# explanatory note completely below both heatmaps
fig.text(
    0.53, 0.045,
    "Heatmaps: ✓ same direction as frozen bulk influenza effect; × opposite direction",
    ha="center", va="bottom",
    fontsize=15.0, fontweight="bold"
)

# ------------------------------------------------------------
# panel headers
# ------------------------------------------------------------
titles = {
    "A": "Localization of the frozen 170-gene architecture",
    "B": "Robustness to gene-detection filtering",
    "C": "Localization of weight-robust priority genes",
    "D": "Localization of response–dependency bridges",
}
axes = {"A":axA, "B":axB, "C":axC, "D":axD}
for lab, ax in axes.items():
    box = ax.get_position()
    fig.text(
        box.x0-0.028, box.y1+0.026, lab,
        ha="left", va="bottom",
        fontsize=26, fontweight="bold"
    )
    fig.text(
        box.x0+0.012, box.y1+0.026, titles[lab],
        ha="left", va="bottom",
        fontsize=21.8, fontweight="bold"
    )

fig.savefig(PDF, bbox_inches="tight")
fig.savefig(PNG, dpi=600, bbox_inches="tight")
plt.close(fig)

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

manifest = PROV / "Supplementary_Figure_S11_SHA256.txt"
with open(manifest, "w") as f:
    for p in [PRIMARY, ROBUST, PRIORITY, BRIDGE, PDF, PNG]:
        f.write(f"{sha256(p)}  {p}\n")

print("=== MANUSCRIPT SUPPLEMENTARY FIGURE S11 GENERATED ===")
print(PDF)
print(PNG)
print(manifest)
