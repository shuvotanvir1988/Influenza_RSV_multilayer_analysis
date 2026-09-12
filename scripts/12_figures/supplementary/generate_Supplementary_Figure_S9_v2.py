#!/usr/bin/env python3
from pathlib import Path
import hashlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path.home() / "Influenza_RSV_Project"
S38 = ROOT / "results/session38_additional_analysis"
FROZEN = S38 / "frozen" / "Session38A"

OUT = ROOT / "results/session42_supplementary_figure_editing/FigureS9_evidence_weighted_priority_robustness/v2"
FIG = OUT / "figures"
PROV = OUT / "provenance"
for d in [FIG, PROV]:
    d.mkdir(parents=True, exist_ok=True)

PDF = FIG / "Supplementary_Figure_S9_v2_evidence_weighted_priority_robustness.pdf"
PNG = FIG / "Supplementary_Figure_S9_v2_evidence_weighted_priority_robustness.png"
SVG = FIG / "Supplementary_Figure_S9_v2_evidence_weighted_priority_robustness.svg"

def resolve(name):
    candidates = [
        FROZEN / "tables" / name,
        S38 / "tables" / name,
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"Could not locate {name}. Checked:\n" +
        "\n".join(str(x) for x in candidates)
    )

MASTER = resolve("Session38A_EVIDENCE_WEIGHTED_GENE_PRIORITY_MASTER_v1.0.tsv")
TOP30 = resolve("Session38A_TOP30_EVIDENCE_WEIGHTED_GENES_v1.0.tsv")
RANKS = resolve("Session38A_WEIGHT_SENSITIVITY_GENE_RANKS_v1.0.tsv")
CORRS = resolve("Session38A_WEIGHT_SENSITIVITY_CORRELATIONS_v1.0.tsv")

master = pd.read_csv(MASTER, sep="\t")
top30 = pd.read_csv(TOP30, sep="\t")
ranks = pd.read_csv(RANKS, sep="\t")
corrs = pd.read_csv(CORRS, sep="\t")

assert len(master) == 170 and master["gene_symbol"].nunique() == 170
assert len(top30) == 30
assert ranks["gene_symbol"].nunique() == 170
assert len(corrs) == 6

priority5 = ["KPNB1", "ISG15", "STAT1", "IFIT3", "IFIH1"]
bridge7 = {"KPNB1","HERC5","CHMP5","TOP2A","OTOF","FCGR1B","HIST2H2AC"}

# -------------------- style --------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 18,
    "font.weight": "bold",
    "axes.labelsize": 22,
    "axes.labelweight": "bold",
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "legend.fontsize": 15,
    "axes.linewidth": 2.4,
    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

C_SHARED = "#69A7D3"
C_FLU = "#F26B38"
C_PRIORITY = "#6B3FA0"
C_BRIDGE = "#E53935"
C_PRIMARY = "#1565C0"
C_IND = "#2EAD4A"
C_FUNC = "#F28C18"
C_RESP = "#9C4DCC"
C_GRID = "#E7E7E7"
C_TEXT = "#222222"

# Wider canvas and stronger inter-panel separation than v1
fig = plt.figure(figsize=(25.0, 20.5))
gs = fig.add_gridspec(
    2, 2,
    left=0.075, right=0.965, top=0.94, bottom=0.080,
    wspace=0.68, hspace=0.56,
    width_ratios=[1.08, 1.12]
)
axA = fig.add_subplot(gs[0,0])
axB = fig.add_subplot(gs[0,1])
axC = fig.add_subplot(gs[1,0])
axD = fig.add_subplot(gs[1,1])

# -------------------- A: complete top-30 ranking --------------------
p = top30.sort_values("session38_rank", ascending=False).copy()
y = np.arange(len(p))
bar_colors = [
    C_FLU if bool(v) else C_SHARED
    for v in p.get("influenza_amplified", pd.Series(False, index=p.index)).fillna(False)
]
axA.barh(
    y,
    p["evidence_score_coverage_adjusted"].astype(float),
    color=bar_colors,
    edgecolor="white",
    linewidth=0.8,
    height=0.72
)
axA.set_yticks(y)
axA.set_yticklabels(p["gene_symbol"], fontsize=18, fontweight="bold")
axA.set_xlabel("Coverage-adjusted evidence score (%)")
axA.grid(axis="x", color=C_GRID, lw=1.15)
axA.set_axisbelow(True)
axA.spines[["top","right"]].set_visible(False)

for yi, (_, r) in enumerate(p.iterrows()):
    g = r["gene_symbol"]
    x = float(r["evidence_score_coverage_adjusted"])
    if g in priority5:
        axA.scatter(x-1.6, yi, s=95, marker="o",
                    facecolor="white", edgecolor=C_PRIORITY,
                    linewidth=2.2, zorder=5)
    if g in bridge7:
        axA.scatter(x-4.2, yi, s=100, marker="D",
                    facecolor=C_BRIDGE, edgecolor="white",
                    linewidth=1.0, zorder=5)

legend = [
    Line2D([0],[0], color=C_FLU, lw=9, label="Influenza-amplified"),
    Line2D([0],[0], color=C_SHARED, lw=9, label="Shared core / other"),
    Line2D([0],[0], marker="o", linestyle="None", markersize=8,
           markerfacecolor="white", markeredgecolor=C_PRIORITY,
           markeredgewidth=2, label="Weight-robust priority-5"),
    Line2D([0],[0], marker="D", linestyle="None", markersize=8,
           markerfacecolor=C_BRIDGE, markeredgecolor="white",
           label="Response–dependency bridge"),
]
# Move legend to unused upper-left plotting space and compact it
legA = axA.legend(
    handles=legend,
    frameon=False,
    fontsize=16.0,
    loc="lower right",
    bbox_to_anchor=(0.995, 0.015),
    borderaxespad=0.0,
    handlelength=1.8,
    handletextpad=0.7,
    labelspacing=0.45
)
for t in legA.get_texts():
    t.set_fontweight("bold")

# -------------------- B: rank trajectories across weighting models --------------------
models = [
    ("primary_rank", "Primary", C_PRIMARY),
    ("independent_heavy_rank", "Independent-heavy", C_IND),
    ("functional_heavy_rank", "Functional-heavy", C_FUNC),
    ("response_heavy_rank", "Response-heavy", C_RESP),
]
plot_genes = ranks.sort_values("consensus_rank").head(15).copy()
xpos = np.arange(4)

# draw trajectories
for _, r in plot_genes.iterrows():
    g = r["gene_symbol"]
    vals = [float(r[c]) for c,_,_ in models]
    is_priority = g in priority5
    axB.plot(
        xpos, vals,
        lw=3.2 if is_priority else 1.8,
        alpha=1.0 if is_priority else 0.42,
        color=C_PRIORITY if is_priority else "#8A8A8A",
        marker="o", markersize=7 if is_priority else 4,
        zorder=4 if is_priority else 2
    )

# deterministic right-side label de-overlap
end = plot_genes[["gene_symbol","response_heavy_rank"]].copy()
end["y_raw"] = end["response_heavy_rank"].astype(float)
end = end.sort_values("y_raw").reset_index(drop=True)

min_gap = 2.0
label_y = []
for i, row in end.iterrows():
    y0 = float(row["y_raw"])
    if i == 0:
        label_y.append(y0)
    else:
        label_y.append(max(y0, label_y[-1] + min_gap))

# If the label stack exceeds the axis, compress gently while preserving order
max_allowed = 48.0
if label_y[-1] > max_allowed:
    shift = label_y[-1] - max_allowed
    label_y = [y - shift for y in label_y]
    for i in range(len(label_y)-2, -1, -1):
        label_y[i] = min(label_y[i], label_y[i+1] - min_gap)

for (idx, row), ylab in zip(end.iterrows(), label_y):
    g = row["gene_symbol"]
    yraw = float(row["y_raw"])
    is_priority = g in priority5
    axB.plot([3.0, 3.12], [yraw, ylab],
             color=C_PRIORITY if is_priority else "#8A8A8A",
             lw=1.6 if is_priority else 1.0,
             alpha=0.9 if is_priority else 0.65,
             clip_on=False)
    axB.text(
        3.15, ylab, g,
        va="center", ha="left",
        fontsize=17.0 if is_priority else 15.0,
        fontweight="bold",
        color=C_PRIORITY if is_priority else "#555555",
        clip_on=False
    )

axB.set_xticks(xpos)
axB.set_xticklabels([x[1] for x in models], rotation=15, ha="right", fontsize=18, fontweight="bold")
axB.set_ylabel("Rank (lower = higher priority)")
axB.invert_yaxis()
axB.grid(axis="y", color=C_GRID, lw=1.05)
axB.set_axisbelow(True)
axB.spines[["top","right"]].set_visible(False)
axB.set_xlim(-0.2, 3.95)
axB.set_ylim(50, -1)

# -------------------- C: rank-correlation matrix --------------------
labels = ["Primary","Independent-heavy","Functional-heavy","Response-heavy"]
keys = ["primary","independent_heavy","functional_heavy","response_heavy"]
M = np.eye(4)

for _, r in corrs.iterrows():
    m1, m2 = str(r["model_1"]), str(r["model_2"])
    if m1 in keys and m2 in keys:
        i, j = keys.index(m1), keys.index(m2)
        rho = float(r["spearman_rho"])
        M[i,j] = rho
        M[j,i] = rho

im = axC.imshow(M, vmin=0.97, vmax=1.0, cmap="Blues", aspect="equal")
axC.set_xticks(range(4))
axC.set_xticklabels(labels, rotation=22, ha="right", fontsize=18, fontweight="bold")
axC.set_yticks(range(4))
axC.set_yticklabels(labels, fontsize=18, fontweight="bold")
for i in range(4):
    for j in range(4):
        axC.text(j, i, f"{M[i,j]:.3f}", ha="center", va="center",
                 fontsize=20.0, fontweight="bold",
                 color="white" if M[i,j] > 0.989 else C_TEXT)
cb = fig.colorbar(im, ax=axC, fraction=0.046, pad=0.04)
cb.set_label("Spearman ρ", fontweight="bold", fontsize=21)
cb.ax.tick_params(labelsize=17, width=1.8)

# -------------------- D: robustness summary --------------------
r = ranks.sort_values("consensus_rank").head(30).copy()
r = r.sort_values(["top10_models_n","top20_models_n","consensus_rank"],
                  ascending=[False,False,True]).head(15)
r = r.sort_values("consensus_rank", ascending=False)

y = np.arange(len(r))
axD.barh(y, r["top20_models_n"], color="#BBD7EA", height=0.72,
         label="Top 20")
axD.barh(y, r["top10_models_n"], color=C_PRIMARY, height=0.72,
         label="Top 10")
axD.set_yticks(y)
axD.set_yticklabels(r["gene_symbol"], fontsize=18, fontweight="bold")
axD.set_xlim(0,4.32)
axD.set_xticks([0,1,2,3,4])
axD.set_xlabel("Weighting models retaining gene")
axD.grid(axis="x", color=C_GRID, lw=1.05)
axD.set_axisbelow(True)
axD.spines[["top","right"]].set_visible(False)

for yi, (_, rr) in enumerate(r.iterrows()):
    if rr["gene_symbol"] in priority5:
        axD.text(4.08, yi, "●", va="center", ha="center",
                 fontsize=23, color=C_PRIORITY, clip_on=False)

legD = axD.legend(frameon=False, loc="lower right", fontsize=16.0)
for t in legD.get_texts():
    t.set_fontweight("bold")

# -------------------- panel headers --------------------
titles = {
    "A": "Full evidence-weighted top-30 ranking",
    "B": "Priority ranks across alternative weights",
    "C": "Rank correlation across weighting schemes",
    "D": "Priority-gene retention across sensitivity models",
}
axes = {"A":axA,"B":axB,"C":axC,"D":axD}
for lab, ax in axes.items():
    box = ax.get_position()
    fig.text(box.x0-0.028, box.y1+0.026, lab,
             ha="left", va="bottom", fontsize=26, fontweight="bold")
    fig.text(box.x0+0.012, box.y1+0.026, titles[lab],
             ha="left", va="bottom", fontsize=24.0, fontweight="bold")

fig.savefig(PDF, bbox_inches="tight")
fig.savefig(PNG, dpi=600, bbox_inches="tight")
fig.savefig(SVG, bbox_inches="tight")
plt.close(fig)

# -------------------- provenance --------------------
def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

manifest = PROV / "Supplementary_Figure_S9_SHA256.txt"
with open(manifest, "w") as f:
    for p in [MASTER, TOP30, RANKS, CORRS, PDF, PNG, SVG]:
        f.write(f"{sha256(p)}  {p}\n")

print("=== MANUSCRIPT SUPPLEMENTARY FIGURE S9 GENERATED ===")
print(PDF)
print(PNG)
print(SVG)
print(manifest)
