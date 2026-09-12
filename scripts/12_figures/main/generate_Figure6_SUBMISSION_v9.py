#!/usr/bin/env python3

from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import TwoSlopeNorm

ROOT = Path.home() / "Influenza_RSV_Project"

MASTER = ROOT / "results/rnaseq_validation/influenza_SIG/validation/SIG_four_season_frozen170_validation_master.tsv"
CROSS = ROOT / "results/rnaseq_validation/influenza_SIG/validation/SIG_four_season_frozen170_crossseason_summary.tsv"
SEASON = ROOT / "results/rnaseq_validation/influenza_SIG/validation/SIG_four_season_frozen170_seasonal_summary.tsv"
CORR = ROOT / "results/rnaseq_validation/influenza_SIG/validation/SIG_four_season_frozen170_effect_correlation.tsv"

OUTDIR = ROOT / "results/figure6_influenza_validation"
FIGDIR = OUTDIR / "figures"
TABLEDIR = OUTDIR / "tables"
LOGDIR = OUTDIR / "logs"
PROVDIR = OUTDIR / "provenance"

for d in [FIGDIR, TABLEDIR, LOGDIR, PROVDIR]:
    d.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(MASTER, sep="\t")
cross = pd.read_csv(CROSS, sep="\t")
season = pd.read_csv(SEASON, sep="\t")
corr = pd.read_csv(CORR, sep="\t")

class_colors = {
    "Shared_core": "#7A7A7A",
    "Influenza_amplified_shared": "#D95F02",
    "RSV_amplified_shared": "#4C78A8",
}

class_labels = {
    "Shared_core": "Shared core",
    "Influenza_amplified_shared": "Influenza-amplified",
    "RSV_amplified_shared": "RSV-amplified",
}

season_map = {
    "2018": ("GSE158592", 74),
    "2019": ("GSE155635", 71),
    "2020": ("GSE196350", 42),
    "2022": ("GSE213168", 21),
}

# Reuse frozen deterministic heatmap selection from v2/v1
selection_file = TABLEDIR / "Figure6_v2_heatmap_gene_selection.tsv"
if not selection_file.exists():
    selection_file = TABLEDIR / "Figure6_v1_heatmap_gene_selection.tsv"

selected = pd.read_csv(selection_file, sep="\t")

selected.to_csv(
    TABLEDIR / "Figure6_SUBMISSION_v9_heatmap_gene_selection.tsv",
    sep="\t",
    index=False
)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 15.0,
    "axes.titlesize": 22.0,
    "axes.labelsize": 18.0,
    "axes.labelweight": "bold",
    "xtick.labelsize": 15.0,
    "ytick.labelsize": 15.0,
    "legend.fontsize": 14.0,
    "axes.linewidth": 1.85,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# No figure-level title/footer in manuscript-facing version
fig = plt.figure(figsize=(18.4, 12.8))

gs = gridspec.GridSpec(
    2, 2,
    width_ratios=[1.02, 1.13],
    height_ratios=[0.82, 1.18],
    left=0.080,
    right=0.955,
    top=0.905,
    bottom=0.135,
    hspace=0.64,
    wspace=0.43
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

# ============================================================
# A — validation design
# ============================================================
axA.axis("off")
# Panel header added later at figure level.

boxes = [
    (0.10, "Discovery\nmicroarray", "#F1F1F1", "#666666"),
    (0.40, "Frozen\n170-gene\narchitecture", "#FFF3E6", "#D95F02"),
    (0.72, "RNA-seq\nvalidation", "#EEF4F8", "#4C78A8"),
]

for x, text, fc, ec in boxes:
    axA.text(
        x, 0.73, text,
        ha="center", va="center",
        fontsize=16.0, fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.44",
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.8
        )
    )

axA.annotate(
    "", xy=(0.275, 0.73), xytext=(0.195, 0.73),
    arrowprops=dict(arrowstyle="->", lw=2.2, color="#444444")
)
axA.annotate(
    "", xy=(0.585, 0.73), xytext=(0.495, 0.73),
    arrowprops=dict(arrowstyle="->", lw=2.2, color="#444444")
)

# Larger season table
y0 = 0.46
dy = 0.095
for i, (yr, (gse, n)) in enumerate(season_map.items()):
    y = y0 - i * dy
    axA.text(
        0.54, y, yr,
        ha="left", va="center",
        fontsize=15.0, fontweight="bold"
    )
    axA.text(
        0.68, y, gse,
        ha="left", va="center",
        fontsize=14.0,
        fontweight="bold"
    )
    axA.text(
        0.99, y, f"n={n}",
        ha="right", va="center",
        fontsize=14.0,
        fontweight="bold"
    )

axA.text(
    0.17, 0.205,
    "Validation manifest\n208 samples\n127 influenza | 81 healthy",
    ha="center", va="center",
    fontsize=14.5,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.40",
        facecolor="white",
        edgecolor="#777777",
        linewidth=1.4
    )
)


axA.set_xlim(0, 1)
axA.set_ylim(0, 1)

# ============================================================
# B — directional replication
# ============================================================
# Panel header added later at figure level.

classes = [
    "ALL_170",
    "Shared_core",
    "Influenza_amplified_shared",
    "RSV_amplified_shared",
]

display = [
    "All frozen genes",
    "Shared core",
    "Influenza-amplified",
    "RSV-amplified",
]

vals, nums = [], []

for c in classes:
    row = cross[cross["class"] == c].iloc[0]
    vals.append(100 * row["direction_match_all_4_fraction"])
    nums.append(
        f'{int(row["direction_match_all_4"])}/'
        f'{int(row["evaluable_all_4_seasons"])}'
    )

bar_colors = [
    "#555555",
    class_colors["Shared_core"],
    class_colors["Influenza_amplified_shared"],
    class_colors["RSV_amplified_shared"],
]

y = np.arange(len(classes))
bars = axB.barh(
    y, vals,
    color=bar_colors,
    height=0.66,
    edgecolor="0.25",
    linewidth=0.55
)

axB.set_yticks(y)
axB.set_yticklabels(display, fontweight="bold")
axB.invert_yaxis()
axB.set_xlim(0, 108)
axB.set_xlabel("Direction-concordant genes in all four seasons (%)", fontweight="bold")

for b, v, n in zip(bars, vals, nums):
    axB.text(
        min(v + 1.2, 103.2),
        b.get_y() + b.get_height()/2,
        f"{n} ({v:.1f}%)",
        va="center",
        fontsize=14.0,
        fontweight="bold"
    )

axB.text(
    0.50, -0.20,
    "RSV-amplified class shown descriptively; n=2 evaluable genes.",
    transform=axB.transAxes,
    ha="center",
    va="top",
    fontsize=12.5,
    fontweight="bold",
    color="#444444",
    clip_on=False
)

axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
axB.tick_params(axis="both", width=1.55, length=5.8)
for tick in axB.get_xticklabels() + axB.get_yticklabels():
    tick.set_fontweight("bold")

# ============================================================
# C — quantitative effect replication
# ============================================================
# Panel header added later at figure level.

plotdf = df[
    df["flu_logFC"].notna() &
    df["mean_RNAseq_logFC"].notna()
].copy()

for cl in [
    "Shared_core",
    "RSV_amplified_shared",
    "Influenza_amplified_shared",
]:
    z = plotdf[plotdf["Figure7A_class"] == cl]

    axC.scatter(
        z["flu_logFC"],
        z["mean_RNAseq_logFC"],
        s=54 if cl == "Influenza_amplified_shared" else 40,
        alpha=0.82,
        color=class_colors[cl],
        edgecolor="white",
        linewidth=0.55,
        label=class_labels[cl]
    )

allcorr = corr[corr["class"] == "ALL_170"].iloc[0]

axC.text(
    0.035, 0.965,
    f'Pearson r = {allcorr["pearson_r_microarray_vs_mean_RNAseq"]:.3f}\n'
    f'Spearman ρ = {allcorr["spearman_rho_microarray_vs_mean_RNAseq"]:.3f}\n'
    f'n = {int(allcorr["evaluable"])} genes',
    transform=axC.transAxes,
    ha="left", va="top",
    fontsize=14.0,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.28",
        facecolor="white",
        edgecolor="#BBBBBB"
    )
)

xv = plotdf["flu_logFC"].to_numpy()
yv = plotdf["mean_RNAseq_logFC"].to_numpy()

coef = np.polyfit(xv, yv, 1)

xmin = xv.min() - 0.15
xmax = xv.max() + 0.15
xx = np.linspace(xmin, xmax, 200)

axC.plot(
    xx,
    coef[0] * xx + coef[1],
    color="#333333",
    lw=2.0,
    zorder=1
)

axC.axhline(0, color="#888888", lw=1.35)
axC.axvline(0, color="#888888", lw=1.35)

# Reduced final label set
label_offsets = {
    "IFI27": (5, 5),
    "OTOF": (5, -2),
    "IFI44L": (8, 7),
    "RSAD2": (8, -8),
    "IRF7": (-24, -9),
    "FCGR1A": (-34, 6),
}

for g, off in label_offsets.items():
    z = plotdf[plotdf["gene_symbol"] == g]
    if len(z) == 1:
        r = z.iloc[0]
        axC.annotate(
            g,
            (r["flu_logFC"], r["mean_RNAseq_logFC"]),
            xytext=off,
            textcoords="offset points",
            fontsize=13.0,
            fontweight="bold",
            arrowprops=dict(
                arrowstyle="-",
                lw=0.85,
                color="#777777"
            )
        )

axC.set_xlim(xmin, xmax)
axC.set_xlabel("Discovery influenza logFC", fontweight="bold")
axC.set_ylabel("Mean validation RNA-seq logFC", fontweight="bold")

legC = axC.legend(
    frameon=True,
    fontsize=13.0,
    loc="lower right",
    facecolor="white",
    edgecolor="0.80",
    framealpha=0.96,
    borderpad=0.38
)
for t in legC.get_texts():
    t.set_fontweight("bold")

axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)
axC.tick_params(axis="both", width=1.55, length=5.8)
for tick in axC.get_xticklabels() + axC.get_yticklabels():
    tick.set_fontweight("bold")

# ============================================================
# D — multiseason heatmap
# ============================================================
# Panel header added later at figure level.

heat_cols = [
    "flu_logFC",
    "logFC_2018",
    "logFC_2019",
    "logFC_2020",
    "logFC_2022",
]

heat_labels = [
    "Discovery",
    "2018",
    "2019",
    "2020",
    "2022",
]

selected = selected.copy()

selected["class_order"] = selected["Figure7A_class"].map({
    "Influenza_amplified_shared": 0,
    "Shared_core": 1,
})

selected = selected.sort_values(
    ["class_order", "flu_logFC"],
    ascending=[True, False]
)

mat = selected[heat_cols].to_numpy()

vmax = 3.0
mat_display = np.clip(mat, -vmax, vmax)
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

im = axD.imshow(
    mat_display,
    aspect="auto",
    cmap="RdBu_r",
    norm=norm,
    interpolation="nearest"
)

axD.set_xticks(np.arange(len(heat_labels)))
axD.set_xticklabels(heat_labels, fontweight="bold")

axD.set_yticks(np.arange(len(selected)))
axD.set_yticklabels(
    selected["gene_symbol"],
    fontsize=13.5,
    fontweight="bold"
)

n_infl = (
    selected["Figure7A_class"] ==
    "Influenza_amplified_shared"
).sum()

axD.axhline(
    n_infl - 0.5,
    color="black",
    lw=1.8
)

cbar = fig.colorbar(
    im,
    ax=axD,
    fraction=0.034,
    pad=0.03
)

cbar.set_label(
    "logFC (display clipped at ±3)",
    fontsize=14.5,
    fontweight="bold"
)
for tick in cbar.ax.get_yticklabels():
    tick.set_fontsize(13.0)
    tick.set_fontweight("bold")

# Compact one-line seasonal correlation annotation
all_season = season[
    season["class"] == "ALL_170"
].copy()

corr_items = [
    f'{int(r.season)} r={r.pearson_r:.2f}'
    for _, r in all_season.iterrows()
]

axD.text(
    0.5, -0.105,
    "Discovery vs season:  " + "  |  ".join(corr_items),
    transform=axD.transAxes,
    ha="center",
    va="top",
    fontsize=12.5,
    fontweight="bold",
    color="#444444"
)

# ============================================================
# PANEL HEADERS — LOCKED FIGURE 2–5 SUBMISSION STYLE
# ============================================================

PANEL_LABEL_SIZE = 25.0
PANEL_TITLE_SIZE = 22.0

panel_titles = {
    "A": "Frozen-signature validation across\nfour influenza seasons",
    "B": "Directional replication across\nall four seasons",
    "C": "Discovery effects quantitatively\nreplicate by RNA-seq",
    "D": "Representative genes retain effect\ndirection across seasons",
}

fig.canvas.draw()
for label, axis in zip(["A", "B", "C", "D"], [axA, axB, axC, axD]):
    box = axis.get_position()
    y = box.y1 + (0.018 if label in ["A", "B"] else 0.028)
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

# ============================================================
# FIGURE-FACING PROVENANCE
# ============================================================

prov_rows = []
for p in [MASTER, CROSS, SEASON, CORR, selection_file]:
    prov_rows.append({
        "file": str(p),
        "sha256": hashlib.sha256(p.read_bytes()).hexdigest()
    })

pd.DataFrame(prov_rows).to_csv(
    TABLEDIR / "Figure6_SUBMISSION_v9_provenance_sha256.tsv",
    sep="\t",
    index=False
)

pdf = FIGDIR / "Figure6_SUBMISSION_v9_influenza_multiseason_validation.pdf"
png = FIGDIR / "Figure6_SUBMISSION_v9_influenza_multiseason_validation.png"
svg = FIGDIR / "Figure6_SUBMISSION_v9_influenza_multiseason_validation.svg"

fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=600, bbox_inches="tight")
fig.savefig(svg, bbox_inches="tight")
plt.close(fig)

print("Written:")
print(pdf)
print(png)
print(svg)
print(TABLEDIR / "Figure6_SUBMISSION_v9_heatmap_gene_selection.tsv")

print("STATUS: FIGURE 6 SUBMISSION v9 COMPLETE")
