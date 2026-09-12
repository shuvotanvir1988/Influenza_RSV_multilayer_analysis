from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

import hashlib

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 16.0,
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

ROOT = Path.home() / "Influenza_RSV_Project"

FLU = (
    ROOT /
    "results/regulatory_driver_analysis/bulk/differential_activity/"
    "GPL6884_InfluenzaA_vs_control_TF_activity_limma_v1.0.tsv"
)

MASTER = (
    ROOT /
    "results/regulatory_driver_analysis/tables/"
    "REGULATORY_DRIVER_EVIDENCE_MASTER_v1.0.tsv"
)

DRIVERS37 = (
    ROOT /
    "results/regulatory_driver_analysis/tables/"
    "REGULATORY_HIGH_CONFIDENCE_DRIVERS_v1.0.tsv"
)

BIASED13 = (
    ROOT /
    "results/regulatory_driver_analysis/tables/"
    "REGULATORY_HIGH_CONFIDENCE_PATHOGEN_BIASED_DRIVERS_v1.0.tsv"
)

OUTDIR = ROOT / "results/figure5_influenza_centered"
FIGDIR = OUTDIR / "figures"
TABDIR = OUTDIR / "tables"

FIGDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)

OUTPDF = FIGDIR / "Figure5_SUBMISSION_v14_influenza_regulatory_architecture.pdf"
OUTPNG = FIGDIR / "Figure5_SUBMISSION_v14_influenza_regulatory_architecture.png"
OUTSVG = FIGDIR / "Figure5_SUBMISSION_v14_influenza_regulatory_architecture.svg"

# ============================================================
# Load frozen data
# ============================================================

flu = pd.read_csv(FLU, sep="\t")
master = pd.read_csv(MASTER, sep="\t")
drivers = pd.read_csv(DRIVERS37, sep="\t")
biased = pd.read_csv(BIASED13, sep="\t")

# ============================================================
# QC
# ============================================================

n_tested = len(flu)
n_flu_sig = int((flu["FDR"] < 0.05).sum())
n_flu_up = int(((flu["FDR"] < 0.05) &
                (flu["delta_activity"] > 0)).sum())
n_flu_down = int(((flu["FDR"] < 0.05) &
                  (flu["delta_activity"] < 0)).sum())
n_direct_sig = int((master["flu_vs_rsv_FDR"] < 0.05).sum())

print("=== FIGURE 5 v3 INPUT CHECK ===")
print("TFs tested:", n_tested)
print("Influenza vs control FDR<0.05:", n_flu_sig)
print("  positive:", n_flu_up)
print("  negative:", n_flu_down)
print("Influenza vs RSV FDR<0.05:", n_direct_sig)
print("Frozen high-confidence drivers:", len(drivers))
print("Frozen pathogen-biased drivers:", len(biased))

assert n_tested == 761
assert n_flu_sig == 294
assert len(drivers) == 37
assert len(biased) == 13

# ============================================================
# Common sets and ordering
# ============================================================

driver37 = set(drivers["TF"])
biased13 = set(biased["TF"])

# Order focused drivers primarily by regulatory-driver score,
# then magnitude of influenza-vs-RSV difference.
focused = biased.copy()
focused["abs_flu_vs_rsv"] = focused["flu_vs_rsv_delta"].abs()

focused = focused.sort_values(
    ["regulatory_driver_score", "abs_flu_vs_rsv", "flu_delta"],
    ascending=[False, False, False]
).reset_index(drop=True)

# Keep biologically central interferon regulators visually prominent
preferred_front = [
    "STAT1",
    "IRF1",
    "STAT3",
    "IRF3",
    "STAT2"
]

front = [x for x in preferred_front if x in set(focused["TF"])]
rest = [x for x in focused["TF"] if x not in front]

order = front + rest

focused = focused.set_index("TF").loc[order].reset_index()

print()
print("Focused 13-driver order:")
for i, tf in enumerate(focused["TF"], 1):
    print(f"{i:2d}. {tf}")

focused.to_csv(
    TABDIR / "Figure5_SUBMISSION_v14_focused13_driver_table.tsv",
    sep="\t",
    index=False
)

# ============================================================
# Figure
# ============================================================

fig = plt.figure(figsize=(20.2, 15.0))

gs = fig.add_gridspec(
    2, 2,
    width_ratios=[1.00, 1.06],
    height_ratios=[1.0, 0.96],
    wspace=0.42,
    hspace=0.56,
    left=0.075,
    right=0.950,
    top=0.92,
    bottom=0.085
)

# ============================================================
# A. Influenza regulatory landscape
# ============================================================

ax = fig.add_subplot(gs[0, 0])

sig = flu["FDR"] < 0.05
y = -np.log10(flu["FDR"].clip(lower=1e-300))

# background
ax.scatter(
    flu.loc[~sig, "delta_activity"],
    y.loc[~sig],
    s=28,
    alpha=0.30,
    rasterized=True
)

# significant
ax.scatter(
    flu.loc[sig, "delta_activity"],
    y.loc[sig],
    s=42,
    alpha=0.72,
    rasterized=True
)

# frozen 37
mask37 = flu["TF"].isin(driver37)

ax.scatter(
    flu.loc[mask37, "delta_activity"],
    y.loc[mask37],
    s=100,
    marker="D",
    facecolors="none",
    edgecolors="black",
    linewidths=1.5,
    label="Frozen 37 drivers"
)

# biased 13
mask13 = flu["TF"].isin(biased13)

ax.scatter(
    flu.loc[mask13, "delta_activity"],
    y.loc[mask13],
    s=130,
    marker="o",
    facecolors="none",
    edgecolors="black",
    linewidths=2.1,
    label="Pathogen-biased drivers"
)

# cleaner annotation set
labels_A = [
    "STAT1",
    "STAT2",
    "STAT3",
    "IRF1",
    "IRF3",
    "IRF7",
    "IRF9",
    "NFKB",
    "SPI1"
]

offsets_A = {
    "STAT1": (6, 3),
    "STAT2": (5, 4),
    "STAT3": (5, 5),
    "IRF1": (4, 4),
    "IRF3": (5, -10),
    "IRF7": (5, 4),
    "IRF9": (5, 4),
    "NFKB": (5, 4),
    "SPI1": (5, 4),
}

for tf in labels_A:
    z = flu.loc[flu["TF"] == tf]
    if len(z) == 1:
        xv = float(z["delta_activity"].iloc[0])
        yv = -np.log10(max(float(z["FDR"].iloc[0]), 1e-300))

        ax.annotate(
            tf,
            (xv, yv),
            xytext=offsets_A.get(tf, (4, 4)),
            textcoords="offset points",
            fontsize=15.0,
            fontweight="bold"
        )

ax.axhline(
    -np.log10(0.05),
    linestyle="--",
    linewidth=1.55,
    color="0.28"
)
ax.axvline(0, linewidth=1.55, color="0.18")

ax.set_xlabel(
    "Δ inferred TF activity\n"
    "Influenza A acute − healthy control",
    fontweight="bold"
)
ax.set_ylabel("−log10(FDR)", fontweight="bold")

# Panel header added later at figure level.
# Panel A summary added later at figure level in a dedicated header band.

# Panel A legend is added later at figure level, outside the plotting axes.
handles_A, labels_A = ax.get_legend_handles_labels()
ax.tick_params(axis="both", width=1.55, length=5.8)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontweight("bold")
ax.spines[["top", "right"]].set_visible(False)

# ============================================================
# B. Influenza–RSV regulatory geometry
# ============================================================

ax = fig.add_subplot(gs[0, 1])

direct = master["flu_vs_rsv_FDR"] < 0.05

# all regulators
ax.scatter(
    master.loc[~direct, "rsv_delta"],
    master.loc[~direct, "flu_delta"],
    s=29,
    alpha=0.29,
    rasterized=True
)

# significant direct contrast
ax.scatter(
    master.loc[direct, "rsv_delta"],
    master.loc[direct, "flu_delta"],
    s=44,
    alpha=0.72,
    label=f"Influenza vs RSV FDR < 0.05 (n={direct.sum()})",
    rasterized=True
)

# frozen 37
m37 = master["TF"].isin(driver37)

ax.scatter(
    master.loc[m37, "rsv_delta"],
    master.loc[m37, "flu_delta"],
    s=100,
    marker="D",
    facecolors="none",
    edgecolors="black",
    linewidths=1.5,
    label="Frozen 37 drivers"
)

# focused 13
m13 = master["TF"].isin(biased13)

ax.scatter(
    master.loc[m13, "rsv_delta"],
    master.loc[m13, "flu_delta"],
    s=132,
    marker="o",
    facecolors="none",
    edgecolors="black",
    linewidths=2.1,
    label="Pathogen-biased drivers"
)

# symmetric limits
xmin = min(
    master["rsv_delta"].min(),
    master["flu_delta"].min()
)
xmax = max(
    master["rsv_delta"].max(),
    master["flu_delta"].max()
)

pad = 0.08 * (xmax - xmin)
lims = [xmin - pad, xmax + pad]

ax.plot(
    lims,
    lims,
    linestyle="--",
    linewidth=1.55,
    color="0.30"
)

ax.set_xlim(lims)
ax.set_ylim(lims)

ax.axhline(0, linewidth=1.45, color="0.18")
ax.axvline(0, linewidth=1.45, color="0.18")

# clearer key regulator labels
labels_B = [
    "STAT1",
    "STAT2",
    "STAT3",
    "IRF1",
    "IRF3",
    "HMGA1",
    "REL",
    "ATF3"
]

offsets_B = {
    "STAT1": (12, 10),
    "STAT2": (-52, 24),
    "STAT3": (24, -16),
    "IRF1": (-58, 15),
    "IRF3": (-52, -24),
    "HMGA1": (-62, 15),
    "REL": (25, 17),
    "ATF3": (-34, -27),
}

for tf in labels_B:
    z = master.loc[master["TF"] == tf]

    if len(z) == 1:
        xv = float(z["rsv_delta"].iloc[0])
        yv = float(z["flu_delta"].iloc[0])

        dx, dy = offsets_B.get(tf, (8, 8))
        ax.annotate(
            tf,
            (xv, yv),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=13.5,
            fontweight="bold",
            ha="left" if dx >= 0 else "right",
            va="bottom" if dy >= 0 else "top",
            arrowprops=dict(
                arrowstyle="-",
                color="0.35",
                lw=0.85,
                shrinkA=1,
                shrinkB=2
            ),
            bbox=dict(
                boxstyle="round,pad=0.08",
                facecolor="white",
                edgecolor="none",
                alpha=0.80
            ),
            zorder=8
        )

# Biological interpretation of diagonal
xpos = lims[0] + 0.07 * (lims[1] - lims[0])
ypos = lims[1] - 0.09 * (lims[1] - lims[0])

ax.text(
    lims[0] + 0.03 * (lims[1] - lims[0]),
    lims[1] - 0.11 * (lims[1] - lims[0]),
    "Influenza-amplified",
    fontsize=15.0,
    fontweight="bold",
    va="top",
    bbox=dict(
        boxstyle="round,pad=0.10",
        facecolor="white",
        edgecolor="none",
        alpha=0.84
    )
)

xpos2 = lims[1] - 0.08 * (lims[1] - lims[0])
ypos2 = lims[0] + 0.24 * (lims[1] - lims[0])

ax.text(
    lims[1] - 0.05 * (lims[1] - lims[0]),
    lims[0] + 0.25 * (lims[1] - lims[0]),
    "RSV-amplified",
    fontsize=15.0,
    fontweight="bold",
    ha="right",
    va="bottom",
    bbox=dict(
        boxstyle="round,pad=0.10",
        facecolor="white",
        edgecolor="none",
        alpha=0.84
    )
)

ax.set_xlabel(
    "RSV acute − healthy control\n"
    "Δ inferred TF activity",
    fontweight="bold"
)

ax.set_ylabel(
    "Influenza A acute − healthy control\n"
    "Δ inferred TF activity",
    fontweight="bold"
)

# Panel header added later at figure level.
# Panel B summary added later at figure level in a dedicated header band.

legB = ax.legend(
    frameon=True,
    fontsize=13.0,
    loc="lower right",
    bbox_to_anchor=(0.985, 0.035),
    facecolor="white",
    edgecolor="0.72",
    framealpha=1.0,
    borderpad=0.40,
    labelspacing=0.30,
    handletextpad=0.40
)
for t in legB.get_texts():
    t.set_fontweight("bold")
ax.tick_params(axis="both", width=1.55, length=5.8)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontweight("bold")
ax.spines[["top", "right"]].set_visible(False)

# ============================================================
# C. Focused 13-driver regulatory contrast heatmap
# ============================================================

ax = fig.add_subplot(gs[1, 0])

contrast = focused.set_index("TF")[
    [
        "flu_delta",
        "rsv_delta",
        "flu_vs_rsv_delta"
    ]
].copy()

contrast.columns = [
    "Influenza\n− control",
    "RSV\n− control",
    "Influenza\n− RSV"
]

vals = contrast.to_numpy(dtype=float)

vmax = np.nanmax(np.abs(vals))
norm = TwoSlopeNorm(
    vmin=-vmax,
    vcenter=0,
    vmax=vmax
)

im = ax.imshow(
    vals,
    aspect="auto",
    interpolation="nearest",
    cmap="coolwarm",
    norm=norm
)

ax.set_xticks(range(3))
ax.set_xticklabels(
    contrast.columns,
    rotation=0,
    ha="center",
    fontsize=15.5,
    fontweight="bold"
)

ax.set_yticks(range(len(contrast)))
ax.set_yticklabels(
    contrast.index,
    fontsize=17.0,
    fontweight="bold"
)

# exact values in every cell
for i in range(vals.shape[0]):
    for j in range(vals.shape[1]):
        val = vals[i, j]

        ax.text(
            j,
            i,
            f"{val:.2f}",
            ha="center",
            va="center",
            fontsize=15.0,
            fontweight="bold"
        )

cb = fig.colorbar(
    im,
    ax=ax,
    fraction=0.034,
    pad=0.025
)

cb.set_label("Δ inferred TF activity", fontsize=17.0, fontweight="bold")

# Panel header added later at figure level.
ax.tick_params(axis="both", width=1.45, length=5.4)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontweight("bold")
for tick in cb.ax.get_yticklabels():
    tick.set_fontsize(15.0)
    tick.set_fontweight("bold")

# ============================================================
# D. Gene-program connectivity, identical row order
# ============================================================

ax = fig.add_subplot(gs[1, 1])

conn = focused.set_index("TF")[
    [
        "shared_core_targets",
        "shared_influenza_amplified_targets",
        "influenza_selective_targets"
    ]
].copy()

conn.columns = [
    "Shared core\n(n=130)",
    "Influenza-amplified\nshared\n(n=37)",
    "Influenza-\nselective"
]

# SAME ORDER AS PANEL C
conn = conn.loc[contrast.index]

cvals = conn.to_numpy(dtype=float)

im2 = ax.imshow(
    np.sqrt(cvals),
    aspect="auto",
    interpolation="nearest",
    cmap="viridis"
)

ax.set_xticks(range(3))
ax.set_xticklabels(
    conn.columns,
    fontsize=15.5,
    fontweight="bold"
)

ax.set_yticks(range(len(conn)))
ax.set_yticklabels(
    conn.index,
    fontsize=17.0,
    fontweight="bold"
)

for i in range(cvals.shape[0]):
    for j in range(cvals.shape[1]):

        value = int(cvals[i, j])

        ax.text(
            j,
            i,
            str(value),
            ha="center",
            va="center",
            fontsize=15.0,
            fontweight="bold"
        )

cb2 = fig.colorbar(
    im2,
    ax=ax,
    fraction=0.034,
    pad=0.045
)

cb2.set_label("√ number of direct CollecTRI targets", fontsize=15.5, fontweight="bold", labelpad=12)

# Panel header added later at figure level.
ax.tick_params(axis="both", width=1.45, length=5.4)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontweight("bold")
for tick in cb2.ax.get_yticklabels():
    tick.set_fontsize(15.0)
    tick.set_fontweight("bold")

# ============================================================
# PANEL HEADERS — LOCKED FIGURE 2/3/4 SUBMISSION STYLE
# ============================================================

PANEL_LABEL_SIZE = 25.0
PANEL_TITLE_SIZE = 23.0

panel_axes = [fig.axes[0], fig.axes[1], fig.axes[2], fig.axes[4]]
panel_titles = {
    "A": "Influenza regulatory landscape",
    "B": "Influenza amplification within a shared\nantiviral regulatory response",
    "C": "Regulatory contrast architecture of the\n13 pathogen-biased drivers",
    "D": "Regulatory connectivity to the\ninfluenza-centered gene architecture",
}

fig.canvas.draw()
for label, axis in zip(["A", "B", "C", "D"], panel_axes):
    box = axis.get_position()
    y = box.y1 + 0.022
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
        va="bottom",
        linespacing=1.05
    )

# Dedicated top-row summary band: outside plotting axes, below panel titles.
# This prevents any collision with TF labels, data, legends, or axes.
boxA = panel_axes[0].get_position()
boxB = panel_axes[1].get_position()

summary_y_A = boxA.y1 + 0.004
summary_y_B = boxB.y1 + 0.004

fig.text(
    boxA.x0 + 0.005,
    summary_y_A,
    f"{n_flu_sig}/{n_tested} TFs significant "
    f"({n_flu_up} higher; {n_flu_down} lower)",
    ha="left",
    va="bottom",
    fontsize=14.5,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.16",
        facecolor="white",
        edgecolor="0.78",
        linewidth=0.65,
        alpha=1.0
    )
)

fig.text(
    boxB.x0 + 0.005,
    summary_y_B,
    f"{n_direct_sig}/{n_tested} TFs differ between influenza and RSV",
    ha="left",
    va="bottom",
    fontsize=14.5,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.16",
        facecolor="white",
        edgecolor="0.78",
        linewidth=0.65,
        alpha=1.0
    )
)

# Panel A legend in dedicated whitespace outside the graph.
# This placement is intentionally independent of data coordinates.
legendA = panel_axes[0].legend(
    handles_A,
    labels_A,
    loc="lower right",
    bbox_to_anchor=(0.985, 0.035),
    bbox_transform=panel_axes[0].transAxes,
    ncol=1,
    fontsize=13.0,
    frameon=True,
    facecolor="white",
    edgecolor="0.72",
    framealpha=1.0,
    borderpad=0.42,
    labelspacing=0.34,
    handletextpad=0.44
)
for t in legendA.get_texts():
    t.set_fontweight("bold")

# No overall Figure X title and no footer.
# Detailed methodological definitions belong in the manuscript legend.

# ============================================================
# FIGURE-FACING PROVENANCE
# ============================================================

focused.to_csv(
    TABDIR / "Figure5_SUBMISSION_v14_focused13_driver_table.tsv",
    sep="\t",
    index=False
)

manifest_rows = []
for p in [FLU, MASTER, DRIVERS37, BIASED13]:
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    manifest_rows.append({"file": str(p), "sha256": h})

pd.DataFrame(manifest_rows).to_csv(
    TABDIR / "Figure5_SUBMISSION_v14_provenance_sha256.tsv",
    sep="\t",
    index=False
)

plt.subplots_adjust(
    top=0.895,
    bottom=0.105
)

fig.savefig(
    OUTPDF,
    bbox_inches="tight"
)

fig.savefig(
    OUTPNG,
    dpi=600,
    bbox_inches="tight"
)

fig.savefig(
    OUTSVG,
    format="svg",
    bbox_inches="tight"
)

print()
print("Written:")
print(OUTPDF)
print(OUTPNG)
print(OUTSVG)

print()
print("STATUS: FIGURE 5 SUBMISSION v14 COMPLETE")
