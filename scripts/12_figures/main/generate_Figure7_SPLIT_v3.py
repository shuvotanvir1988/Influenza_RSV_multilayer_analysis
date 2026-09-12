from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm

# ============================================================
# PATHS
# ============================================================

ROOT = Path.home() / "Influenza_RSV_Project"

OLDTAB = ROOT / "results/figure7_influenza_integration/tables"
S38 = ROOT / "results/session38_additional_analysis/frozen"

OUTDIR = ROOT / "results/session42_main_figure_editing/Figure7_SUBMISSION_SPLIT_v3"
FIGDIR = OUTDIR / "figures"
TABDIR = OUTDIR / "tables"
PROVDIR = OUTDIR / "provenance"

for d in [FIGDIR, TABDIR, PROVDIR]:
    d.mkdir(parents=True, exist_ok=True)

PDF = FIGDIR / "Figure7_SUBMISSION_SPLIT_v3_influenza_multilayer_integration.pdf"
PNG = FIGDIR / "Figure7_SUBMISSION_SPLIT_v3_influenza_multilayer_integration.png"
SVG = FIGDIR / "Figure7_SUBMISSION_SPLIT_v3_influenza_multilayer_integration.svg"

OLD_ORDERED = OLDTAB / "Figure7_SUBMISSION_v7_ordered_170_gene_evidence.tsv"
OLD_SUMMARY = OLDTAB / "Figure7_SUBMISSION_v7_evidence_summary.tsv"

A_TOP30 = (
    S38 / "Session38A/tables/"
    "Session38A_TOP30_EVIDENCE_WEIGHTED_GENES_v1.0.tsv"
)

B_STATS = (
    S38 / "Session38B/tables/"
    "Session38B_STATISTICAL_TESTS_v1.0.tsv"
)
B_BRIDGE = (
    S38 / "Session38B/tables/"
    "Session38B_CRISPR_BRIDGE_MOLECULAR_CONTEXT_v1.0.tsv"
)

C_PRIMARY = (
    S38 / "Session38C/tables/"
    "Session38C_ARCHITECTURE_LOCALIZATION_v1.0.tsv"
)
C_ROBUST = (
    S38 / "Session38C/tables/"
    "Session38C_ARCHITECTURE_DETECTION_ROBUSTNESS_v1.0.tsv"
)

D_CLASS = (
    S38 / "Session38D/tables/"
    "Session38D_STABILITY_CLASS_SUMMARY_v1.0.tsv"
)
D_ARCH = (
    S38 / "Session38D/tables/"
    "Session38D_ARCHITECTURE_CLASS_STABILITY_SUMMARY_v1.0.tsv"
)

inputs = [
    OLD_ORDERED, OLD_SUMMARY, A_TOP30, B_STATS, B_BRIDGE,
    C_PRIMARY, C_ROBUST, D_CLASS, D_ARCH,
]
for p in inputs:
    if not p.exists():
        raise FileNotFoundError(p)

# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(OLD_ORDERED, sep="\t")
summary = pd.read_csv(OLD_SUMMARY, sep="\t")
top30 = pd.read_csv(A_TOP30, sep="\t")
b_stats = pd.read_csv(B_STATS, sep="\t")
bridge = pd.read_csv(B_BRIDGE, sep="\t")
loc_primary = pd.read_csv(C_PRIMARY, sep="\t")
loc_robust = pd.read_csv(C_ROBUST, sep="\t")
stability = pd.read_csv(D_CLASS, sep="\t")
arch_stab = pd.read_csv(D_ARCH, sep="\t")

assert len(df) == 170 and df["gene_symbol"].nunique() == 170
assert len(top30) == 30
assert len(bridge) == 7
assert len(loc_primary) == 5
assert len(loc_robust) == 5
assert int(stability["n_genes"].sum()) == 170

priority5 = {"IFIH1", "IFIT3", "ISG15", "KPNB1", "STAT1"}
bridge7 = {"CHMP5", "FCGR1B", "HERC5", "HIST2H2AC", "KPNB1", "OTOF", "TOP2A"}

# ============================================================
# STYLE
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 20.0,
    "axes.titlesize": 23.0,
    "axes.labelsize": 22.0,
    "axes.labelweight": "bold",
    "xtick.labelsize": 19.0,
    "ytick.labelsize": 19.0,
    "legend.fontsize": 17.0,
    "axes.linewidth": 1.9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

C_DISC = "#4C78A8"
C_RNA = "#59A14F"
C_PROT = "#F28E2B"
C_REG = "#B279A2"
C_CRISPR = "#E15759"
C_SHARED = "#8FA6B8"
C_FLU = "#D76B45"
C_RSV = "#B8B8B8"

LIGHT = "#ECECEC"
MID = "#BDBDBD"
DARK = "#333333"

fig = plt.figure(figsize=(27.0, 19.5))
gs = fig.add_gridspec(
    3, 15,
    height_ratios=[1.24, 1.62, 1.34],
    left=0.070,
    right=0.955,
    top=0.955,
    bottom=0.080,
    wspace=1.25,
    hspace=0.78,
)

# Split v3 top row:
# Panel A gets slightly more width/height for larger boxes.
# Two effective columns of gutter remain between A and B.
# Panel B still starts far enough right for its y-axis labels.
axA = fig.add_subplot(gs[0, 0:6])
axB = fig.add_subplot(gs[0, 7:15])
axC = fig.add_subplot(gs[1, 0:15])
axD = fig.add_subplot(gs[2, 0:15])

# Dummy invisible axes allow the frozen E/F construction code below to
# execute unchanged; those panels are intentionally excluded from this figure.
axE = fig.add_axes([0.001, 0.001, 0.001, 0.001])
axF = fig.add_axes([0.001, 0.001, 0.001, 0.001])
axE.set_visible(False)
axF.set_visible(False)

# Extra top-row clearance: shift B slightly right without changing its height.
_posB = axB.get_position()
axB.set_position([
    _posB.x0 + 0.018,
    _posB.y0,
    _posB.width - 0.018,
    _posB.height,
])

# Reserve room beneath C and D for enlarged labels/legends.
for _ax, _dy in [(axC, 0.010), (axD, 0.008)]:
    _pos = _ax.get_position()
    _ax.set_position([_pos.x0, _pos.y0 + _dy, _pos.width, _pos.height - _dy])

# ============================================================
# PANEL A — PRESERVED FROZEN v7 FRAMEWORK
# ============================================================

axA.axis("off")
axA.set_xlim(0, 1)
axA.set_ylim(0, 1)

def box(ax, xy, w, h, text, fc, fontsize=15.0, weight="normal"):
    x, y = xy
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.015,rounding_size=0.025",
        linewidth=1.65,
        edgecolor="white",
        facecolor=fc,
    )
    ax.add_patch(p)
    ax.text(
        x + w / 2, y + h / 2, text,
        ha="center", va="center",
        fontsize=fontsize, weight=weight, color="white",
    )

box(
    axA, (0.06, 0.71), 0.88, 0.25,
    "Discovery architecture\n170 genes",
    C_DISC, 19.5, "bold",
)

layers = [
    ("4-season\nRNA-seq\n155/161 supported\n≥1 season", C_RNA),
    ("SomaScan\nproteomics\n22/84 strong", C_PROT),
    ("Regulatory\nconnections\n84/170 genes", C_REG),
    ("IAV CRISPR\ndependency\n7/152 genes", C_CRISPR),
]
evidence_positions = [
    (0.02, 0.37),
    (0.515, 0.37),
    (0.02, 0.015),
    (0.515, 0.015),
]
box_w, box_h = 0.465, 0.305

for (label, color), (x, y) in zip(layers, evidence_positions):
    box(axA, (x, y), box_w, box_h, label, color, 17.5, "bold")

left_center = evidence_positions[0][0] + box_w / 2
right_center = evidence_positions[1][0] + box_w / 2
spine_x = 0.50
upper_branch_y = 0.685
lower_branch_y = 0.335

axA.plot([spine_x, spine_x], [0.71, lower_branch_y],
         color="#777777", lw=1.9, solid_capstyle="round",
         clip_on=False, zorder=1)
axA.plot([left_center, right_center], [upper_branch_y, upper_branch_y],
         color="#777777", lw=1.9, solid_capstyle="round",
         clip_on=False, zorder=1)
for x, y in evidence_positions[:2]:
    xc = x + box_w / 2
    axA.annotate(
        "", xy=(xc, y + box_h + 0.003), xytext=(xc, upper_branch_y),
        arrowprops=dict(arrowstyle="-|>", lw=1.9, color="#777777",
                        shrinkA=0, shrinkB=2),
        zorder=1,
    )
axA.plot([left_center, right_center], [lower_branch_y, lower_branch_y],
         color="#777777", lw=1.9, solid_capstyle="round",
         clip_on=False, zorder=1)
for x, y in evidence_positions[2:]:
    xc = x + box_w / 2
    axA.annotate(
        "", xy=(xc, y + box_h + 0.003), xytext=(xc, lower_branch_y),
        arrowprops=dict(arrowstyle="-|>", lw=1.9, color="#777777",
                        shrinkA=0, shrinkB=2),
        zorder=1,
    )

# ============================================================
# PANEL B — PRESERVED FROZEN v7 170-GENE EVIDENCE LANDSCAPE
# ============================================================

matrix_cols = ["RNAseq", "Proteomics", "Regulatory", "CRISPR"]
eval_map = {
    "RNAseq": "RNAseq_evaluable",
    "Proteomics": "Proteomics_evaluable",
    "Regulatory": "Regulatory_evaluable",
    "CRISPR": "CRISPR_evaluable",
}

M = np.zeros((len(df), len(matrix_cols)))
for i, row in df.iterrows():
    for j, col in enumerate(matrix_cols):
        if not bool(row[eval_map[col]]):
            M[i, j] = -1
        elif bool(row[col]):
            M[i, j] = 1
        else:
            M[i, j] = 0

cmap = ListedColormap(["#F2F2F2", "#FFFFFF", "#3A3A3A"])
norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)

axB.imshow(
    M.T, aspect="auto", interpolation="none",
    cmap=cmap, norm=norm, extent=[0, len(df), 4.28, 0.28],
)

arch_colors = []
for x in df["proteomics_architecture_class"]:
    if x == "influenza_amplified":
        arch_colors.append(C_FLU)
    elif x == "shared_core":
        arch_colors.append(C_SHARED)
    else:
        arch_colors.append(C_RSV)

for i, color in enumerate(arch_colors):
    axB.add_patch(
        plt.Rectangle((i, 0.02), 1, 0.14,
                      facecolor=color, edgecolor="none", clip_on=False)
    )

axB.set_yticks(np.arange(4) + 0.78)
axB.set_yticklabels([
    "4-season RNA-seq",
    "Proteomics FDR",
    "Regulatory connection",
    "IAV CRISPR dependency",
], fontweight="bold")
axB.set_xticks([])
axB.tick_params(axis="y", width=1.5, length=5.5)
axB.set_xlim(0, len(df))
axB.set_ylim(4.35, -0.10)

classes = df["proteomics_architecture_class"].to_numpy()
for i in range(1, len(classes)):
    if classes[i] != classes[i - 1]:
        axB.axvline(i, color="#666666", lw=1.35)

legend_elements = [
    Line2D([0], [0], marker="s", linestyle="None",
           markerfacecolor="#3A3A3A", markeredgecolor="none",
           markersize=7, label="Supported"),
    Line2D([0], [0], marker="s", linestyle="None",
           markerfacecolor="#FFFFFF", markeredgecolor="#AAAAAA",
           markersize=7, label="Evaluable, not supported"),
    Line2D([0], [0], marker="s", linestyle="None",
           markerfacecolor="#F2F2F2", markeredgecolor="#DDDDDD",
           markersize=7, label="Not evaluable"),
]
legB = axB.legend(
    handles=legend_elements,
    frameon=True, fontsize=10.4,
    loc="lower center", bbox_to_anchor=(0.5, -0.24),
    ncol=3, facecolor="white", edgecolor="0.80",
    framealpha=1.0, borderpad=0.38,
    columnspacing=1.0, handletextpad=0.42,
)
for t in legB.get_texts():
    t.set_fontweight("bold")

# ============================================================
# PANEL C — SESSION 38A EVIDENCE-WEIGHTED PRIORITIZATION
# Ranked component-evidence profile; frozen session38_rank preserved.
# ============================================================

plotC = top30.sort_values("session38_rank").head(15).copy()
plotC = plotC.sort_values("session38_rank", ascending=True).reset_index(drop=True)

component_specs = [
    ("RNA", "score_rna", 4.0, C_RNA),
    ("Protein", "score_proteomics", 3.0, C_PROT),
    ("CRISPR", "score_crispr", 3.0, C_CRISPR),
    ("Regulatory", "score_regulatory", 2.0, C_REG),
    ("Leading edge", "score_leading_edge", 2.0, "#6C8EAD"),
    ("Architecture", "score_architecture", 1.0, C_FLU),
]

yC = np.arange(len(plotC))
xC = np.arange(len(component_specs))

# Light guide grid.
for x in xC:
    axC.axvline(x, color="#EEEEEE", lw=0.9, zorder=0)
for y in yC:
    axC.axhline(y, color="#F4F4F4", lw=0.7, zorder=0)

for yi, (_, r) in enumerate(plotC.iterrows()):
    for xi, (lab, col, vmax, color) in enumerate(component_specs):
        val = r[col]
        if pd.isna(val):
            axC.scatter(
                xi, yi, s=72, marker="s",
                facecolor="#F2F2F2", edgecolor="#D0D0D0",
                linewidth=0.9, zorder=3,
            )
            continue

        frac = max(0.0, min(float(val) / vmax, 1.0))
        if frac == 0:
            axC.scatter(
                xi, yi, s=52, marker="o",
                facecolor="white", edgecolor="#C8C8C8",
                linewidth=0.9, zorder=3,
            )
        else:
            axC.scatter(
                xi, yi,
                s=58 + 150 * frac,
                marker="o",
                facecolor=color, edgecolor="white",
                linewidth=0.8, alpha=0.35 + 0.65 * frac,
                zorder=4,
            )

    score = float(r["evidence_score_coverage_adjusted"])
    axC.text(
        6.62, yi, f"{score:.1f}",
        ha="center", va="center",
        fontsize=18.5, fontweight="bold", color="#444444",
    )

    if r["gene_symbol"] in priority5:
        axC.scatter(
            7.42, yi, s=55, marker="o",
            facecolor="white", edgecolor="#111111",
            linewidth=1.5, zorder=5,
        )
    if r["gene_symbol"] in bridge7:
        axC.scatter(
            8.02, yi, s=60, marker="D",
            facecolor=C_CRISPR, edgecolor="white",
            linewidth=0.8, zorder=5,
        )

labelsC = [
    f'#{int(r.session38_rank)}  {r.gene_symbol}'
    for _, r in plotC.iterrows()
]
axC.set_yticks(yC)
axC.set_yticklabels(labelsC, fontweight="bold", fontsize=19.5)

axC.set_xticks(list(xC) + [6.62, 7.42, 8.02])
axC.set_xticklabels(
    [x[0] for x in component_specs] +
    ["Adjusted\nscore", "Priority-5", "Bridge"],
    rotation=18, ha="right", fontweight="bold", fontsize=18.0,
)

axC.set_xlim(-0.55, 8.45)
axC.set_ylim(-0.5, len(plotC) - 0.5)
axC.invert_yaxis()
axC.tick_params(axis="x", length=0)
axC.tick_params(axis="y", length=0)
axC.spines[["top", "right", "bottom"]].set_visible(False)
axC.spines["left"].set_visible(False)

axC.text(
    0.0, 1.012,
    "Marker size/intensity reflects the frozen component score; "
    "rank is not determined by adjusted score alone.",
    transform=axC.transAxes,
    ha="left", va="bottom",
    fontsize=17.5, fontweight="bold", color="#555555",
)

handlesC = [
    Line2D([0], [0], marker="o", linestyle="None", markersize=6,
           markerfacecolor="white", markeredgecolor="#111111",
           label="Weight-robust priority-5"),
    Line2D([0], [0], marker="D", linestyle="None", markersize=6,
           markerfacecolor=C_CRISPR, markeredgecolor="white",
           label="Response–dependency bridge"),
]
legC = axC.legend(
    handles=handlesC, frameon=False,
    loc="upper right", bbox_to_anchor=(1.0, -0.14),
    fontsize=9.4, ncol=2,
)
for t in legC.get_texts():
    t.set_fontweight("bold")

# ============================================================
# PANEL D — SESSION 38B MOLECULAR RESPONSE VS DEPENDENCY
# ============================================================

bridge_order = [
    "KPNB1", "CHMP5", "TOP2A", "FCGR1B",
    "HIST2H2AC", "OTOF", "HERC5",
]
plotD = bridge.set_index("gene_symbol").loc[bridge_order].reset_index()
yD = np.arange(len(plotD))

axD.set_xlim(-0.55, 4.25)
axD.set_ylim(len(plotD) - 0.5, -0.5)

cols = ["Discovery RNA", "Independent RNA", "Protein", "CRISPR"]
for yi, (_, r) in enumerate(plotD.iterrows()):
    # Discovery RNA
    axD.scatter(0, yi, s=112, facecolor=C_RNA, edgecolor="white", linewidth=0.8)
    # Independent RNA
    axD.scatter(1, yi, s=112, facecolor=C_RNA, edgecolor="white", linewidth=0.8)

    # Protein
    if not bool(r["proteomics_evaluable"]):
        pfc, pec = LIGHT, MID
    elif str(r["proteomics_validation_category"]).startswith("PROTEIN_CONCORDANT"):
        pfc, pec = C_PROT, C_PROT
    else:
        pfc, pec = "white", C_PROT
    axD.scatter(2, yi, s=112, facecolor=pfc, edgecolor=pec, linewidth=1.5)

    # CRISPR
    axD.scatter(3, yi, s=112, facecolor=C_CRISPR, edgecolor="white", linewidth=0.8)
    if bool(r["crispr_iav_dependency_high_confidence"]):
        axD.scatter(
            3, yi, s=165, facecolor="none",
            edgecolor="#8B1A1A", linewidth=2.0, zorder=5,
        )

    # concise molecular context
    prot_txt = "NA"
    if bool(r["proteomics_evaluable"]):
        prot_txt = (
            "concordant" if str(r["proteomics_validation_category"]).startswith("PROTEIN_CONCORDANT")
            else "discordant"
        )
    axD.text(
        3.42, yi,
        prot_txt,
        va="center", ha="left",
        fontsize=18.0, fontweight="bold", color="#555555",
    )

axD.set_xticks(range(4))
axD.set_xticklabels(cols, rotation=5, ha="center", fontweight="bold", fontsize=19.0)
axD.set_yticks(yD)
axD.set_yticklabels(plotD["gene_symbol"], fontweight="bold", fontsize=19.5)
axD.grid(axis="x", color="#EEEEEE", lw=1.0)
axD.spines[["top", "right"]].set_visible(False)

def get_stat(name, test=None):
    z = b_stats[b_stats["analysis"] == name]
    if test is not None:
        z = z[z["test"] == test]
    if len(z) != 1:
        return None
    return z.iloc[0]

r_disc = get_stat("Discovery_RNA_vs_protein", "Pearson")
r_ind = get_stat("Independent_RNA_vs_protein", "Pearson")
r_cr = get_stat("Abs_discovery_RNA_vs_CRISPR_dependency", "Spearman")

stat_text = (
    f'RNA–protein: r={r_disc["effect"]:.2f} (discovery), '
    f'r={r_ind["effect"]:.2f} (independent)\n'
    f'RNA–CRISPR: ρ={r_cr["effect"]:.2f}, FDR={r_cr["FDR"]:.2f}'
)
axD.text(
    0.99, 0.985, stat_text,
    transform=axD.transAxes,
    ha="right", va="top",
    fontsize=12.5, fontweight="bold", color="#555555",
    bbox=dict(
        boxstyle="round,pad=0.28",
        facecolor="white", edgecolor="#D7D7D7",
        linewidth=0.8, alpha=0.95
    ),
    zorder=7,
)

handlesD = [
    Line2D([0], [0], marker="o", linestyle="None",
           markerfacecolor=C_RNA, markeredgecolor="white",
           markersize=7, label="RNA response"),
    Line2D([0], [0], marker="o", linestyle="None",
           markerfacecolor=C_PROT, markeredgecolor=C_PROT,
           markersize=7, label="Protein concordant"),
    Line2D([0], [0], marker="o", linestyle="None",
           markerfacecolor="white", markeredgecolor=C_PROT,
           markersize=7, label="Protein discordant"),
    Line2D([0], [0], marker="o", linestyle="None",
           markerfacecolor=LIGHT, markeredgecolor=MID,
           markersize=7, label="Protein not evaluable"),
]
legD = axD.legend(
    handles=handlesD, frameon=False,
    loc="upper center", bbox_to_anchor=(0.50, -0.16),
    fontsize=9.6, ncol=2,
)
for t in legD.get_texts():
    t.set_fontweight("bold")

# ============================================================
# PANEL E — SESSION 38C CELL-STATE LOCALIZATION
# ============================================================

E = loc_primary.merge(
    loc_robust[
        ["cell_type", "robust_spearman_rho", "robust_spearman_FDR"]
    ],
    on="cell_type",
    how="left",
    validate="one_to_one",
)

E = E.sort_values("spearman_rho", ascending=True).reset_index(drop=True)
yE = np.arange(len(E))

short_names = {
    "natural killer cell": "Natural killer",
    "CD8-positive, alpha-beta T cell": "CD8 T",
    "effector CD8-positive, alpha-beta T cell": "Effector CD8 T",
    "effector CD4-positive, alpha-beta T cell": "Effector CD4 T",
    "classical monocyte": "Classical monocyte",
}
labelsE = [short_names.get(x, x) for x in E["cell_type"]]

for yi, (_, r) in enumerate(E.iterrows()):
    axE.plot(
        [r["spearman_rho"], r["robust_spearman_rho"]],
        [yi, yi],
        color="#A0A0A0", lw=2.0, zorder=1,
    )
    sig = float(r["spearman_FDR"]) < 0.05
    axE.scatter(
        r["spearman_rho"], yi,
        s=105, marker="o",
        facecolor=C_DISC if sig else "white",
        edgecolor=C_DISC, linewidth=1.5, zorder=3,
    )
    axE.scatter(
        r["robust_spearman_rho"], yi,
        s=82, marker="D",
        facecolor="white", edgecolor="#444444",
        linewidth=1.4, zorder=4,
    )
    axE.text(
        0.315, yi,
        f'{100*r["direction_concordance"]:.1f}%',
        ha="left", va="center",
        fontsize=13.0, fontweight="bold",
    )

axE.axvline(0, color="#777777", lw=1.3)
axE.set_yticks(yE)
axE.set_yticklabels(labelsE, fontweight="bold")
axE.set_xlabel("Correlation with frozen bulk influenza effects (Spearman ρ)",
               fontweight="bold")
axE.set_xlim(-0.25, 0.44)
axE.grid(axis="x", color="#EEEEEE", lw=0.9)
axE.set_axisbelow(True)
axE.spines[["top", "right"]].set_visible(False)

axE.text(
    0.315, len(E) - 0.18,
    "Direction\nconcordance",
    ha="left", va="bottom",
    fontsize=12.5, fontweight="bold", color="#555555",
    clip_on=False,
)

handlesE = [
    Line2D([0], [0], marker="o", linestyle="None",
           markerfacecolor=C_DISC, markeredgecolor=C_DISC,
           markersize=7, label="Primary"),
    Line2D([0], [0], marker="D", linestyle="None",
           markerfacecolor="white", markeredgecolor="#444444",
           markersize=7, label="Detection-filtered robustness"),
]
legE = axE.legend(
    handles=handlesE, frameon=False,
    loc="upper center", bbox_to_anchor=(0.46, -0.18),
    fontsize=9.8, ncol=2
)
for t in legE.get_texts():
    t.set_fontweight("bold")

# ============================================================
# PANEL F — SESSION 38D CROSS-SEASON STABILITY
# ============================================================

plotF = arch_stab[
    arch_stab["architecture_class"].isin([
        "Influenza_amplified_shared",
        "Shared_core",
        "RSV_amplified_shared",
    ])
].copy()

label_map = {
    "Influenza_amplified_shared": "Influenza-amplified",
    "Shared_core": "Shared core",
    "RSV_amplified_shared": "RSV-amplified",
}
plotF["label"] = plotF["architecture_class"].map(label_map)
orderF = ["RSV-amplified", "Shared core", "Influenza-amplified"]
plotF["ord"] = plotF["label"].map({x: i for i, x in enumerate(orderF)})
plotF = plotF.sort_values("ord").reset_index(drop=True)

cats = [
    ("Persistent", "n_stable_persistent", "#3B7A57"),
    ("Predominant", "n_stable_predominant", "#8FB996"),
    ("Intermittent", "n_stable_intermittent", "#D9B26F"),
    ("Directionally variable", "n_directionally_variable", "#C75C5C"),
]

yF = np.arange(len(plotF))
left = np.zeros(len(plotF))

for lab, col, color in cats:
    vals = 100 * plotF[col].to_numpy() / plotF["n_evaluable_ge3"].to_numpy()
    axF.barh(
        yF, vals, left=left,
        color=color, edgecolor="white", linewidth=0.7,
        height=0.62, label=lab,
    )
    left += vals

axF.set_yticks(yF)
axF.set_yticklabels(plotF["label"], fontweight="bold")
axF.set_xlabel("Genes evaluable in ≥3 seasons (%)", fontweight="bold")
axF.set_xlim(0, 100)
axF.grid(axis="x", color="#EEEEEE", lw=0.9)
axF.set_axisbelow(True)
axF.spines[["top", "right"]].set_visible(False)

for yi, (_, r) in enumerate(plotF.iterrows()):
    axF.text(
        101.0, yi,
        f'{int(r["n_evaluable_ge3"])}/{int(r["n_genes"])} evaluable',
        ha="left", va="center",
        fontsize=17.0, fontweight="bold", color="#555555",
        clip_on=False,
    )

stable150 = int(
    stability.loc[
        stability["stability_class"].isin(
            ["STABLE_PERSISTENT", "STABLE_PREDOMINANT"]
        ),
        "n_genes",
    ].sum()
)

axF.text(
    0.01, 1.02,
    f"{stable150}/170 genes (88.2%) persistent or predominant across seasons",
    transform=axF.transAxes,
    ha="left", va="bottom",
    fontsize=13.5, fontweight="bold", color="#555555",
)

legF = axF.legend(
    frameon=False, loc="upper center",
    bbox_to_anchor=(0.5, -0.18),
    ncol=2, fontsize=9.8,
)
for t in legF.get_texts():
    t.set_fontweight("bold")

# ============================================================
# PANEL HEADERS — LOCKED SUBMISSION STYLE
# ============================================================

PANEL_LABEL_SIZE = 29.0
PANEL_TITLE_SIZE = 26.0

panel_titles = {
    "A": "Multilayer validation framework",
    "B": "Evidence landscape across the frozen 170-gene architecture",
    "C": "Evidence-weighted influenza host-factor prioritization",
    "D": "Molecular response and functional dependency",
    "E": "Cell-state localization of the influenza program",
    "F": "Stable and context-dependent programs across seasons",
}

fig.canvas.draw()

# Panel B architecture-class labels in figure coordinates.
boxB = axB.get_position()
class_specs = [
    ("influenza_amplified", "Influenza-amplified (n=37)"),
    ("shared_core", "Shared core (n=130)"),
    ("RSV_amplified", "RSV-amplified (n=3)"),
]
label_y = boxB.y1 + 0.005
for cls, label_txt in class_specs:
    idx = np.where(classes == cls)[0]
    if len(idx):
        center_frac = ((idx.min() + idx.max() + 1) / 2) / len(df)
        x_fig = boxB.x0 + center_frac * boxB.width
        fig.text(
            x_fig, label_y, label_txt,
            ha="center", va="bottom",
            fontsize=13.0, fontweight="bold", color="#555555",
        )

for label, axis in zip(
    ["A", "B", "C", "D"],
    [axA, axB, axC, axD],
):
    boxp = axis.get_position()
    ytitle = boxp.y1 + 0.018
    fig.text(
        boxp.x0 - 0.034, ytitle,
        label,
        fontsize=PANEL_LABEL_SIZE,
        fontweight="bold",
        ha="left", va="bottom",
    )
    fig.text(
        boxp.x0, ytitle,
        panel_titles[label],
        fontsize=PANEL_TITLE_SIZE,
        fontweight="bold",
        ha="left", va="bottom",
        linespacing=1.04,
    )

# ============================================================
# FIGURE-FACING PROVENANCE
# ============================================================

prov = pd.DataFrame([
    {
        "file": str(p),
        "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
    }
    for p in inputs
])
prov.to_csv(
    PROVDIR / "Figure7_SUBMISSION_SPLIT_v3_input_sha256.tsv",
    sep="\t", index=False,
)

# Save exact displayed top-15 ordering.
plotC.sort_values("session38_rank").to_csv(
    TABDIR / "Figure7_SUBMISSION_SPLIT_v3_panelC_top15.tsv",
    sep="\t", index=False,
)
plotD.to_csv(
    TABDIR / "Figure7_SUBMISSION_SPLIT_v3_panelD_bridge7.tsv",
    sep="\t", index=False,
)
E.to_csv(
    TABDIR / "Figure7_SUBMISSION_SPLIT_v3_panelE_localization.tsv",
    sep="\t", index=False,
)
plotF.to_csv(
    TABDIR / "Figure7_SUBMISSION_SPLIT_v3_panelF_stability.tsv",
    sep="\t", index=False,
)

# No overall title; no global footer.

# ============================================================
# SESSION 41 v6 — FINAL READABILITY POLISH
# Presentation only; no scientific values or selections changed.
# ============================================================

# ------------------------------------------------------------
# Panel A — enlarge workflow-box/internal text
# ------------------------------------------------------------
for _t in axA.texts:
    _t.set_fontsize(max(float(_t.get_fontsize()), 14.0))
    _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel B — evidence matrix labels + legend
# ------------------------------------------------------------
for _t in axB.get_yticklabels():
    _t.set_fontsize(19.0)
    _t.set_fontweight("bold")

for _t in axB.get_xticklabels():
    _t.set_fontsize(max(float(_t.get_fontsize()), 12.5))
    _t.set_fontweight("bold")

_leg = axB.get_legend()
if _leg is not None:
    for _t in _leg.get_texts():
        _t.set_fontsize(17.5)
        _t.set_fontweight("bold")
    for _h in _leg.legend_handles:
        try:
            _h.set_markersize(9)
        except Exception:
            pass

# ------------------------------------------------------------
# Panel C — prioritization matrix
# Cleaner multiline column labels and stronger text.
# ------------------------------------------------------------
_c_labels = [t.get_text() for t in axC.get_xticklabels()]

_c_map = {
    "Leading edge": "Leading\nedge",
    "Adjusted score": "Adjusted\nscore",
    "Priority-5": "Priority-\n5",
    "Bridge": "Bridge",
}

_c_labels = [_c_map.get(x, x) for x in _c_labels]

axC.set_xticks(axC.get_xticks())
axC.set_xticklabels(
    _c_labels,
    rotation=10,
    ha="center",
    fontsize=12.7,
    fontweight="bold"
)

for _t in axC.get_yticklabels():
    _t.set_fontsize(13.5)
    _t.set_fontweight("bold")

# Increase numerical/annotation text while preserving positions.
for _t in axC.texts:
    _txt = _t.get_text()
    _t.set_fontsize(max(float(_t.get_fontsize()), 12.0))

    # Pull the Direction-concordance annotation slightly upward
    # so it no longer hangs into the C/E inter-panel whitespace.
    if "Direction" in _txt and "concordance" in _txt:
        _x, _y = _t.get_position()
        _t.set_position((_x, _y + 0.10))
        _t.set_fontsize(16.0)
        _t.set_fontweight("bold")

_leg = axC.get_legend()
if _leg is not None:
    for _t in _leg.get_texts():
        _t.set_fontsize(16.0)
        _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel D — bridge/convergence labels, statistics, legend
# ------------------------------------------------------------
for _t in axD.get_yticklabels():
    _t.set_fontsize(17.5)
    _t.set_fontweight("bold")

for _t in axD.get_xticklabels():
    _t.set_fontsize(19.0)
    _t.set_fontweight("bold")

for _t in axD.texts:
    _t.set_fontsize(max(float(_t.get_fontsize()), 12.3))

_leg = axD.get_legend()
if _leg is not None:
    for _t in _leg.get_texts():
        _t.set_fontsize(12.7)
        _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel E — cell-state localization
# ------------------------------------------------------------
for _t in axE.get_yticklabels():
    _t.set_fontsize(13.5)
    _t.set_fontweight("bold")

for _t in axE.get_xticklabels():
    _t.set_fontsize(16.5)
    _t.set_fontweight("bold")

axE.xaxis.label.set_fontsize(20.0)
axE.xaxis.label.set_fontweight("bold")

for _t in axE.texts:
    _t.set_fontsize(max(float(_t.get_fontsize()), 12.5))
    _t.set_fontweight("bold")

_leg = axE.get_legend()
if _leg is not None:
    for _t in _leg.get_texts():
        _t.set_fontsize(16.5)
        _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel F — cross-season stability
# ------------------------------------------------------------
for _t in axF.get_yticklabels():
    _t.set_fontsize(13.5)
    _t.set_fontweight("bold")

for _t in axF.get_xticklabels():
    _t.set_fontsize(16.5)
    _t.set_fontweight("bold")

axF.xaxis.label.set_fontsize(20.0)
axF.xaxis.label.set_fontweight("bold")

for _t in axF.texts:
    _t.set_fontsize(max(float(_t.get_fontsize()), 12.5))
    _t.set_fontweight("bold")

_leg = axF.get_legend()
if _leg is not None:
    for _t in _leg.get_texts():
        _t.set_fontsize(16.5)
        _t.set_fontweight("bold")



# ============================================================
# SESSION 41 v7 — FINAL PANEL A/C CORRECTION
# Scientific content unchanged.
# ============================================================

# ------------------------------------------------------------
# PANEL A
# Slightly reduce workflow-box text from the v6 forced minimum.
# This prevents text/connector crowding while remaining readable.
# ------------------------------------------------------------
for _t in axA.texts:
    _txt = _t.get_text()

    if (
        "4-season" in _txt
        or "SomaScan" in _txt
        or "Regulatory" in _txt
        or "IAV CRISPR" in _txt
        or "Discovery architecture" in _txt
    ):
        _t.set_fontsize(16.5)
        _t.set_fontweight("bold")


# ------------------------------------------------------------
# PANEL C
# Restore a clean single-line categorical axis.
# The previous multiline labels created excessive congestion.
# ------------------------------------------------------------
_c_labels = [
    "RNA",
    "Protein",
    "CRISPR",
    "Regulatory",
    "Leading edge",
    "Architecture",
    "Adjusted score",
    "Priority-5",
    "Bridge",
]

_c_ticks = axC.get_xticks()

# Only apply if the axis has the expected nine columns.
if len(_c_ticks) == len(_c_labels):
    axC.set_xticks(_c_ticks)
    axC.set_xticklabels(
        _c_labels,
        rotation=24,
        ha="right",
        rotation_mode="anchor",
        fontsize=17.0,
        fontweight="bold"
    )

# Preserve strong gene labels.
for _t in axC.get_yticklabels():
    _t.set_fontsize(19.0)
    _t.set_fontweight("bold")

# Keep numerical score text readable but not oversized.
for _t in axC.texts:
    _txt = _t.get_text()

    if _txt.strip().replace(".", "", 1).isdigit():
        _t.set_fontsize(17.0)

    # Move the Direction-concordance annotation completely
    # upward into Panel C rather than the C/E inter-panel gap.
    if "Direction" in _txt and "concordance" in _txt:
        _x, _y = _t.get_position()
        _t.set_position((_x, _y + 0.24))
        _t.set_fontsize(11.2)
        _t.set_fontweight("bold")

# Slightly reduce Panel C legend typography to avoid competing
# with Panel E below.
_leg = axC.get_legend()
if _leg is not None:
    for _t in _leg.get_texts():
        _t.set_fontsize(16.5)
        _t.set_fontweight("bold")


# ============================================================
# SESSION 41 v8 — FINAL SURGICAL LAYOUT PATCH
# Scientific content unchanged.
# Only Panel A connector/text collision and Panel C
# Direction-concordance annotation are modified.
# ============================================================

# ------------------------------------------------------------
# PANEL A
# Remove the visible collision between the vertical workflow
# connector and the first line of the green 4-season RNA-seq box.
#
# Keep the workflow boxes and their scientific content unchanged.
# Shift the green-box text slightly downward inside its box.
# ------------------------------------------------------------
for _t in axA.texts:
    _txt = _t.get_text()

    if "4-season" in _txt and "RNA-seq" in _txt:
        _x, _y = _t.get_position()
        _t.set_position((_x, _y - 0.018))
        _t.set_fontsize(16.0)
        _t.set_fontweight("bold")


# ------------------------------------------------------------
# PANEL C
# Pull "Direction concordance" completely back into Panel C.
# Do not allow it to occupy the C/E inter-panel whitespace.
# ------------------------------------------------------------
for _t in axC.texts:
    _txt = _t.get_text()

    if "Direction" in _txt and "concordance" in _txt:
        # Use axes coordinates so placement is deterministic.
        _t.set_transform(axC.transAxes)
        _t.set_position((0.86, -0.105))
        _t.set_ha("center")
        _t.set_va("top")
        _t.set_fontsize(10.5)
        _t.set_fontweight("bold")
        _t.set_clip_on(False)


# ============================================================
# SESSION 41 v9 — FINAL PANEL C ANNOTATION CORRECTION
# Scientific content unchanged.
# ============================================================

# Remove the Direction-concordance annotation from the
# C/E inter-panel whitespace and place it inside Panel C.
for _t in axC.texts:
    _txt = _t.get_text()

    if "Direction" in _txt and "concordance" in _txt:
        _t.set_transform(axC.transAxes)

        # Inside lower-right portion of Panel C.
        _t.set_position((0.88, 0.035))

        _t.set_ha("center")
        _t.set_va("bottom")
        _t.set_fontsize(9.8)
        _t.set_fontweight("bold")
        _t.set_clip_on(True)


# ============================================================
# SESSION 41 v10 — FINAL FIGURE 7 PRESENTATION POLISH
#
# Presentation-only changes:
#   1. Remove redundant Direction-concordance annotation in C.
#   2. Remove Panel D grid/guide lines through markers.
#   3. Ensure Panel D points render above all axis elements.
#   4. Rebuild Panel D legend with exact plotted colors/styles.
#
# Scientific values, classifications, genes and statistics
# remain unchanged.
# ============================================================


# ------------------------------------------------------------
# PANEL C
# Remove redundant Direction-concordance annotation.
# ------------------------------------------------------------

for _t in list(axC.texts):
    _txt = _t.get_text().strip().lower()

    if "direction" in _txt and "concordance" in _txt:
        _t.set_visible(False)


# ------------------------------------------------------------
# PANEL D
# Remove vertical guide/grid lines that visually split markers.
# ------------------------------------------------------------

axD.grid(False)
axD.xaxis.grid(False)
axD.yaxis.grid(False)

# Make sure plotted scientific markers sit above axes/background.
axD.set_axisbelow(True)

for _collection in axD.collections:
    try:
        _collection.set_zorder(5)
    except Exception:
        pass


# ------------------------------------------------------------
# PANEL D — rebuild legend so symbols exactly match plot.
# ------------------------------------------------------------

from matplotlib.lines import Line2D

# Remove legacy legend.
_old_leg = axD.get_legend()
if _old_leg is not None:
    _old_leg.remove()

# Colors already defined by the canonical Session 38 script:
# C_RNA    = RNA response
# C_PROT   = proteomics
# C_CRISPR = CRISPR dependency

legend_handles_D = [

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=9,
        markerfacecolor=C_RNA,
        markeredgecolor="white",
        markeredgewidth=0.9,
        label="RNA response"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=9,
        markerfacecolor=C_PROT,
        markeredgecolor="white",
        markeredgewidth=0.9,
        label="Protein concordant"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=9,
        markerfacecolor="white",
        markeredgecolor=C_PROT,
        markeredgewidth=1.4,
        label="Protein discordant"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=9,
        markerfacecolor="#F2F2F2",
        markeredgecolor="#B5B5B5",
        markeredgewidth=1.1,
        label="Protein not evaluable"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=9,
        markerfacecolor=C_CRISPR,
        markeredgecolor="white",
        markeredgewidth=0.9,
        label="CRISPR dependency"
    ),
]

legD_final = axD.legend(
    handles=legend_handles_D,
    loc="upper center",
    bbox_to_anchor=(0.50, -0.16),
    ncol=3,
    frameon=False,
    fontsize=11.5,
    columnspacing=1.35,
    handletextpad=0.45
)

for _txt in legD_final.get_texts():
    _txt.set_fontweight("bold")


# ------------------------------------------------------------
# PANEL D — make circular points visually intact.
# Ensure each scatter collection has a visible clean edge and
# remains above the axis background.
# ------------------------------------------------------------

for _collection in axD.collections:
    try:
        _collection.set_zorder(6)
    except Exception:
        pass



# ============================================================
# SESSION 41 v11 — FINAL PANEL C/D CORRECTION
#
# Scientific interpretation:
#   KPNB1 = RNA + RNA + concordant protein + CRISPR
#   Black CRISPR ring = high-confidence CRISPR dependency
#
# No analytical values or classifications are changed.
# ============================================================

from matplotlib.lines import Line2D


# ------------------------------------------------------------
# PANEL C
# Remove redundant Direction-concordance annotation everywhere.
# ------------------------------------------------------------

for _ax in fig.axes:
    for _t in list(_ax.texts):
        _txt = " ".join(_t.get_text().split()).lower()
        if "direction" in _txt and "concordance" in _txt:
            _t.set_visible(False)

for _t in list(fig.texts):
    _txt = " ".join(_t.get_text().split()).lower()
    if "direction" in _txt and "concordance" in _txt:
        _t.set_visible(False)


# ------------------------------------------------------------
# PANEL D
# Remove grid/guide lines through scientific markers.
# ------------------------------------------------------------

axD.grid(False)
axD.xaxis.grid(False)
axD.yaxis.grid(False)


# ------------------------------------------------------------
# PANEL D
# Move the RNA-protein / RNA-CRISPR statistics annotation
# COMPLETELY ABOVE the data region.
# ------------------------------------------------------------

for _t in axD.texts:
    _txt = _t.get_text()

    if "RNA" in _txt and "protein" in _txt and "CRISPR" in _txt:

        _t.set_transform(axD.transAxes)
        _t.set_position((0.98, 1.015))

        _t.set_ha("right")
        _t.set_va("bottom")
        _t.set_clip_on(False)

        _t.set_fontsize(16.5)
        _t.set_fontweight("bold")

        # Remove the opaque box so no point can be concealed.
        _t.set_bbox(
            dict(
                facecolor="none",
                edgecolor="none",
                pad=0.0
            )
        )


# ------------------------------------------------------------
# PANEL D
# Redraw bridge-gene evidence ABOVE all previous artists.
#
# x = 0  Discovery RNA
# x = 1  Independent RNA
# x = 2  Protein
# x = 3  CRISPR
# ------------------------------------------------------------

# Match the existing Panel D y-axis order exactly.
_y_by_gene = {
    str(gene): yi
    for yi, gene in enumerate(plotD["gene_symbol"])
}

for _, _r in plotD.iterrows():

    _gene = str(_r["gene_symbol"])
    _yi = _y_by_gene[_gene]

    # -------------------------
    # Discovery RNA
    # -------------------------
    axD.scatter(
        0, _yi,
        s=125,
        marker="o",
        facecolor=C_RNA,
        edgecolor="white",
        linewidth=1.0,
        zorder=20
    )

    # -------------------------
    # Independent RNA
    # -------------------------
    axD.scatter(
        1, _yi,
        s=125,
        marker="o",
        facecolor=C_RNA,
        edgecolor="white",
        linewidth=1.0,
        zorder=20
    )

    # -------------------------
    # Protein
    # -------------------------
    _protein_eval = bool(_r["proteomics_evaluable"])
    _protein_cat = str(_r["proteomics_validation_category"])

    if not _protein_eval:

        # Not evaluable
        axD.scatter(
            2, _yi,
            s=125,
            marker="o",
            facecolor="#F2F2F2",
            edgecolor="#AFAFAF",
            linewidth=1.2,
            zorder=20
        )

    elif _protein_cat.startswith("PROTEIN_CONCORDANT"):

        # Concordant protein evidence — filled orange
        axD.scatter(
            2, _yi,
            s=125,
            marker="o",
            facecolor=C_PROT,
            edgecolor="white",
            linewidth=1.0,
            zorder=20
        )

    else:

        # Discordant protein evidence — orange outline
        axD.scatter(
            2, _yi,
            s=125,
            marker="o",
            facecolor="white",
            edgecolor=C_PROT,
            linewidth=1.8,
            zorder=20
        )


    # -------------------------
    # CRISPR
    # -------------------------
    _crispr_hc = bool(
        _r["crispr_iav_dependency_high_confidence"]
    )

    if _crispr_hc:

        # High-confidence dependency:
        # red filled point + BLACK outer border
        axD.scatter(
            3, _yi,
            s=140,
            marker="o",
            facecolor=C_CRISPR,
            edgecolor="black",
            linewidth=1.6,
            zorder=21
        )

    else:

        # Standard CRISPR dependency
        axD.scatter(
            3, _yi,
            s=125,
            marker="o",
            facecolor=C_CRISPR,
            edgecolor="white",
            linewidth=1.0,
            zorder=20
        )


# ------------------------------------------------------------
# PANEL D
# Rebuild explanatory legend.
# ------------------------------------------------------------

_old_leg = axD.get_legend()
if _old_leg is not None:
    _old_leg.remove()

_handles = [

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=8,
        markerfacecolor=C_RNA,
        markeredgecolor="white",
        label="RNA response"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=8,
        markerfacecolor=C_PROT,
        markeredgecolor="white",
        label="Protein concordant"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=8,
        markerfacecolor="white",
        markeredgecolor=C_PROT,
        markeredgewidth=1.5,
        label="Protein discordant"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=8,
        markerfacecolor="#F2F2F2",
        markeredgecolor="#AFAFAF",
        label="Protein not evaluable"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=8,
        markerfacecolor=C_CRISPR,
        markeredgecolor="white",
        label="CRISPR dependency"
    ),

    Line2D(
        [0], [0],
        marker="o",
        linestyle="none",
        markersize=8,
        markerfacecolor=C_CRISPR,
        markeredgecolor="black",
        markeredgewidth=1.4,
        label="High-confidence CRISPR"
    ),
]

_leg_final = axD.legend(
    handles=_handles,
    loc="upper center",
    bbox_to_anchor=(0.50, -0.16),
    ncol=3,
    frameon=False,
    fontsize=16.5,
    columnspacing=1.45,
    handletextpad=0.45
)

for _t in _leg_final.get_texts():
    _t.set_fontweight("bold")



# ============================================================
# SESSION 41 v12 — FINAL PANEL D STATISTICS POSITION
# Presentation only. Scientific content unchanged.
# ============================================================

for _t in axD.texts:
    _txt = _t.get_text()

    if (
        "RNA" in _txt
        and "protein" in _txt
        and "CRISPR" in _txt
    ):
        # Put statistics inside the upper-right plot area,
        # below the panel title and above the KPNB1 row.
        _t.set_transform(axD.transAxes)
        _t.set_position((0.985, 0.965))
        _t.set_ha("right")
        _t.set_va("top")
        _t.set_clip_on(False)

        _t.set_fontsize(9.6)
        _t.set_fontweight("bold")
        _t.set_color("#555555")

        _t.set_bbox(
            dict(
                boxstyle="round,pad=0.18",
                facecolor="white",
                edgecolor="none",
                alpha=0.88
            )
        )


# ============================================================
# SESSION 41 v13 — FINAL PANEL D STATISTICS POLISH
# Presentation only. Scientific content unchanged.
# ============================================================

for _t in axD.texts:
    _txt = _t.get_text()

    if (
        "RNA" in _txt
        and "protein" in _txt
        and "CRISPR" in _txt
    ):
        _t.set_transform(axD.transAxes)

        # Upper-right internal whitespace:
        # below the panel heading but above the KPNB1 row.
        _t.set_position((0.975, 0.915))

        _t.set_ha("right")
        _t.set_va("top")
        _t.set_clip_on(False)

        _t.set_fontsize(16.5)
        _t.set_fontweight("bold")
        _t.set_color("#4A4A4A")

        _t.set_bbox(
            dict(
                boxstyle="round,pad=0.20",
                facecolor="white",
                edgecolor="#D0D0D0",
                linewidth=0.7,
                alpha=0.94
            )
        )


# ============================================================
# SESSION 41 v14 — FINAL PANEL D STATISTICS PLACEMENT
# Presentation only; scientific content unchanged.
# ============================================================

for _t in axD.texts:
    _txt = _t.get_text()

    if (
        "RNA" in _txt
        and "protein" in _txt
        and "CRISPR" in _txt
    ):
        # Place annotation entirely ABOVE the plotting area.
        _t.set_transform(axD.transAxes)
        _t.set_position((0.98, 1.055))

        _t.set_ha("right")
        _t.set_va("bottom")
        _t.set_clip_on(False)

        _t.set_fontsize(16.0)
        _t.set_fontweight("bold")
        _t.set_color("#555555")

        # No box: keep the title/annotation area visually clean.
        _t.set_bbox(
            dict(
                facecolor="none",
                edgecolor="none",
                pad=0.0
            )
        )


# ============================================================
# SESSION 42 v17 — TARGETED AXIS TYPOGRAPHY + PANEL C CLEANUP
# Presentation only; scientific content unchanged.
# ============================================================

# ------------------------------------------------------------
# Panel C — larger x-axis typography.
# Keep scientific columns unchanged, but make labels readable.
# ------------------------------------------------------------
_c_ticklabels = list(axC.get_xticklabels())
_c_text = [x.get_text() for x in _c_ticklabels]
_c_pos = axC.get_xticks()

# First six evidence columns: larger, bold, modest rotation.
for i, _t in enumerate(_c_ticklabels):
    if i < 6:
        _t.set_fontsize(19.5)
        _t.set_fontweight("bold")
        _t.set_rotation(18)
        _t.set_ha("right")
        _t.set_rotation_mode("anchor")

# Hide the final three crowded bottom labels and redraw them as
# horizontal headers above their respective columns.
# Cache the Text objects first; repeatedly calling get_xticklabels()
# after hiding labels can return a shorter list in Matplotlib.
for i in range(6, min(9, len(_c_ticklabels))):
    _c_ticklabels[i].set_visible(False)

_header_names = ["Adjusted score", "Priority-5", "Bridge"]
for _x, _lab in zip(_c_pos[6:9], _header_names):
    axC.text(
        _x, -0.72, _lab,
        ha="center", va="bottom",
        fontsize=19.0,
        fontweight="bold",
        clip_on=False,
    )

# Larger Panel C gene names.
for _t in axC.get_yticklabels():
    _t.set_fontsize(19.5)
    _t.set_fontweight("bold")

# Move Panel C legend lower to separate it from the enlarged
# column labels/headers.
_legC = axC.get_legend()
if _legC is not None:
    _legC.set_bbox_to_anchor((0.70, -0.22))
    try:
        _legC.set_loc("upper center")
    except Exception:
        pass
    for _t in _legC.get_texts():
        _t.set_fontsize(15.5)
        _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel D — larger x and y axis typography.
# ------------------------------------------------------------
for _t in axD.get_xticklabels():
    _t.set_fontsize(21.0)
    _t.set_fontweight("bold")

for _t in axD.get_yticklabels():
    _t.set_fontsize(21.0)
    _t.set_fontweight("bold")

# Give Panel D legend extra clearance from the x-axis.
_legD = axD.get_legend()
if _legD is not None:
    _legD.set_bbox_to_anchor((0.50, -0.22))
    try:
        _legD.set_loc("upper center")
    except Exception:
        pass
    for _t in _legD.get_texts():
        _t.set_fontsize(16.0)
        _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel E — larger y-axis cell-state labels.
# ------------------------------------------------------------
for _t in axE.get_yticklabels():
    _t.set_fontsize(21.0)
    _t.set_fontweight("bold")

# Slightly larger x ticks to balance the larger y labels.
for _t in axE.get_xticklabels():
    _t.set_fontsize(18.0)
    _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel F — larger y-axis architecture labels.
# ------------------------------------------------------------
for _t in axF.get_yticklabels():
    _t.set_fontsize(21.0)
    _t.set_fontweight("bold")

for _t in axF.get_xticklabels():
    _t.set_fontsize(18.0)
    _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel A — modest final increase in workflow-box text.
# ------------------------------------------------------------
for _t in axA.texts:
    _t.set_fontsize(max(float(_t.get_fontsize()), 17.0))
    _t.set_fontweight("bold")

# ------------------------------------------------------------
# Panel B — slightly larger architecture class headings and legend.
# ------------------------------------------------------------
for _t in axB.get_yticklabels():
    _t.set_fontsize(19.0)
    _t.set_fontweight("bold")

_legB = axB.get_legend()
if _legB is not None:
    for _t in _legB.get_texts():
        _t.set_fontsize(16.5)
        _t.set_fontweight("bold")


# ============================================================
# SESSION 42 v18 — FINAL NARROW POLISH
# Presentation only; scientific content unchanged.
# Lock D/E/F typography and all frozen science exactly as v17.
# ============================================================

# ------------------------------------------------------------
# Panel A — enlarged workflow footprint with fitted internal text.
# Presentation only; values and labels are unchanged.
# ------------------------------------------------------------
for _t in axA.texts:
    _txt = _t.get_text().strip()
    if not _txt:
        continue
    if _txt.startswith("Discovery architecture"):
        _t.set_fontsize(19.5)
    else:
        _t.set_fontsize(17.5)
    _t.set_fontweight("bold")
    try:
        _t.set_linespacing(1.08)
    except Exception:
        pass

# ------------------------------------------------------------
# Panel B — slightly enlarge secondary typography.
# Keep the evidence matrix and science untouched.
# ------------------------------------------------------------
# Row labels
for _t in axB.get_yticklabels():
    _t.set_fontsize(max(float(_t.get_fontsize()), 19.5))
    _t.set_fontweight("bold")

# Architecture-class headings and other free text above matrix.
for _t in axB.texts:
    _txt = _t.get_text().strip()
    if _txt:
        _t.set_fontsize(max(float(_t.get_fontsize()), 15.5))
        _t.set_fontweight("bold")

# Legend
_legB = axB.get_legend()
if _legB is not None:
    for _t in _legB.get_texts():
        _t.set_fontsize(17.0)
        _t.set_fontweight("bold")
    try:
        _legB.set_bbox_to_anchor((0.5, -0.10))
    except Exception:
        pass

# ------------------------------------------------------------
# Panel C — retain v17 axis typography; only add a little more
# vertical clearance below the matrix for the symbol legend.
# ------------------------------------------------------------
_legC = axC.get_legend()
if _legC is not None:
    _legC.set_bbox_to_anchor((0.70, -0.26))
    try:
        _legC.set_loc("upper center")
    except Exception:
        pass
    for _t in _legC.get_texts():
        _t.set_fontsize(max(float(_t.get_fontsize()), 15.5))
        _t.set_fontweight("bold")

# Keep v17 D/E/F axes unchanged by design.


# ============================================================
# SPLIT FIGURE FINALIZATION
# Panels E/F are intentionally moved to the new Figure 8.
# Keep their frozen computations available in-script but exclude
# all of their artists from this four-panel Figure 7.
# ============================================================
axE.set_visible(False)
axF.set_visible(False)

fig.savefig(PDF, bbox_inches="tight")
fig.savefig(PNG, dpi=600, bbox_inches="tight")
fig.savefig(SVG, format="svg", bbox_inches="tight")
plt.close(fig)

print("=== FIGURE 7 SUBMISSION SPLIT v3 PANEL A BOX ENLARGEMENT COMPLETE ===")
print("PDF:", PDF)
print("PNG:", PNG)
print("SVG:", SVG)
print()
print("Panel C genes:")
print(
    plotC.sort_values("session38_rank", ascending=True)[
        ["session38_rank", "gene_symbol",
         "evidence_score_coverage_adjusted",
         "response_dependency_bridge"]
    ].to_string(index=False)
)
print()
print("Panel D bridge genes:")
print(
    plotD[
        ["gene_symbol", "session38B_convergence_class",
         "crispr_iav_dependency_high_confidence"]
    ].to_string(index=False)
)
print()
print("Panel E localization:")
print(
    E[
        ["cell_type", "spearman_rho", "spearman_FDR",
         "robust_spearman_rho"]
    ].to_string(index=False)
)
print()
print("STATUS: FIGURE 7 SUBMISSION SPLIT v3 PANEL A BOX ENLARGEMENT CANDIDATE")
