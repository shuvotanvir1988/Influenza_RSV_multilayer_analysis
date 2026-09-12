#!/usr/bin/env python3

from pathlib import Path
import hashlib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

# ============================================================
# FIGURE 1 SUBMISSION v15
# FINAL SCRIPT-GENERATED DESIGN CANDIDATE
# ============================================================

ROOT = Path.home() / "Influenza_RSV_Project"

OUTDIR = ROOT / "results/session41_figure_editing/Figure1_SUBMISSION_v15"
FIGDIR = OUTDIR / "figures"
TABDIR = OUTDIR / "tables"
PROVDIR = OUTDIR / "provenance"

for d in (FIGDIR, TABDIR, PROVDIR):
    d.mkdir(parents=True, exist_ok=True)

PDF = FIGDIR / "Figure1_SUBMISSION_v15_study_design.pdf"
PNG = FIGDIR / "Figure1_SUBMISSION_v15_study_design.png"
SVG = FIGDIR / "Figure1_SUBMISSION_v15_study_design.svg"

# ============================================================
# VERIFIED MANUSCRIPT VALUES
# ============================================================

V = {
    "GPL10558_control": 8,
    "GPL10558_RSV": 28,

    "GPL6884_control": 31,
    "GPL6884_influenza": 16,
    "GPL6884_RSV": 107,

    "genes_tested": 25440,
    "influenza_FDR05": 3501,

    "frozen_total": 170,
    "shared_core": 130,
    "influenza_amplified": 37,
    "RSV_amplified": 3,

    "fourseason_evaluable": 157,
    "fourseason_concordant": 154,

    "proteomics_validation": 84,
    "singlecell_types": 5,
}

# ============================================================
# STYLE
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 14,
    "font.weight": "bold",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

C = {
    "navy": "#123A60",
    "blue": "#246FA5",
    "teal": "#25857F",
    "orange": "#C77800",
    "red": "#B83E3E",

    "ink": "#20262D",
    "gray": "#64717C",
    "line": "#C5D1DA",
    "white": "#FFFFFF",

    "pale_blue": "#EAF3F8",
    "pale_teal": "#EAF5F3",
    "pale_orange": "#FBF3E4",
    "pale_red": "#F9ECEC",
    "pale_gray": "#F3F5F7",
}

# Strong dark arrow color requested
ARROW = "#134D85"

fig = plt.figure(figsize=(18.5, 11.0), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")


def rounded(x, y, w, h, fc, ec, lw=1.6, radius=0.012, z=2):
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.007,rounding_size={radius}",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        zorder=z,
    )
    ax.add_patch(p)
    return p


def txt(x, y, s, fs=13, color=None, ha="left",
        va="center", weight="bold", linespacing=1.08, z=5):
    ax.text(
        x, y, s,
        fontsize=fs,
        color=color or C["ink"],
        ha=ha,
        va=va,
        fontweight=weight,
        family="DejaVu Sans",
        linespacing=linespacing,
        zorder=z,
    )


def arrow(x1, y1, x2, y2, color=ARROW, lw=2.0, ms=13, z=4):
    p = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>",
        mutation_scale=ms,
        linewidth=lw,
        color=color,
        shrinkA=1,
        shrinkB=1,
        zorder=z,
    )
    ax.add_patch(p)
    return p


def panel_header(letter, title, x, y):
    txt(x, y, letter, fs=25, va="top")
    txt(x + 0.031, y, title, fs=21.5, va="top")


# ============================================================
# PANEL A
# ============================================================

panel_header(
    "A",
    "Discovery cohort and construction of the frozen architecture",
    0.020, 0.975
)

rounded(
    0.020, 0.390, 0.610, 0.525,
    C["white"], C["line"], lw=1.5
)

txt(0.080, 0.885, "Discovery cohort",
    fs=17.5, color=C["navy"])

txt(0.350, 0.885, "Influenza-centered derivation",
    fs=17.5, color=C["navy"])

# ------------------------------------------------------------
# Discovery cohort
# ------------------------------------------------------------

rounded(
    0.040, 0.500, 0.250, 0.340,
    C["pale_blue"], C["blue"], lw=1.6
)

txt(
    0.065, 0.805,
    "GSE38900 whole-blood microarray",
    fs=15.3, color=C["navy"]
)

platform_x = 0.065
group_x = 0.142
n_x = 0.278

# GPL10558
txt(platform_x, 0.750, "GPL10558",
    fs=15, color=C["blue"])

txt(group_x, 0.750, "Healthy control", fs=13.2)
txt(n_x, 0.750, f"n={V['GPL10558_control']}",
    fs=13.2, ha="right")

txt(group_x, 0.705, "RSV acute",
    fs=13.2, color=C["gray"])
txt(n_x, 0.705, f"n={V['GPL10558_RSV']}",
    fs=13.2, color=C["gray"], ha="right")

ax.plot(
    [0.060, 0.270],
    [0.665, 0.665],
    color=C["blue"],
    linewidth=1.1,
    linestyle=":",
    zorder=1
)

# GPL6884
txt(platform_x, 0.620, "GPL6884",
    fs=15, color=C["blue"])

txt(group_x, 0.620, "Healthy control", fs=13.2)
txt(n_x, 0.620, f"n={V['GPL6884_control']}",
    fs=13.2, ha="right")

txt(group_x, 0.575, "Influenza A acute",
    fs=13.1, color=C["red"])
txt(n_x, 0.575, f"n={V['GPL6884_influenza']}",
    fs=13.1, color=C["red"], ha="right")

txt(group_x, 0.530, "RSV acute",
    fs=13.2, color=C["gray"])
txt(n_x, 0.530, f"n={V['GPL6884_RSV']}",
    fs=13.2, color=C["gray"], ha="right")

# Primary contrast BELOW cohort card
rounded(
    0.050, 0.420, 0.230, 0.060,
    C["pale_red"], C["red"], lw=1.4
)

txt(
    0.165, 0.458,
    "Primary contrast",
    fs=12.2, color=C["red"], ha="center"
)

txt(
    0.165, 0.433,
    "Influenza A acute vs healthy control",
    fs=10.8, ha="center"
)

# ------------------------------------------------------------
# Frozen architecture derivation
# ------------------------------------------------------------

FLOW_X = 0.315
FLOW_W = 0.280

rounded(
    FLOW_X, 0.765, FLOW_W, 0.060,
    C["pale_blue"], C["blue"], lw=1.5
)
txt(
    FLOW_X + FLOW_W/2, 0.795,
    f"{V['genes_tested']:,} genes tested",
    fs=15.5, color=C["blue"], ha="center"
)

rounded(
    FLOW_X, 0.680, FLOW_W, 0.060,
    C["pale_orange"], C["orange"], lw=1.5
)
txt(
    FLOW_X + FLOW_W/2, 0.710,
    f"{V['influenza_FDR05']:,} genes at FDR < 0.05",
    fs=15.5, color=C["orange"], ha="center"
)

arrow(
    FLOW_X + FLOW_W/2, 0.765,
    FLOW_X + FLOW_W/2, 0.742
)

rounded(
    FLOW_X, 0.595, FLOW_W, 0.060,
    C["pale_teal"], C["teal"], lw=1.5
)
txt(
    FLOW_X + FLOW_W/2, 0.625,
    "Cross-pathogen / replication framework",
    fs=14.0, color=C["teal"], ha="center"
)

arrow(
    FLOW_X + FLOW_W/2, 0.680,
    FLOW_X + FLOW_W/2, 0.657
)

rounded(
    0.335, 0.500, 0.240, 0.068,
    C["navy"], C["navy"], lw=1.7
)

txt(
    0.455, 0.540,
    "FROZEN 170-GENE ARCHITECTURE",
    fs=15.2, color=C["white"], ha="center"
)

txt(
    0.455, 0.516,
    "shared-response framework",
    fs=10.6, color=C["white"], ha="center"
)

arrow(
    FLOW_X + FLOW_W/2, 0.595,
    FLOW_X + FLOW_W/2, 0.570
)

# 130 / 37 / 3
class_y = 0.408
class_w = 0.085
class_h = 0.060

class_x = [0.315, 0.410, 0.505]

class_data = [
    (130, "Shared core", C["pale_blue"], C["blue"]),
    (37, "Influenza-\namplified", C["pale_red"], C["red"]),
    (3, "RSV-\namplified", C["pale_gray"], C["gray"]),
]

for x, (n, lab, fc, ec) in zip(class_x, class_data):
    rounded(x, class_y, class_w, class_h, fc, ec, lw=1.4)

    txt(
        x + class_w/2,
        class_y + 0.041,
        str(n),
        fs=17, ha="center"
    )

    txt(
        x + class_w/2,
        class_y + 0.014,
        lab,
        fs=9.4, ha="center",
        linespacing=0.95
    )

# Strong visible arrows
for x in class_x:
    arrow(
        0.455, 0.500,
        x + class_w/2,
        class_y + class_h,
        lw=1.8, ms=11
    )


# ============================================================
# PANEL B
# ============================================================

panel_header(
    "B",
    "Multilayer biological characterization",
    0.020, 0.360
)

rounded(
    0.020, 0.050, 0.610, 0.260,
    C["white"], C["line"], lw=1.5
)

module_y = 0.190
module_h = 0.070
module_w = 0.130
module_x = [0.040, 0.185, 0.330, 0.475]

mods = [
    ("Pathway\narchitecture", C["pale_teal"], C["teal"]),
    ("Regulatory /\nTF activity", C["pale_orange"], C["orange"]),
    ("Cellular\ncontext", C["pale_blue"], C["blue"]),
    ("Gene / leading-edge\narchitecture", C["pale_red"], C["red"]),
]

for x, (lab, fc, ec) in zip(module_x, mods):
    rounded(
        x, module_y, module_w, module_h,
        fc, ec, lw=1.4
    )

    txt(
        x + module_w/2,
        module_y + module_h/2,
        lab,
        fs=11.3, ha="center"
    )

summary_x = 0.140
summary_y = 0.075
summary_w = 0.370
summary_h = 0.075

rounded(
    summary_x, summary_y,
    summary_w, summary_h,
    C["navy"], C["navy"], lw=1.7
)

txt(
    summary_x + summary_w/2,
    summary_y + 0.050,
    "INTERFERON-CENTERED INFLUENZA",
    fs=14.2, color=C["white"], ha="center"
)

txt(
    summary_x + summary_w/2,
    summary_y + 0.023,
    "HOST-RESPONSE ARCHITECTURE",
    fs=14.2, color=C["white"], ha="center"
)

# Publication-style convergence:
# each module points downward to a common horizontal rail,
# followed by one central arrow into the final architecture box.

center_x = summary_x + summary_w / 2

rail_y = 0.166

# Horizontal collection rail.
left_rail = module_x[0] + module_w / 2
right_rail = module_x[-1] + module_w / 2

ax.plot(
    [left_rail, right_rail],
    [rail_y, rail_y],
    color=ARROW,
    linewidth=1.75,
    solid_capstyle="round",
    zorder=2
)

# Four short vertical arrows from modules to rail.
for x in module_x:
    mx = x + module_w / 2

    arrow(
        mx,
        module_y,
        mx,
        rail_y,
        color=ARROW,
        lw=1.65,
        ms=9
    )

# One central arrow from rail to summary box.
arrow(
    center_x,
    rail_y,
    center_x,
    summary_y + summary_h,
    color=ARROW,
    lw=1.9,
    ms=11
)


# ============================================================
# PANEL C
# ============================================================

panel_header(
    "C",
    "Independent validation and orthogonal extension",
    0.655, 0.975
)

rounded(
    0.655, 0.050, 0.320, 0.865,
    C["white"], C["line"], lw=1.5
)

card_x = 0.680
card_w = 0.270
card_h = 0.110

cards = [
    (
        0.750,
        "Four-season RNA-seq",
        f"{V['fourseason_concordant']}/{V['fourseason_evaluable']} direction-concordant\nin all four seasons",
        C["pale_blue"], C["blue"]
    ),
    (
        0.605,
        "SomaScan proteomics",
        f"{V['proteomics_validation']}-gene frozen validation set",
        C["pale_teal"], C["teal"]
    ),
    (
        0.460,
        "IAV CRISPR",
        "Functional host-dependency evidence",
        C["pale_red"], C["red"]
    ),
    (
        0.315,
        "Single-cell localization",
        f"{V['singlecell_types']} prespecified cell states",
        C["pale_orange"], C["orange"]
    ),
]

for y, title, body, fc, ec in cards:
    rounded(
        card_x, y, card_w, card_h,
        fc, ec, lw=1.5
    )

    txt(
        card_x + 0.018,
        y + 0.073,
        title,
        fs=14.2, color=ec
    )

    txt(
        card_x + 0.018,
        y + 0.035,
        body,
        fs=11.2,
        linespacing=1.05
    )

# Integration box
INT_X = 0.700
INT_Y = 0.145
INT_W = 0.230
INT_H = 0.100

rounded(
    INT_X, INT_Y,
    INT_W, INT_H,
    C["navy"], C["navy"], lw=1.7
)

txt(
    INT_X + INT_W/2,
    INT_Y + 0.063,
    "EVIDENCE-WEIGHTED",
    fs=13.7, color=C["white"], ha="center"
)

txt(
    INT_X + INT_W/2,
    INT_Y + 0.038,
    "INTEGRATION",
    fs=13.7, color=C["white"], ha="center"
)

txt(
    INT_X + INT_W/2,
    INT_Y + 0.015,
    "prioritization • convergence • stability",
    fs=9.6, color=C["white"], ha="center"
)

# Minimal publication-style Panel C connector.
# The four evidence cards are intentionally presented as
# independent orthogonal evidence layers.
# No arrows connect the evidence cards to one another.
# One short arrow indicates transition from the evidence stack
# to evidence-weighted integration.

INTEGRATION_CENTER_X = INT_X + INT_W / 2
INTEGRATION_TOP_Y = INT_Y + INT_H

# Bottom of the final evidence card (Single-cell localization).
single_cell_bottom_y = cards[-1][0]

# Leave visible whitespace below the card before the arrow starts.
arrow_start_y = single_cell_bottom_y - 0.022

# Leave visible whitespace above the integration box.
arrow_end_y = INTEGRATION_TOP_Y + 0.012

arrow(
    INTEGRATION_CENTER_X,
    arrow_start_y,
    INTEGRATION_CENTER_X,
    arrow_end_y,
    color=ARROW,
    lw=1.9,
    ms=11
)

# ============================================================
# SOURCE DATA
# ============================================================

source = pd.DataFrame([
    ["GPL10558 healthy control", V["GPL10558_control"]],
    ["GPL10558 RSV acute", V["GPL10558_RSV"]],

    ["GPL6884 healthy control", V["GPL6884_control"]],
    ["GPL6884 influenza A acute", V["GPL6884_influenza"]],
    ["GPL6884 RSV acute", V["GPL6884_RSV"]],

    ["Genes tested", V["genes_tested"]],
    ["Influenza FDR < 0.05", V["influenza_FDR05"]],

    ["Frozen architecture", V["frozen_total"]],
    ["Shared core", V["shared_core"]],
    ["Influenza-amplified", V["influenza_amplified"]],
    ["RSV-amplified", V["RSV_amplified"]],

    ["Four-season evaluable all four",
     V["fourseason_evaluable"]],

    ["Four-season direction concordant all four",
     V["fourseason_concordant"]],

    ["Proteomics frozen validation genes",
     V["proteomics_validation"]],

    ["Single-cell primary cell types",
     V["singlecell_types"]],
], columns=["item", "value"])

source_path = TABDIR / "Figure1_SUBMISSION_v15_source_data.tsv"
source.to_csv(source_path, sep="\t", index=False)

# ============================================================
# PROVENANCE
# ============================================================

script_path = Path(__file__)

pd.DataFrame([{
    "file": str(script_path),
    "sha256": hashlib.sha256(
        script_path.read_bytes()
    ).hexdigest()
}]).to_csv(
    PROVDIR / "Figure1_SUBMISSION_v15_script_sha256.tsv",
    sep="\t",
    index=False
)

# ============================================================
# EXPORT
# ============================================================

fig.savefig(
    PDF,
    bbox_inches="tight",
    pad_inches=0.10
)

fig.savefig(
    PNG,
    dpi=600,
    bbox_inches="tight",
    pad_inches=0.10
)

fig.savefig(
    SVG,
    bbox_inches="tight",
    pad_inches=0.10
)

plt.close(fig)

print("=== FIGURE 1 SUBMISSION v15 COMPLETE ===")
print("PDF:", PDF)
print("PNG:", PNG)
print("SVG:", SVG)
print("SOURCE:", source_path)
