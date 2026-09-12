#!/usr/bin/env python3

from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT = Path.home() / "Influenza_RSV_Project"
BASE = ROOT / "results/DS001_GSE38900/pathway_analysis"
FIGDIR = ROOT / "results/figure3_influenza_centered/figures"
TABDIR = ROOT / "results/figure3_influenza_centered/tables"
FIGDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)

HALLMARK_INF = BASE / "tables/fgsea_hallmark/GPL6884_InfluenzaA_vs_control_Hallmark_fgsea.tsv"
HALLMARK_DIRECT = BASE / "tables/fgsea_hallmark/GPL6884_InfluenzaA_vs_RSVacute_Hallmark_fgsea.tsv"
DOMAIN_REP = BASE / "tables/preregistered_domains/DS001_GSE38900_preregistered_domain_representative_pathways.tsv"
DOMAIN_SUM = BASE / "tables/preregistered_domains/DS001_GSE38900_preregistered_domain_summary.tsv"
LE_REC = BASE / "tables/robustness_audit/DS001_GSE38900_corrected_leading_edge_gene_recurrence.tsv"

for f in [HALLMARK_INF, HALLMARK_DIRECT, DOMAIN_REP, DOMAIN_SUM, LE_REC]:
    if not f.exists():
        raise FileNotFoundError(f)

hall = pd.read_csv(HALLMARK_INF, sep="\t")
direct = pd.read_csv(HALLMARK_DIRECT, sep="\t")
rep = pd.read_csv(DOMAIN_REP, sep="\t")
dom = pd.read_csv(DOMAIN_SUM, sep="\t")
rec = pd.read_csv(LE_REC, sep="\t")

# ============================================================
# SUBMISSION TYPOGRAPHY — LOCKED TO FIGURE 2 v5 STANDARD
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 16.0,
    "axes.titlesize": 22.0,
    "axes.titleweight": "bold",
    "axes.labelsize": 19.0,
    "axes.labelweight": "bold",
    "xtick.labelsize": 16.0,
    "ytick.labelsize": 16.0,
    "legend.fontsize": 16.0,
    "axes.linewidth": 2.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

C_GRAY = "#B8BDC5"
C_DARK = "#333333"
C_BLUE = "#4C72B0"
C_RED = "#C44E52"

HALLMARK_LABELS = {
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α response",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ response",
    "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNFα signaling via NF-κB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": "IL-6/JAK/STAT3 signaling",
    "HALLMARK_COMPLEMENT": "Complement",
    "HALLMARK_IL2_STAT5_SIGNALING": "IL-2/STAT5 signaling",
    "HALLMARK_APOPTOSIS": "Apoptosis",
    "HALLMARK_G2M_CHECKPOINT": "G2/M checkpoint",
    "HALLMARK_E2F_TARGETS": "E2F targets",
    "HALLMARK_MITOTIC_SPINDLE": "Mitotic spindle",
    "HALLMARK_MTORC1_SIGNALING": "mTORC1 signaling",
}

def hallmark_label(x):
    return HALLMARK_LABELS.get(
        x,
        x.replace("HALLMARK_", "").replace("_", " ").title()
    )

# ============================================================
# LAYOUT — COMPACT COLUMNS, SAFE ROW SPACING
# ============================================================

fig = plt.figure(figsize=(19.5, 13.5))
gs = GridSpec(
    2, 2,
    figure=fig,
    width_ratios=[1.06, 1.10],
    height_ratios=[1.0, 1.0],
    wspace=0.32,
    hspace=0.56,
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

# ============================================================
# PANEL A — HALLMARK LANDSCAPE
# ============================================================

sig = hall[hall["significant_FDR05"].astype(str).str.upper().eq("TRUE")].copy()
sig["abs_NES"] = sig["NES"].abs()

immune_programs = [
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_APOPTOSIS",
]

immune = sig[sig["pathway"].isin(immune_programs)].copy()
context = (
    sig[~sig["pathway"].isin(immune_programs)]
    .sort_values("abs_NES", ascending=False)
    .head(4)
)

a = (
    pd.concat([immune, context], ignore_index=True)
    .drop_duplicates("pathway")
    .sort_values("NES", ascending=True)
)

y = np.arange(len(a))
axA.hlines(y, 0, a["NES"], color="0.52", linewidth=2.9, zorder=1)
maxabs = max(abs(a["NES"].min()), abs(a["NES"].max()))
axA.scatter(
    a["NES"], y,
    s=np.clip(-np.log10(a["padj"]) * 11.5, 70, 170),
    c=a["NES"], cmap="coolwarm", vmin=-maxabs, vmax=maxabs,
    edgecolor="black", linewidth=0.7, zorder=3
)
axA.axvline(0, color="0.15", linewidth=1.8)
axA.set_yticks(y)
axA.set_yticklabels([hallmark_label(x) for x in a["pathway"]])
for tick in axA.get_yticklabels():
    tick.set_fontsize(17.0)
    tick.set_fontweight("bold")
axA.set_xlabel("Normalized enrichment score (NES)")
axA.spines[["top", "right"]].set_visible(False)
axA.grid(axis="x", alpha=0.24, linewidth=1.05)

# ============================================================
# PANEL B — PREREGISTERED DOMAIN MATRIX
# ============================================================

domain_order = [
    "Interferon",
    "Innate_Cytokine",
    "Antigen_Presentation",
    "Neutrophil_Monocyte",
    "Adaptive_Immunity",
    "Metabolism",
    "Translation",
]

domain_labels = {
    "Interferon": "Interferon",
    "Innate_Cytokine": "Innate / cytokine",
    "Antigen_Presentation": "Antigen presentation",
    "Neutrophil_Monocyte": "Neutrophil / monocyte",
    "Adaptive_Immunity": "Adaptive immunity",
    "Metabolism": "Metabolism",
    "Translation": "Translation",
}

d = dom.set_index("preregistered_domain").reindex(domain_order).reset_index()
matrix = d[["median_NES_Influenza_GPL6884", "median_NES_Influenza_vs_RSV"]].to_numpy(dtype=float)
vmax = np.nanmax(np.abs(matrix))

im = axB.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
axB.set_xticks([0, 1])
axB.set_xticklabels(["Influenza A\nvs control", "Influenza A\nvs RSV"])
axB.set_yticks(np.arange(len(d)))
axB.set_yticklabels([domain_labels[x] for x in d["preregistered_domain"]])
for tick in axB.get_yticklabels():
    tick.set_fontsize(17.0)
    tick.set_fontweight("bold")
for tick in axB.get_xticklabels():
    tick.set_fontsize(16.0)
    tick.set_fontweight("bold")

for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        value = matrix[i, j]
        text_color = "white" if abs(value) > 0.58 * vmax else "black"
        axB.text(
            j, i, f"{value:.2f}",
            ha="center", va="center", color=text_color,
            fontsize=17.0, fontweight="bold"
        )

for i, row in d.iterrows():
    axB.text(
        1.70, i, f'n={int(row["pathways_in_domain"])}',
        ha="left", va="center", fontsize=16.5,
        fontweight="bold", color="0.22", clip_on=False
    )

axB.set_xlim(-0.5, 2.22)
for spine in axB.spines.values():
    spine.set_visible(False)

cbar = fig.colorbar(im, ax=axB, fraction=0.040, pad=0.065)
cbar.set_label("Median NES", fontweight="bold")
cbar.ax.tick_params(
    width=1.5,
    length=5,
    labelsize=16.5
)
for tick in cbar.ax.get_yticklabels():
    tick.set_fontweight("bold")

# ============================================================
# PANEL C — DIRECT INFLUENZA VS RSV ANTIVIRAL AMPLIFICATION
# ============================================================

selected_hallmarks = [
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]

c1 = direct[direct["pathway"].isin(selected_hallmarks)][["pathway", "NES", "padj"]].copy()
reactome_targets = [
    "REACTOME_MODULATION_OF_HOST_RESPONSES_BY_IFN_STIMULATED_GENES",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
    "REACTOME_INTERFERON_SIGNALING",
]

c2 = rep[rep["pathway"].isin(reactome_targets)][
    ["pathway", "NES_Influenza_vs_RSV_GPL6884", "padj_Influenza_vs_RSV_GPL6884"]
].copy()
c2 = c2.rename(columns={
    "NES_Influenza_vs_RSV_GPL6884": "NES",
    "padj_Influenza_vs_RSV_GPL6884": "padj",
})
c = pd.concat([c1, c2], ignore_index=True)

label_map_c = {
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-α response",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ response",
    "REACTOME_MODULATION_OF_HOST_RESPONSES_BY_IFN_STIMULATED_GENES": "IFN-stimulated host-response\nmodulation",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING": "IFN-α/β signaling",
    "REACTOME_INTERFERON_SIGNALING": "Interferon signaling",
}

c["display"] = c["pathway"].map(label_map_c)
c = c.sort_values("NES", ascending=True)
y = np.arange(len(c))

axC.hlines(y, 0, c["NES"], color="0.50", linewidth=2.9, zorder=1)
axC.scatter(
    c["NES"], y,
    s=np.clip(-np.log10(c["padj"]) * 15.5, 80, 170),
    color=C_BLUE, edgecolor="black", linewidth=0.7, zorder=3
)
axC.axvline(0, color="0.15", linewidth=1.8)
axC.set_yticks(y)
axC.set_yticklabels(c["display"])
axC.margins(y=0.035)
for tick in axC.get_yticklabels():
    tick.set_fontsize(17.0)
    tick.set_fontweight("bold")
axC.set_xlabel("NES, influenza A vs RSV")
axC.spines[["top", "right"]].set_visible(False)
axC.grid(axis="x", alpha=0.22, linewidth=1.05)

x_pad = 0.07
for yi, val in zip(y, c["NES"]):
    axC.text(
        val + x_pad, yi, f"{val:.2f}",
        va="center", ha="left",
        fontsize=17.0, fontweight="bold"
    )

axC.set_xlim(min(-0.12, c["NES"].min() - 0.10), c["NES"].max() + 0.48)

# ============================================================
# PANEL D — RECURRENT LEADING-EDGE MACHINERY
# ============================================================

r = rec[
    (rec["contrast"] == "Influenza_GPL6884") &
    (rec["theme"] == "Interferon") &
    (rec["pathways_containing_gene"] >= 4)
].copy()

r = r.sort_values(
    ["pathways_containing_gene", "median_pathway_NES", "leading_edge_gene"],
    ascending=[True, True, True]
)

y = np.arange(len(r))
axD.barh(y, r["pathways_containing_gene"], height=0.74, color=C_BLUE, edgecolor="0.22", linewidth=0.50)
axD.set_yticks(y)
axD.set_yticklabels(r["leading_edge_gene"])
for tick in axD.get_yticklabels():
    tick.set_fontsize(18.0)
    tick.set_fontweight("bold")
axD.set_xlabel("Influenza interferon pathways containing gene")
axD.spines[["top", "right"]].set_visible(False)
axD.set_xlim(0, r["pathways_containing_gene"].max() + 0.95)
axD.grid(axis="x", alpha=0.22, linewidth=1.05)

for yi, n in zip(y, r["pathways_containing_gene"]):
    axD.text(
        n + 0.10, yi, f"{int(n)}",
        va="center", fontsize=18.0, fontweight="bold"
    )

# ============================================================
# PANEL HEADINGS — SAME BASELINE, SAME BOLD TITLE SIZE
# ============================================================

PANEL_LABEL_SIZE = 25.0
PANEL_TITLE_SIZE = 23.0
PANEL_HEADING_Y = 1.080
PANEL_LABEL_X = -0.095
PANEL_TITLE_X = 0.00

def add_panel_heading(ax, letter, title):
    ax.text(
        PANEL_LABEL_X, PANEL_HEADING_Y, letter,
        transform=ax.transAxes,
        fontsize=PANEL_LABEL_SIZE, fontweight="bold",
        ha="left", va="bottom", clip_on=False
    )
    ax.text(
        PANEL_TITLE_X, PANEL_HEADING_Y, title,
        transform=ax.transAxes,
        fontsize=PANEL_TITLE_SIZE, fontweight=700,
        ha="left", va="bottom", clip_on=False,
        linespacing=1.05
    )

add_panel_heading(axA, "A", "Influenza A Hallmark program landscape")
add_panel_heading(axB, "B", "Preregistered biological-domain architecture")
add_panel_heading(axC, "C", "Influenza amplification of antiviral programs")
add_panel_heading(axD, "D", "Recurrent leading-edge antiviral machinery")

# No figure-level title/footer in manuscript-facing submission version.
fig.subplots_adjust(top=0.94, bottom=0.075, left=0.095, right=0.985)

# ============================================================
# TYPOGRAPHY LOCK — AXES / TICKS
# ============================================================

for ax in [axA, axB, axC, axD]:
    ax.tick_params(axis="both", width=1.8, length=6.5)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")

for ax in [axA, axC, axD]:
    ax.xaxis.label.set_fontweight("bold")
    ax.yaxis.label.set_fontweight("bold")

# ============================================================
# SCIENTIFIC QC / FROZEN OUTPUT TABLES
# ============================================================

# Save exactly what is displayed in Panels B-D for figure provenance.
d.to_csv(TABDIR / "Figure3_SUBMISSION_v6_domain_matrix.tsv", sep="\t", index=False)
c[["pathway", "display", "NES", "padj"]].to_csv(
    TABDIR / "Figure3_SUBMISSION_v6_antiviral_programs.tsv", sep="\t", index=False
)
r[["leading_edge_gene", "pathways_containing_gene", "median_pathway_NES"]].to_csv(
    TABDIR / "Figure3_SUBMISSION_v6_leading_edge_recurrence.tsv", sep="\t", index=False
)

# ============================================================
# SAVE
# ============================================================

pdf = FIGDIR / "Figure3_SUBMISSION_v6_influenza_pathway_architecture.pdf"
png = FIGDIR / "Figure3_SUBMISSION_v6_influenza_pathway_architecture.png"
svg = FIGDIR / "Figure3_SUBMISSION_v6_influenza_pathway_architecture.svg"

fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=600, bbox_inches="tight")
fig.savefig(svg, format="svg", bbox_inches="tight")
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

prov_rows = []
for path, role in [
    (HALLMARK_INF, "Hallmark influenza A vs control enrichment"),
    (HALLMARK_DIRECT, "Hallmark influenza A vs RSV enrichment"),
    (DOMAIN_REP, "Preregistered representative pathways"),
    (DOMAIN_SUM, "Preregistered biological-domain summary"),
    (LE_REC, "Corrected leading-edge recurrence audit"),
    (pdf, "Figure 3 submission v6 PDF"),
    (png, "Figure 3 submission v6 PNG"),
    (svg, "Figure 3 submission v6 editable SVG"),
]:
    prov_rows.append({
        "file": str(path.relative_to(ROOT)),
        "sha256": sha256(path),
        "role": role,
    })

prov = pd.DataFrame(prov_rows)
prov_path = TABDIR / "Figure3_SUBMISSION_v6_provenance_sha256.tsv"
prov.to_csv(prov_path, sep="\t", index=False)

# ============================================================
# REPORT QC
# ============================================================

print("=== FIGURE 3 SUBMISSION v6 QC ===")
print("Hallmark programs displayed:", len(a))
print("Preregistered domains displayed:", len(d))
print("Direct antiviral programs displayed:", len(c))
print("Leading-edge genes retained (recurrence >=4):", len(r))
print()
print("Panel B median NES:")
print(d[["preregistered_domain", "median_NES_Influenza_GPL6884", "median_NES_Influenza_vs_RSV"]].to_string(index=False))
print()
print("Panel C direct influenza-vs-RSV NES:")
print(c[["display", "NES"]].sort_values("NES", ascending=False).to_string(index=False))
print()
print("Panel D recurrence:")
print(r[["leading_edge_gene", "pathways_containing_gene"]].sort_values(
    ["pathways_containing_gene", "leading_edge_gene"], ascending=[False, True]
).to_string(index=False))
print()
print("=== WRITTEN ===")
print(pdf)
print(png)
print(svg)
print(TABDIR / "Figure3_SUBMISSION_v6_domain_matrix.tsv")
print(TABDIR / "Figure3_SUBMISSION_v6_antiviral_programs.tsv")
print(TABDIR / "Figure3_SUBMISSION_v6_leading_edge_recurrence.tsv")
print(prov_path)
print("\nSTATUS: FIGURE 3 SUBMISSION v6 COMPLETE")
