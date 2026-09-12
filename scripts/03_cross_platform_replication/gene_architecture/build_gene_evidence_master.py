#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("results/DS001_GSE38900")

DE = ROOT / "differential_expression/tables"
PA = ROOT / "pathway_analysis/tables/robustness_audit"
OUT = ROOT / "signature_analysis/tables"
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    "flu":
        DE / "GPL6884_InfluenzaA_vs_control_full_DE.tsv",
    "rsv6884":
        DE / "GPL6884_RSV_acute_vs_control_full_DE.tsv",
    "rsv10558":
        DE / "GPL10558_RSV_acute_vs_control_full_DE.tsv",
    "direct":
        DE / "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv",
    "replication":
        DE / "DS001_GSE38900_RSV_cross_platform_replication.tsv",
    "leading_edge":
        PA / "DS001_GSE38900_corrected_leading_edge_genes.tsv.gz",
}

FDR = 0.05


def load_de(path, prefix):
    df = pd.read_csv(path, sep="\t")

    keep = [
        "gene_symbol",
        "entrez_gene_id",
        "logFC",
        "P.Value",
        "adj.P.Val",
    ]

    df = df[keep].copy()

    df = df.rename(columns={
        "entrez_gene_id": f"{prefix}_entrez_gene_id",
        "logFC": f"{prefix}_logFC",
        "P.Value": f"{prefix}_P",
        "adj.P.Val": f"{prefix}_FDR",
    })

    return df


# ---------------------------------------------------------
# 1. Primary GPL6884 universe
# ---------------------------------------------------------

flu = load_de(FILES["flu"], "flu")
rsv6884 = load_de(FILES["rsv6884"], "rsv6884")
direct = load_de(FILES["direct"], "flu_vs_rsv")

master = (
    flu
    .merge(rsv6884, on="gene_symbol", how="outer", validate="one_to_one")
    .merge(direct, on="gene_symbol", how="outer", validate="one_to_one")
)

# ---------------------------------------------------------
# 2. GPL10558 RSV evidence
# ---------------------------------------------------------

rsv10558 = load_de(FILES["rsv10558"], "rsv10558")

master = master.merge(
    rsv10558,
    on="gene_symbol",
    how="left",
    validate="one_to_one"
)

# ---------------------------------------------------------
# 3. Existing Session 16 cross-platform replication
# ---------------------------------------------------------

rep = pd.read_csv(FILES["replication"], sep="\t")

rep_keep = [
    "gene_symbol",
    "significant_GPL10558",
    "significant_GPL6884",
    "significant_both",
    "same_direction",
    "minimum_abs_logFC",
]

rep = rep[rep_keep].copy()

rep = rep.rename(columns={
    "significant_GPL10558": "rsv_rep_sig_GPL10558",
    "significant_GPL6884": "rsv_rep_sig_GPL6884",
    "significant_both": "rsv_rep_significant_both",
    "same_direction": "rsv_rep_same_direction",
    "minimum_abs_logFC": "rsv_rep_minimum_abs_logFC",
})

master = master.merge(
    rep,
    on="gene_symbol",
    how="left",
    validate="one_to_one"
)

master["rsv_crossplatform_evaluable"] = \
    master["rsv_rep_same_direction"].notna()

master["rsv_replicated"] = (
    master["rsv_rep_significant_both"].fillna(False).astype(bool)
    &
    master["rsv_rep_same_direction"].fillna(False).astype(bool)
)

# ---------------------------------------------------------
# 4. Primary statistical indicators
# ---------------------------------------------------------

master["flu_significant"] = master["flu_FDR"] < FDR
master["rsv6884_significant"] = master["rsv6884_FDR"] < FDR
master["direct_significant"] = master["flu_vs_rsv_FDR"] < FDR

master["infection_same_direction"] = (
    np.sign(master["flu_logFC"]) ==
    np.sign(master["rsv6884_logFC"])
)

master["infection_opposite_direction"] = (
    np.sign(master["flu_logFC"]) ==
    -np.sign(master["rsv6884_logFC"])
)

# Shared on the common GPL6884 cohort/platform
master["shared_GPL6884"] = (
    master["flu_significant"]
    &
    master["rsv6884_significant"]
    &
    master["infection_same_direction"]
)

# Strongest shared category: RSV additionally reproduces on GPL10558
master["shared_replicated"] = (
    master["shared_GPL6884"]
    &
    master["rsv_replicated"]
)

# ---------------------------------------------------------
# 5. Direct pathogen-bias evidence
#
# Direct contrast is:
# Influenza A acute vs RSV acute
#
# Positive logFC = influenza-biased
# Negative logFC = RSV-biased
# ---------------------------------------------------------

master["influenza_biased_direct"] = (
    master["direct_significant"]
    &
    (master["flu_vs_rsv_logFC"] > 0)
)

master["rsv_biased_direct"] = (
    master["direct_significant"]
    &
    (master["flu_vs_rsv_logFC"] < 0)
)

# Same-direction infection response but significantly
# different magnitude between pathogens
master["shared_quantitatively_divergent"] = (
    master["shared_GPL6884"]
    &
    master["direct_significant"]
)

# ---------------------------------------------------------
# 6. Corrected Session 17 leading-edge evidence
# ---------------------------------------------------------

le = pd.read_csv(FILES["leading_edge"], sep="\t")

# One gene can occur in many pathways/themes/contrasts.
le_summary = (
    le.groupby("leading_edge_gene")
      .agg(
          leading_edge_occurrences=("leading_edge_gene", "size"),
          leading_edge_contrasts=("contrast",
                                  lambda x: "|".join(sorted(set(map(str, x))))),
          leading_edge_themes=("theme",
                               lambda x: "|".join(sorted(set(map(str, x))))),
          leading_edge_pathway_count=("pathway", "nunique"),
      )
      .reset_index()
      .rename(columns={"leading_edge_gene": "gene_symbol"})
)

master = master.merge(
    le_summary,
    on="gene_symbol",
    how="left",
    validate="one_to_one"
)

master["session17_leading_edge"] = \
    master["leading_edge_occurrences"].notna()

master["leading_edge_occurrences"] = \
    master["leading_edge_occurrences"].fillna(0).astype(int)

master["leading_edge_pathway_count"] = \
    master["leading_edge_pathway_count"].fillna(0).astype(int)

# ---------------------------------------------------------
# 7. Evidence-class assignment
#
# Classification is deliberately hierarchical.
# Pathway membership does NOT determine statistical class.
# ---------------------------------------------------------

def classify(row):

    if row["shared_replicated"]:
        if row["direct_significant"]:
            if row["flu_vs_rsv_logFC"] > 0:
                return "shared_replicated_influenza_biased"
            elif row["flu_vs_rsv_logFC"] < 0:
                return "shared_replicated_RSV_biased"
        return "shared_replicated"

    if row["shared_GPL6884"]:
        if not row["rsv_crossplatform_evaluable"]:
            return "shared_GPL6884_not_crossplatform_evaluable"

        if row["direct_significant"]:
            if row["flu_vs_rsv_logFC"] > 0:
                return "shared_nonreplicated_influenza_biased"
            elif row["flu_vs_rsv_logFC"] < 0:
                return "shared_nonreplicated_RSV_biased"

        return "shared_GPL6884_nonreplicated"

    if row["influenza_biased_direct"]:
        return "influenza_biased"

    if row["rsv_biased_direct"]:
        return "RSV_biased"

    return "other"


master["signature_class"] = master.apply(classify, axis=1)

# ---------------------------------------------------------
# 8. Stable ordering
# ---------------------------------------------------------

master = master.sort_values("gene_symbol").reset_index(drop=True)

# ---------------------------------------------------------
# 9. Write master table
# ---------------------------------------------------------

master_file = OUT / "DS001_GSE38900_gene_evidence_master.tsv"
master.to_csv(master_file, sep="\t", index=False)

# ---------------------------------------------------------
# 10. Summary
# ---------------------------------------------------------

summary_rows = [
    ("master_gene_universe", len(master)),
    ("flu_significant_FDR05", int(master["flu_significant"].sum())),
    ("rsv6884_significant_FDR05", int(master["rsv6884_significant"].sum())),
    ("direct_significant_FDR05", int(master["direct_significant"].sum())),
    ("rsv_crossplatform_evaluable",
        int(master["rsv_crossplatform_evaluable"].sum())),
    ("rsv_replicated",
        int(master["rsv_replicated"].sum())),
    ("shared_GPL6884",
        int(master["shared_GPL6884"].sum())),
    ("shared_replicated",
        int(master["shared_replicated"].sum())),
    ("influenza_biased_direct",
        int(master["influenza_biased_direct"].sum())),
    ("rsv_biased_direct",
        int(master["rsv_biased_direct"].sum())),
    ("shared_quantitatively_divergent",
        int(master["shared_quantitatively_divergent"].sum())),
    ("session17_leading_edge",
        int(master["session17_leading_edge"].sum())),
]

summary = pd.DataFrame(summary_rows, columns=["metric", "value"])

summary.to_csv(
    OUT / "DS001_GSE38900_gene_evidence_master_summary.tsv",
    sep="\t",
    index=False
)

class_summary = (
    master["signature_class"]
    .value_counts(dropna=False)
    .rename_axis("signature_class")
    .reset_index(name="n_genes")
)

class_summary.to_csv(
    OUT / "DS001_GSE38900_signature_class_counts.tsv",
    sep="\t",
    index=False
)

print("MASTER TABLE:", master_file)
print("Genes:", len(master))

print("\nSUMMARY")
print(summary.to_string(index=False))

print("\nSIGNATURE CLASSES")
print(class_summary.to_string(index=False))
