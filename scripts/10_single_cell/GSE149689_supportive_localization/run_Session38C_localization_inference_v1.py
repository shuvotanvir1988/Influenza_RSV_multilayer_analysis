from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest

ROOT = Path.home() / "Influenza_RSV_Project"

TAB = (
    ROOT /
    "results/session38_additional_analysis/Session38C/tables"
)

DESIGN = (
    ROOT /
    "results/session38_additional_analysis/Session38C/design"
)

FROZEN170 = (
    ROOT /
    "results/figure2_influenza_centered/"
    "frozen_submission_v5/tables/"
    "Figure2C_v8_replicated_170_gene_architecture.tsv"
)

OUT = TAB
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# INPUTS
# ============================================================

long = pd.read_csv(
    TAB / "Session38C_PSEUDOBULK_RAW_COUNTS_LONG_v1.0.tsv",
    sep="\t"
)

libs = pd.read_csv(
    TAB / "Session38C_PSEUDOBULK_FULL_LIBRARY_SIZES_v1.0.tsv",
    sep="\t"
)

targets = pd.read_csv(
    DESIGN / "Session38C_RESOLVED_TARGET_MANIFEST_v1.0.tsv",
    sep="\t"
)

celltypes = pd.read_csv(
    DESIGN / "Session38C_PRIMARY_CELLTYPES_FROZEN_v1.0.tsv",
    sep="\t"
)

bulk = pd.read_csv(
    FROZEN170,
    sep="\t"
)

assert len(bulk) == 170

bulk = bulk[
    [
        "gene_symbol",
        "flu_logFC",
        "architecture_class",
    ]
].copy()

targets = targets[
    targets["analysis_measurable"]
].copy()

assert len(targets) == 161

eligible_celltypes = (
    celltypes.loc[
        celltypes["primary_eligible"],
        "cell_type"
    ]
    .astype(str)
    .tolist()
)

assert len(eligible_celltypes) == 5

flu_donors = [
    "Flu 1", "Flu 2", "Flu 3",
    "Flu 4", "Flu 5",
]

normal_donors = [
    "Normal 1", "Normal 2",
    "Normal 3", "Normal 4",
]

all_donors = flu_donors + normal_donors

priority5 = {
    "IFIH1",
    "IFIT3",
    "ISG15",
    "KPNB1",
    "STAT1",
}

bridge5 = {
    "CHMP5",
    "HERC5",
    "KPNB1",
    "OTOF",
    "TOP2A",
}

# ============================================================
# MERGE FULL LIBRARY SIZES
# ============================================================

libs_use = libs[
    [
        "donor_id",
        "disease",
        "cell_type",
        "full_raw_library_size",
    ]
].copy()

dat = long.merge(
    libs_use,
    on=[
        "donor_id",
        "disease",
        "cell_type",
    ],
    how="left",
    validate="many_to_one"
)

assert dat["full_raw_library_size"].notna().all()

dat["logCPM"] = np.log2(
    (
        dat["raw_count"]
        / dat["full_raw_library_size"]
    )
    * 1e6
    + 1
)

# Attach frozen bulk architecture.
dat = dat.merge(
    bulk,
    on="gene_symbol",
    how="left",
    validate="many_to_one"
)

assert dat["flu_logFC"].notna().all()

dat["bulk_sign"] = np.sign(
    dat["flu_logFC"]
)

# ============================================================
# EXACT PERMUTATION UTILITIES
# ============================================================

all_assignments = list(
    combinations(
        range(9),
        5
    )
)

assert len(all_assignments) == 126

def exact_perm_p(values, labels):
    """
    Exact two-sided permutation test for difference in means.
    labels must contain 5 influenza and 4 normal.
    """

    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)

    obs_flu = labels == "influenza"
    obs_norm = labels == "normal"

    observed = (
        values[obs_flu].mean()
        - values[obs_norm].mean()
    )

    perm_effects = []

    for idx in all_assignments:

        mask = np.zeros(9, dtype=bool)
        mask[list(idx)] = True

        eff = (
            values[mask].mean()
            - values[~mask].mean()
        )

        perm_effects.append(eff)

    perm_effects = np.asarray(
        perm_effects
    )

    p = np.mean(
        np.abs(perm_effects)
        >= abs(observed) - 1e-12
    )

    return observed, p

def bh_adjust(pvalues):
    p = np.asarray(
        pvalues,
        dtype=float
    )

    n = len(p)
    order = np.argsort(p)
    ranked = p[order]

    q_ranked = (
        ranked
        * n
        / np.arange(1, n + 1)
    )

    q_ranked = np.minimum.accumulate(
        q_ranked[::-1]
    )[::-1]

    q_ranked = np.minimum(
        q_ranked,
        1.0
    )

    q = np.empty(n)
    q[order] = q_ranked

    return q

# ============================================================
# GENE-LEVEL CELL-TYPE EFFECTS
# ============================================================

gene_effect_rows = []

for ct in eligible_celltypes:

    z = dat[
        dat["cell_type"] == ct
    ].copy()

    for gene in targets["gene_symbol"]:

        g = z[
            z["gene_symbol"] == gene
        ].copy()

        g = (
            g.set_index("donor_id")
            .reindex(all_donors)
            .reset_index()
        )

        assert len(g) == 9
        assert g["logCPM"].notna().all()

        labels = np.array(
            [
                "influenza"
                if d in flu_donors
                else "normal"
                for d in all_donors
            ]
        )

        effect, p = exact_perm_p(
            g["logCPM"].to_numpy(),
            labels
        )

        bulk_fc = float(
            g["flu_logFC"].iloc[0]
        )

        gene_effect_rows.append({
            "cell_type": ct,
            "gene_symbol": gene,
            "sc_effect_logCPM":
                effect,
            "bulk_flu_logFC":
                bulk_fc,
            "same_direction":
                bool(
                    np.sign(effect)
                    == np.sign(bulk_fc)
                ),
            "exact_perm_p":
                p,
            "in_priority5":
                gene in priority5,
            "in_bridge5":
                gene in bridge5,
        })

gene_effects = pd.DataFrame(
    gene_effect_rows
)

assert len(gene_effects) == (
    5 * 161
)

# Full-gene FDR within cell type.
gene_effects["gene_FDR_within_celltype"] = np.nan

for ct in eligible_celltypes:

    idx = (
        gene_effects["cell_type"]
        == ct
    )

    gene_effects.loc[
        idx,
        "gene_FDR_within_celltype"
    ] = bh_adjust(
        gene_effects.loc[
            idx,
            "exact_perm_p"
        ].to_numpy()
    )

# ============================================================
# ARCHITECTURE CORRELATION + DIRECTION CONCORDANCE
# ============================================================

architecture_rows = []

for ct in eligible_celltypes:

    z = gene_effects[
        gene_effects["cell_type"] == ct
    ].copy()

    rho, p_rho = spearmanr(
        z["sc_effect_logCPM"],
        z["bulk_flu_logFC"]
    )

    same = int(
        z["same_direction"].sum()
    )

    n = len(z)

    bt = binomtest(
        same,
        n,
        p=0.5,
        alternative="two-sided"
    )

    architecture_rows.append({
        "cell_type": ct,
        "n_genes": n,
        "spearman_rho":
            rho,
        "spearman_p":
            p_rho,
        "direction_concordant_n":
            same,
        "direction_total_n":
            n,
        "direction_concordance":
            same / n,
        "direction_binomial_p":
            bt.pvalue,
    })

architecture = pd.DataFrame(
    architecture_rows
)

architecture["spearman_FDR"] = bh_adjust(
    architecture[
        "spearman_p"
    ].to_numpy()
)

architecture["direction_binomial_FDR"] = bh_adjust(
    architecture[
        "direction_binomial_p"
    ].to_numpy()
)

# ============================================================
# DONOR-LEVEL SIGNED ARCHITECTURE SCORES
# ============================================================

score_rows = []

for ct in eligible_celltypes:

    z = dat[
        dat["cell_type"] == ct
    ].copy()

    mat = (
        z.pivot(
            index="donor_id",
            columns="gene_symbol",
            values="logCPM"
        )
        .reindex(
            index=all_donors,
            columns=targets["gene_symbol"]
        )
    )

    assert mat.shape == (9, 161)
    assert mat.notna().all().all()

    # Standardize each gene across 9 donors.
    mu = mat.mean(axis=0)
    sd = mat.std(
        axis=0,
        ddof=1
    )

    usable = sd > 0

    matz = (
        mat.loc[:, usable]
        - mu[usable]
    ) / sd[usable]

    bulk_sign = (
        bulk.set_index("gene_symbol")
        .loc[
            matz.columns,
            "flu_logFC"
        ]
        .apply(np.sign)
    )

    signed = (
        matz
        * bulk_sign
    )

    scores = signed.mean(
        axis=1
    )

    for donor in all_donors:

        score_rows.append({
            "cell_type": ct,
            "donor_id": donor,
            "disease":
                "influenza"
                if donor in flu_donors
                else "normal",
            "architecture_score":
                float(scores.loc[donor]),
            "n_genes_used":
                int(usable.sum()),
        })

scores = pd.DataFrame(
    score_rows
)

score_test_rows = []

for ct in eligible_celltypes:

    z = (
        scores[
            scores["cell_type"] == ct
        ]
        .set_index("donor_id")
        .reindex(all_donors)
        .reset_index()
    )

    labels = z[
        "disease"
    ].to_numpy()

    effect, p = exact_perm_p(
        z["architecture_score"]
        .to_numpy(),
        labels
    )

    score_test_rows.append({
        "cell_type": ct,
        "architecture_score_difference":
            effect,
        "exact_perm_p":
            p,
    })

score_tests = pd.DataFrame(
    score_test_rows
)

score_tests["FDR"] = bh_adjust(
    score_tests[
        "exact_perm_p"
    ].to_numpy()
)

# ============================================================
# PRIORITY 5
# ============================================================

priority = gene_effects[
    gene_effects[
        "gene_symbol"
    ].isin(priority5)
].copy()

assert len(priority) == 25

priority["priority_FDR_25"] = bh_adjust(
    priority[
        "exact_perm_p"
    ].to_numpy()
)

# ============================================================
# BRIDGE 5
# ============================================================

bridge = gene_effects[
    gene_effects[
        "gene_symbol"
    ].isin(bridge5)
].copy()

assert len(bridge) == 25

bridge["bridge_FDR_25"] = bh_adjust(
    bridge[
        "exact_perm_p"
    ].to_numpy()
)

# ============================================================
# DONOR-LEVEL PRIORITY/BRIDGE LOGCPM
# ============================================================

target_donor = dat[
    dat["gene_symbol"].isin(
        priority5 | bridge5
    )
][
    [
        "donor_id",
        "disease",
        "cell_type",
        "gene_symbol",
        "logCPM",
        "raw_count",
        "full_raw_library_size",
    ]
].copy()

# ============================================================
# OUTPUTS
# ============================================================

gene_effects.to_csv(
    OUT /
    "Session38C_GENE_LEVEL_CELLTYPE_EFFECTS_v1.0.tsv",
    sep="\t",
    index=False
)

architecture.to_csv(
    OUT /
    "Session38C_ARCHITECTURE_LOCALIZATION_v1.0.tsv",
    sep="\t",
    index=False
)

scores.to_csv(
    OUT /
    "Session38C_DONOR_ARCHITECTURE_SCORES_v1.0.tsv",
    sep="\t",
    index=False
)

score_tests.to_csv(
    OUT /
    "Session38C_ARCHITECTURE_SCORE_TESTS_v1.0.tsv",
    sep="\t",
    index=False
)

priority.to_csv(
    OUT /
    "Session38C_PRIORITY5_CELLTYPE_EFFECTS_v1.0.tsv",
    sep="\t",
    index=False
)

bridge.to_csv(
    OUT /
    "Session38C_BRIDGE5_CELLTYPE_EFFECTS_v1.0.tsv",
    sep="\t",
    index=False
)

target_donor.to_csv(
    OUT /
    "Session38C_TARGET_DONOR_LOGCPM_v1.0.tsv",
    sep="\t",
    index=False
)

# ============================================================
# TERMINAL REPORT
# ============================================================

print("=== SESSION 38C LOCALIZATION INFERENCE COMPLETE ===")
print()

print("=== ARCHITECTURE LOCALIZATION ===")
print(
    architecture
    .sort_values(
        "spearman_rho",
        ascending=False
    )
    .to_string(index=False)
)

print()
print("=== ARCHITECTURE SCORE TESTS ===")
print(
    score_tests
    .sort_values(
        "architecture_score_difference",
        ascending=False
    )
    .to_string(index=False)
)

print()
print("=== PRIORITY 5 ===")
print(
    priority[
        [
            "cell_type",
            "gene_symbol",
            "sc_effect_logCPM",
            "bulk_flu_logFC",
            "same_direction",
            "exact_perm_p",
            "priority_FDR_25",
        ]
    ]
    .sort_values(
        [
            "priority_FDR_25",
            "cell_type",
            "gene_symbol",
        ]
    )
    .to_string(index=False)
)

print()
print("=== BRIDGE 5 ===")
print(
    bridge[
        [
            "cell_type",
            "gene_symbol",
            "sc_effect_logCPM",
            "bulk_flu_logFC",
            "same_direction",
            "exact_perm_p",
            "bridge_FDR_25",
        ]
    ]
    .sort_values(
        [
            "bridge_FDR_25",
            "cell_type",
            "gene_symbol",
        ]
    )
    .to_string(index=False)
)

print()
print("Written outputs:")
for f in [
    "Session38C_GENE_LEVEL_CELLTYPE_EFFECTS_v1.0.tsv",
    "Session38C_ARCHITECTURE_LOCALIZATION_v1.0.tsv",
    "Session38C_DONOR_ARCHITECTURE_SCORES_v1.0.tsv",
    "Session38C_ARCHITECTURE_SCORE_TESTS_v1.0.tsv",
    "Session38C_PRIORITY5_CELLTYPE_EFFECTS_v1.0.tsv",
    "Session38C_BRIDGE5_CELLTYPE_EFFECTS_v1.0.tsv",
    "Session38C_TARGET_DONOR_LOGCPM_v1.0.tsv",
]:
    print(OUT / f)
