import anndata as ad
import pandas as pd
from pathlib import Path

H5AD = Path(
    "data/scrnaseq_validation/GSE283746/processed/"
    "GSE283746_COVIDRSV_scRNA.h5ad"
)

OUTDIR = Path(
    "results/scrnaseq_validation/GSE283746/isg_state/design"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

PAIRS = [
    ("Naive_CD4_T", "ISG+ Naive CD4+ T", "Naive CD4+ T"),
    ("Naive_CD8_T", "ISG+ Naive CD8+ T", "Naive CD8+ T"),
    ("Naive_B", "ISG+ Naive B", "Naive B"),
    ("TrB", "ISG+ TrB", "TrB"),
    ("Naive_Treg", "ISG+ Naive Treg", "Naive Treg"),
    ("Tfh_like", "ISG+ Tfh-like", "Tfh-like"),
    ("CD14_Monocyte", "ISG+ CD14+ Monocyte", "CD14+ Monocyte"),
    ("CD16_Monocyte", "ISG+ CD16+ Monocyte", "CD16+ Monocyte"),
]

MIN_CELLS = 25
MIN_PARTICIPANTS = 8

a = ad.read_h5ad(H5AD, backed="r")
obs = a.obs.copy()

target = obs[
    obs["Combined Condition"].isin(["RSV", "Healthy"])
].copy()

counts = (
    target.groupby(
        ["Sample", "Combined Condition", "Annotation"],
        observed=True
    )
    .size()
    .reset_index(name="n_cells")
)

all_rows = []
summary_rows = []

for pair_id, isg_state, ref_state in PAIRS:

    print("\n========================================")
    print(pair_id)
    print("========================================")
    print("ISG:", isg_state)
    print("REF:", ref_state)

    samples = (
        target[
            ["Sample", "Combined Condition"]
        ]
        .drop_duplicates()
        .copy()
    )

    isg = (
        counts[counts["Annotation"] == isg_state]
        [["Sample", "n_cells"]]
        .rename(columns={"n_cells": "ISG_cells"})
    )

    ref = (
        counts[counts["Annotation"] == ref_state]
        [["Sample", "n_cells"]]
        .rename(columns={"n_cells": "Reference_cells"})
    )

    d = (
        samples
        .merge(isg, on="Sample", how="left")
        .merge(ref, on="Sample", how="left")
    )

    d["ISG_cells"] = d["ISG_cells"].fillna(0).astype(int)
    d["Reference_cells"] = d["Reference_cells"].fillna(0).astype(int)

    d["pair_total_cells"] = (
        d["ISG_cells"] +
        d["Reference_cells"]
    )

    d["ISG_fraction"] = (
        d["ISG_cells"] /
        d["pair_total_cells"].replace(0, pd.NA)
    )

    d["ISG_pseudobulk_eligible"] = (
        d["ISG_cells"] >= MIN_CELLS
    )

    d["Reference_pseudobulk_eligible"] = (
        d["Reference_cells"] >= MIN_CELLS
    )

    d["pair_id"] = pair_id
    d["ISG_state"] = isg_state
    d["Reference_state"] = ref_state

    all_rows.append(d)

    print("\nCell-count summary:")
    print(
        d.groupby("Combined Condition")[
            ["ISG_cells", "Reference_cells"]
        ]
        .describe()
        .to_string()
    )

    print("\nParticipants with >=25 ISG-state cells:")
    print(
        d.groupby("Combined Condition")[
            "ISG_pseudobulk_eligible"
        ]
        .sum()
        .to_string()
    )

    print("\nParticipants with >=25 reference-state cells:")
    print(
        d.groupby("Combined Condition")[
            "Reference_pseudobulk_eligible"
        ]
        .sum()
        .to_string()
    )

    rsv_isg = int(
        d.loc[
            d["Combined Condition"] == "RSV",
            "ISG_pseudobulk_eligible"
        ].sum()
    )

    hc_isg = int(
        d.loc[
            d["Combined Condition"] == "Healthy",
            "ISG_pseudobulk_eligible"
        ].sum()
    )

    rsv_ref = int(
        d.loc[
            d["Combined Condition"] == "RSV",
            "Reference_pseudobulk_eligible"
        ].sum()
    )

    hc_ref = int(
        d.loc[
            d["Combined Condition"] == "Healthy",
            "Reference_pseudobulk_eligible"
        ].sum()
    )

    primary_expression_eligible = (
        rsv_isg >= MIN_PARTICIPANTS and
        hc_isg >= MIN_PARTICIPANTS and
        rsv_ref >= MIN_PARTICIPANTS and
        hc_ref >= MIN_PARTICIPANTS
    )

    summary_rows.append({
        "pair_id": pair_id,
        "ISG_state": isg_state,
        "Reference_state": ref_state,
        "RSV_ISG_ge25": rsv_isg,
        "Healthy_ISG_ge25": hc_isg,
        "RSV_reference_ge25": rsv_ref,
        "Healthy_reference_ge25": hc_ref,
        "primary_expression_eligible":
            primary_expression_eligible
    })

detail = pd.concat(all_rows, ignore_index=True)
summary = pd.DataFrame(summary_rows)

detail.to_csv(
    OUTDIR /
    "GSE283746_ISG_state_participant_counts_FROZEN_v1.0.tsv",
    sep="\t",
    index=False
)

summary.to_csv(
    OUTDIR /
    "GSE283746_ISG_state_eligibility_FROZEN_v1.0.tsv",
    sep="\t",
    index=False
)

print("\n========================================")
print("FINAL ISG-STATE ELIGIBILITY")
print("========================================")

print(summary.to_string(index=False))

print("\nEligible primary expression pairs:")
for x in summary.loc[
    summary["primary_expression_eligible"],
    "pair_id"
]:
    print(x)

a.file.close()
