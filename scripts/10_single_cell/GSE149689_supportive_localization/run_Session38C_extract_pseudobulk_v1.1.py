
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import cellxgene_census

ROOT = Path.home() / "Influenza_RSV_Project"

CENSUS_VERSION = "2025-11-08"
DATASET_ID = "de2c780c-1747-40bd-9ccf-9588ec186cee"

DESIGN = (
    ROOT /
    "results/session38_additional_analysis/Session38C/design"
)

OUT = (
    ROOT /
    "results/session38_additional_analysis/Session38C/tables"
)

OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Frozen design
# ------------------------------------------------------------

targets = pd.read_csv(
    DESIGN / "Session38C_RESOLVED_TARGET_MANIFEST_v1.0.tsv",
    sep="\t"
)

targets = targets[
    targets["analysis_measurable"]
].copy()

assert len(targets) == 161
assert targets["gene_symbol"].nunique() == 161
assert targets["analysis_feature_symbol"].nunique() == 161

celltypes = pd.read_csv(
    DESIGN / "Session38C_PRIMARY_CELLTYPES_FROZEN_v1.0.tsv",
    sep="\t"
)

celltypes = celltypes.loc[
    celltypes["primary_eligible"],
    "cell_type"
].astype(str).tolist()

assert len(celltypes) == 5

influenza_donors = [
    "Flu 1",
    "Flu 2",
    "Flu 3",
    "Flu 4",
    "Flu 5",
]

normal_donors = [
    "Normal 1",
    "Normal 2",
    "Normal 3",
    "Normal 4",
]

donors = influenza_donors + normal_donors

feature_symbols = (
    targets["analysis_feature_symbol"]
    .astype(str)
    .tolist()
)

print("=== SESSION 38C PSEUDOBULK EXTRACTION ===")
print("Census version:", CENSUS_VERSION)
print("Dataset:", DATASET_ID)
print("Frozen measurable targets:", len(feature_symbols))
print("Frozen cell types:", len(celltypes))
print("Eligible donors:", len(donors))
print()

# ------------------------------------------------------------
# Query Census
# ------------------------------------------------------------

with cellxgene_census.open_soma(
    census_version=CENSUS_VERSION
) as census:

    exp = census["census_data"]["homo_sapiens"]

    obs = (
        exp["obs"]
        .read(
            value_filter=f'dataset_id == "{DATASET_ID}"'
        )
        .concat()
        .to_pandas()
    )

    obs = obs[
        obs["donor_id"].isin(donors)
        & obs["cell_type"].isin(celltypes)
        & obs["disease"].isin(["influenza", "normal"])
    ].copy()

    print("Selected cells:", len(obs))

    assert set(obs["donor_id"]) == set(donors)
    assert set(obs["cell_type"]) == set(celltypes)

    # --------------------------------------------------------
    # Frozen Ensembl-level feature lookup
    # --------------------------------------------------------

    feature_mapping = pd.read_csv(
        DESIGN /
        "Session38C_FROZEN_ENSEMBL_FEATURE_MAPPING_v1.0.tsv",
        sep="\\t"
    )

    assert len(feature_mapping) == 161
    assert feature_mapping["frozen_gene_symbol"].nunique() == 161
    assert feature_mapping["soma_joinid"].nunique() == 161

    # Ensure mapping order follows resolved target manifest.
    feature_mapping = (
        targets[
            ["gene_symbol", "analysis_feature_symbol"]
        ]
        .merge(
            feature_mapping,
            left_on=[
                "gene_symbol",
                "analysis_feature_symbol"
            ],
            right_on=[
                "frozen_gene_symbol",
                "analysis_feature_symbol"
            ],
            how="left",
            validate="one_to_one"
        )
    )

    assert feature_mapping["soma_joinid"].notna().all()

    var_joinids = (
        feature_mapping["soma_joinid"]
        .astype(np.int64)
        .to_numpy()
    )

    obs_joinids = (
        obs["soma_joinid"]
        .astype(np.int64)
        .to_numpy()
    )

    print("Selected features:", len(var_joinids))

    # --------------------------------------------------------
    # Retrieve targeted raw-count matrix
    # --------------------------------------------------------

    raw = exp["ms"]["RNA"]["X"]["raw"]

    coo = (
        raw
        .read(
            coords=(
                obs_joinids,
                var_joinids
            )
        )
        .tables()
        .concat()
        .to_pandas()
    )

print("Retrieved nonzero entries:", len(coo))

# ------------------------------------------------------------
# Convert global SOMA joinids to compact matrix indices
# ------------------------------------------------------------

obs_order = (
    obs[
        [
            "soma_joinid",
            "donor_id",
            "disease",
            "cell_type",
        ]
    ]
    .drop_duplicates("soma_joinid")
    .reset_index(drop=True)
)

obs_index = {
    int(x): i
    for i, x in enumerate(obs_order["soma_joinid"])
}

var_index = {
    int(joinid): i
    for i, joinid in enumerate(var_joinids)
}

row = coo["soma_dim_0"].map(obs_index).to_numpy()
col = coo["soma_dim_1"].map(var_index).to_numpy()
data = coo["soma_data"].to_numpy()

if pd.isna(row).any() or pd.isna(col).any():
    raise RuntimeError(
        "Unexpected SOMA coordinate outside requested dimensions."
    )

X = sp.coo_matrix(
    (
        data,
        (
            row.astype(int),
            col.astype(int)
        )
    ),
    shape=(
        len(obs_order),
        len(feature_symbols)
    )
).tocsr()

print("Compact matrix shape:", X.shape)
print("Compact matrix nnz:", X.nnz)

# ------------------------------------------------------------
# Cell-level QC metadata
# ------------------------------------------------------------

cell_counts = (
    obs_order
    .groupby(
        ["donor_id", "disease", "cell_type"],
        observed=True
    )
    .size()
    .reset_index(name="n_cells")
)

expected = pd.MultiIndex.from_product(
    [
        donors,
        celltypes
    ],
    names=[
        "donor_id",
        "cell_type"
    ]
).to_frame(index=False)

disease_map = {
    **{x: "influenza" for x in influenza_donors},
    **{x: "normal" for x in normal_donors},
}

expected["disease"] = expected["donor_id"].map(disease_map)

cell_counts = expected.merge(
    cell_counts,
    on=[
        "donor_id",
        "cell_type",
        "disease"
    ],
    how="left"
)

cell_counts["n_cells"] = (
    cell_counts["n_cells"]
    .fillna(0)
    .astype(int)
)

assert len(cell_counts) == 45
assert (cell_counts["n_cells"] >= 20).all()

# ------------------------------------------------------------
# Pseudobulk aggregation
# ------------------------------------------------------------

pb_rows = []

for _, meta in cell_counts.iterrows():

    donor = meta["donor_id"]
    disease = meta["disease"]
    cell_type = meta["cell_type"]

    mask = (
        (obs_order["donor_id"] == donor)
        & (obs_order["cell_type"] == cell_type)
    ).to_numpy()

    idx = np.where(mask)[0]

    assert len(idx) == meta["n_cells"]

    sums = np.asarray(
        X[idx, :].sum(axis=0)
    ).ravel()

    row_dict = {
        "donor_id": donor,
        "disease": disease,
        "cell_type": cell_type,
        "n_cells": int(meta["n_cells"]),
        "target_library_size": float(sums.sum()),
        "n_targets_detected": int((sums > 0).sum()),
    }

    for gene, value in zip(
        feature_symbols,
        sums
    ):
        row_dict[gene] = value

    pb_rows.append(row_dict)

pb = pd.DataFrame(pb_rows)

assert len(pb) == 45

# ------------------------------------------------------------
# Rename analysis symbols back to frozen symbols
# ------------------------------------------------------------

analysis_to_frozen = dict(
    zip(
        targets["analysis_feature_symbol"],
        targets["gene_symbol"]
    )
)

rename_map = {
    x: analysis_to_frozen[x]
    for x in feature_symbols
}

pb_frozen = pb.rename(columns=rename_map)

# ------------------------------------------------------------
# Long-form pseudobulk table
# ------------------------------------------------------------

metadata_cols = [
    "donor_id",
    "disease",
    "cell_type",
    "n_cells",
    "target_library_size",
    "n_targets_detected",
]

long = pb.melt(
    id_vars=metadata_cols,
    value_vars=feature_symbols,
    var_name="analysis_feature_symbol",
    value_name="raw_count"
)

mapping = targets[
    [
        "gene_symbol",
        "analysis_feature_symbol",
        "mapping_status",
        "in_priority5",
        "in_bridge7",
    ]
].copy()

long = long.merge(
    mapping,
    on="analysis_feature_symbol",
    how="left",
    validate="many_to_one"
)

assert long["gene_symbol"].notna().all()

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

cell_counts.to_csv(
    OUT /
    "Session38C_PSEUDOBULK_CELL_COUNTS_v1.0.tsv",
    sep="\t",
    index=False
)

pb_frozen.to_csv(
    OUT /
    "Session38C_PSEUDOBULK_RAW_COUNTS_WIDE_v1.0.tsv",
    sep="\t",
    index=False
)

long.to_csv(
    OUT /
    "Session38C_PSEUDOBULK_RAW_COUNTS_LONG_v1.0.tsv",
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

print()
print("=== PSEUDOBULK QC ===")

print("Pseudobulk samples:", len(pb_frozen))
print(
    "Expected donor x cell-type samples:",
    len(donors) * len(celltypes)
)

print(
    "Minimum cells/sample:",
    cell_counts["n_cells"].min()
)

print(
    "Maximum cells/sample:",
    cell_counts["n_cells"].max()
)

print(
    "Minimum targets detected:",
    pb["n_targets_detected"].min()
)

print(
    "Maximum targets detected:",
    pb["n_targets_detected"].max()
)

print()
print("=== CELL COUNTS ===")

print(
    cell_counts
    .sort_values(
        ["cell_type", "donor_id"]
    )
    .to_string(index=False)
)

print()
print("=== TARGET DETECTION QC ===")

print(
    pb[
        [
            "donor_id",
            "disease",
            "cell_type",
            "n_cells",
            "target_library_size",
            "n_targets_detected",
        ]
    ]
    .sort_values(
        ["cell_type", "donor_id"]
    )
    .to_string(index=False)
)

print()
print("SESSION 38C PSEUDOBULK EXTRACTION COMPLETE")
