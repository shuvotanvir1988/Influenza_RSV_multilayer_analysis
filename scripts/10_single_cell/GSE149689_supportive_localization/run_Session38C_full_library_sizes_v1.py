from pathlib import Path
import pandas as pd
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

celltypes = pd.read_csv(
    DESIGN / "Session38C_PRIMARY_CELLTYPES_FROZEN_v1.0.tsv",
    sep="\t"
)

celltypes = (
    celltypes.loc[
        celltypes["primary_eligible"],
        "cell_type"
    ]
    .astype(str)
    .tolist()
)

flu = [
    "Flu 1", "Flu 2", "Flu 3",
    "Flu 4", "Flu 5",
]

normal = [
    "Normal 1", "Normal 2",
    "Normal 3", "Normal 4",
]

donors = flu + normal

with cellxgene_census.open_soma(
    census_version=CENSUS_VERSION
) as census:

    obs = (
        census["census_data"]["homo_sapiens"]["obs"]
        .read(
            value_filter=f'dataset_id == "{DATASET_ID}"',
            column_names=[
                "soma_joinid",
                "donor_id",
                "disease",
                "cell_type",
                "raw_sum",
                "n_measured_vars",
            ]
        )
        .concat()
        .to_pandas()
    )

obs = obs[
    obs["donor_id"].isin(donors)
    & obs["cell_type"].isin(celltypes)
    & obs["disease"].isin(["influenza", "normal"])
].copy()

assert len(obs) == 21623

summary = (
    obs.groupby(
        ["donor_id", "disease", "cell_type"],
        observed=True
    )
    .agg(
        n_cells=("soma_joinid", "size"),
        full_raw_library_size=("raw_sum", "sum"),
        median_cell_raw_sum=("raw_sum", "median"),
        min_measured_vars=("n_measured_vars", "min"),
        max_measured_vars=("n_measured_vars", "max"),
    )
    .reset_index()
)

assert len(summary) == 45
assert (summary["n_cells"] >= 20).all()
assert (summary["full_raw_library_size"] > 0).all()

old = pd.read_csv(
    OUT / "Session38C_PSEUDOBULK_CELL_COUNTS_v1.0.tsv",
    sep="\t"
)

check = old.merge(
    summary[
        [
            "donor_id",
            "disease",
            "cell_type",
            "n_cells",
        ]
    ],
    on=["donor_id", "disease", "cell_type"],
    suffixes=("_old", "_new"),
    validate="one_to_one",
)

assert (
    check["n_cells_old"]
    == check["n_cells_new"]
).all()

summary.to_csv(
    OUT /
    "Session38C_PSEUDOBULK_FULL_LIBRARY_SIZES_v1.0.tsv",
    sep="\t",
    index=False
)

print("=== SESSION 38C FULL LIBRARY SIZE QC ===")
print()
print("Cells:", len(obs))
print("Pseudobulks:", len(summary))
print()

print(
    summary.sort_values(
        ["cell_type", "donor_id"]
    ).to_string(index=False)
)

print()
print("Library-size range:")
print(
    summary["full_raw_library_size"]
    .describe()
    .to_string()
)

print()
print("ALL FULL-LIBRARY QC CHECKS PASSED")
