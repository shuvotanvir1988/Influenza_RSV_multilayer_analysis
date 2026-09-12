from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path.home() / "Influenza_RSV_Project"

INFILE = (
    ROOT /
    "results/rnaseq_validation/influenza_SIG/validation/"
    "SIG_four_season_frozen170_validation_master.tsv"
)

OUT = (
    ROOT /
    "results/session38_additional_analysis/Session38D/tables"
)

OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INFILE, sep="\t")

assert len(df) == 170
assert df["gene_symbol"].nunique() == 170

seasons = [2018, 2019, 2020, 2022]

priority5 = {
    "IFIH1",
    "IFIT3",
    "ISG15",
    "KPNB1",
    "STAT1",
}

bridge7 = {
    "CHMP5",
    "FCGR1B",
    "HERC5",
    "HIST2H2AC",
    "KPNB1",
    "OTOF",
    "TOP2A",
}

# ------------------------------------------------------------
# DERIVED SEASONAL COUNTS
# ------------------------------------------------------------

evaluable_cols = [
    f"evaluable_{s}" for s in seasons
]

direction_cols = [
    f"direction_match_{s}" for s in seasons
]

support_cols = [
    f"FDR05_direction_match_{s}" for s in seasons
]

logfc_cols = [
    f"logFC_{s}" for s in seasons
]

df["n_seasons_evaluable_recalc"] = (
    df[evaluable_cols]
    .astype(bool)
    .sum(axis=1)
)

df["n_seasons_direction_match_recalc"] = (
    df[direction_cols]
    .fillna(0)
    .astype(float)
    .sum(axis=1)
)

df["n_supported_seasons"] = (
    df[support_cols]
    .astype(bool)
    .sum(axis=1)
)

df["direction_fraction"] = (
    df["n_seasons_direction_match_recalc"]
    / df["n_seasons_evaluable_recalc"]
)

df["supported_fraction"] = (
    df["n_supported_seasons"]
    / df["n_seasons_evaluable_recalc"]
)

# ------------------------------------------------------------
# PRIMARY CLASSIFICATION
# ------------------------------------------------------------

def classify(row):

    n_eval = int(
        row["n_seasons_evaluable_recalc"]
    )

    direction_fraction = float(
        row["direction_fraction"]
    )

    supported_fraction = float(
        row["supported_fraction"]
    )

    if n_eval < 3:
        return "INSUFFICIENT_EVALUATION"

    if direction_fraction < 1.0:
        return "DIRECTIONALLY_VARIABLE"

    if supported_fraction == 1.0:
        return "STABLE_PERSISTENT"

    if supported_fraction >= 0.75:
        return "STABLE_PREDOMINANT"

    if supported_fraction > 0:
        return "STABLE_INTERMITTENT"

    return "STABLE_DIRECTION_ONLY"


df["session38D_stability_class"] = (
    df.apply(classify, axis=1)
)

# ------------------------------------------------------------
# MAGNITUDE STABILITY
# ------------------------------------------------------------

def magnitude_metrics(row):

    vals = []

    for s in seasons:

        if bool(row[f"evaluable_{s}"]):
            x = row[f"logFC_{s}"]

            if pd.notna(x):
                vals.append(float(x))

    if len(vals) == 0:
        return pd.Series({
            "season_mean_logFC": np.nan,
            "season_median_logFC": np.nan,
            "season_SD_logFC": np.nan,
            "season_min_logFC": np.nan,
            "season_max_logFC": np.nan,
            "effect_range": np.nan,
            "absolute_mean_logFC": np.nan,
            "CV_abs": np.nan,
        })

    a = np.asarray(vals, dtype=float)

    mean = a.mean()

    sd = (
        a.std(ddof=1)
        if len(a) > 1
        else np.nan
    )

    absolute_mean = abs(mean)

    cv = (
        sd / absolute_mean
        if (
            pd.notna(sd)
            and absolute_mean > 0
        )
        else np.nan
    )

    return pd.Series({
        "season_mean_logFC":
            mean,
        "season_median_logFC":
            np.median(a),
        "season_SD_logFC":
            sd,
        "season_min_logFC":
            a.min(),
        "season_max_logFC":
            a.max(),
        "effect_range":
            a.max() - a.min(),
        "absolute_mean_logFC":
            absolute_mean,
        "CV_abs":
            cv,
    })


mag = df.apply(
    magnitude_metrics,
    axis=1
)

df = pd.concat(
    [df, mag],
    axis=1
)

# ------------------------------------------------------------
# ARCHITECTURE CLASS
# ------------------------------------------------------------

def architecture_class(row):

    if bool(
        row[
            "shared_replicated_influenza_amplified"
        ]
    ):
        return "Influenza_amplified_shared"

    if bool(
        row[
            "shared_replicated_RSV_amplified"
        ]
    ):
        return "RSV_amplified_shared"

    if bool(
        row["shared_replicated_core"]
    ):
        return "Shared_core"

    return str(
        row.get(
            "Figure7A_class",
            "OTHER"
        )
    )


df["architecture_class_38D"] = (
    df.apply(
        architecture_class,
        axis=1
    )
)

df["in_priority5"] = (
    df["gene_symbol"]
    .isin(priority5)
)

df["in_bridge7"] = (
    df["gene_symbol"]
    .isin(bridge7)
)

# ------------------------------------------------------------
# QA AGAINST PREEXISTING SUMMARY FIELDS
# ------------------------------------------------------------

assert (
    df[
        "n_seasons_evaluable_recalc"
    ]
    == df[
        "n_seasons_evaluable"
    ]
).all()

assert np.allclose(
    df[
        "n_seasons_direction_match_recalc"
    ],
    df[
        "n_seasons_direction_match"
    ],
    equal_nan=True,
)

# ------------------------------------------------------------
# SUMMARY TABLES
# ------------------------------------------------------------

class_summary = (
    df[
        "session38D_stability_class"
    ]
    .value_counts()
    .rename_axis(
        "stability_class"
    )
    .reset_index(
        name="n_genes"
    )
)

class_summary[
    "fraction_of_170"
] = (
    class_summary["n_genes"]
    / len(df)
)

architecture_summary_rows = []

for cls, z in df.groupby(
    "architecture_class_38D",
    observed=True
):

    architecture_summary_rows.append({
        "architecture_class":
            cls,
        "n_genes":
            len(z),
        "n_evaluable_ge3":
            int(
                (
                    z[
                        "n_seasons_evaluable_recalc"
                    ]
                    >= 3
                ).sum()
            ),
        "n_stable_persistent":
            int(
                (
                    z[
                        "session38D_stability_class"
                    ]
                    == "STABLE_PERSISTENT"
                ).sum()
            ),
        "n_stable_predominant":
            int(
                (
                    z[
                        "session38D_stability_class"
                    ]
                    == "STABLE_PREDOMINANT"
                ).sum()
            ),
        "n_stable_intermittent":
            int(
                (
                    z[
                        "session38D_stability_class"
                    ]
                    == "STABLE_INTERMITTENT"
                ).sum()
            ),
        "n_stable_direction_only":
            int(
                (
                    z[
                        "session38D_stability_class"
                    ]
                    == "STABLE_DIRECTION_ONLY"
                ).sum()
            ),
        "n_directionally_variable":
            int(
                (
                    z[
                        "session38D_stability_class"
                    ]
                    == "DIRECTIONALLY_VARIABLE"
                ).sum()
            ),
        "median_supported_fraction":
            z[
                "supported_fraction"
            ].median(),
        "median_CV_abs":
            z[
                "CV_abs"
            ].median(),
        "median_effect_range":
            z[
                "effect_range"
            ].median(),
    })


architecture_summary = pd.DataFrame(
    architecture_summary_rows
)

# Priority and bridge subsets are generated AFTER
# magnitude-stability ranks are added below.

# ------------------------------------------------------------
# MAGNITUDE RANKS
# ------------------------------------------------------------

eligible = (
    df[
        "n_seasons_evaluable_recalc"
    ]
    >= 3
)

df["CV_abs_rank_stability"] = np.nan
df["effect_range_rank_stability"] = np.nan

df.loc[
    eligible,
    "CV_abs_rank_stability"
] = (
    df.loc[
        eligible,
        "CV_abs"
    ]
    .rank(
        method="min",
        ascending=True
    )
)

df.loc[
    eligible,
    "effect_range_rank_stability"
] = (
    df.loc[
        eligible,
        "effect_range"
    ]
    .rank(
        method="min",
        ascending=True
    )
)

# Generate subsets only after all derived ranking columns exist.
priority = df[
    df["in_priority5"]
].copy()

bridge = df[
    df["in_bridge7"]
].copy()

assert len(priority) == 5
assert set(priority["gene_symbol"]) == {
    "IFIH1", "IFIT3", "ISG15", "KPNB1", "STAT1"
}

assert len(bridge) == 7
assert set(bridge["gene_symbol"]) == {
    "CHMP5", "FCGR1B", "HERC5",
    "HIST2H2AC", "KPNB1", "OTOF", "TOP2A"
}

# ------------------------------------------------------------
# OUTPUTS
# ------------------------------------------------------------

master_cols = [
    "gene_symbol",
    "architecture_class_38D",
    "flu_logFC",
    "logFC_2018",
    "FDR_2018",
    "logFC_2019",
    "FDR_2019",
    "logFC_2020",
    "FDR_2020",
    "logFC_2022",
    "FDR_2022",
    "n_seasons_evaluable_recalc",
    "n_seasons_direction_match_recalc",
    "direction_fraction",
    "n_supported_seasons",
    "supported_fraction",
    "session38D_stability_class",
    "season_mean_logFC",
    "season_median_logFC",
    "season_SD_logFC",
    "season_min_logFC",
    "season_max_logFC",
    "effect_range",
    "absolute_mean_logFC",
    "CV_abs",
    "CV_abs_rank_stability",
    "effect_range_rank_stability",
    "in_priority5",
    "in_bridge7",
]

master = df[
    master_cols
].copy()

master.to_csv(
    OUT /
    "Session38D_CROSSSEASON_STABILITY_MASTER_v1.0.tsv",
    sep="\t",
    index=False
)

class_summary.to_csv(
    OUT /
    "Session38D_STABILITY_CLASS_SUMMARY_v1.0.tsv",
    sep="\t",
    index=False
)

architecture_summary.to_csv(
    OUT /
    "Session38D_ARCHITECTURE_CLASS_STABILITY_SUMMARY_v1.0.tsv",
    sep="\t",
    index=False
)

priority[
    master_cols
].to_csv(
    OUT /
    "Session38D_PRIORITY5_STABILITY_v1.0.tsv",
    sep="\t",
    index=False
)

bridge[
    master_cols
].to_csv(
    OUT /
    "Session38D_BRIDGE7_STABILITY_v1.0.tsv",
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# TERMINAL REPORT
# ------------------------------------------------------------

print(
    "=== SESSION 38D CROSS-SEASON STABILITY COMPLETE ==="
)

print()
print(
    "=== STABILITY CLASS SUMMARY ==="
)

print(
    class_summary
    .to_string(
        index=False
    )
)

print()
print(
    "=== ARCHITECTURE CLASS SUMMARY ==="
)

print(
    architecture_summary
    .to_string(
        index=False
    )
)

print()
print(
    "=== PRIORITY 5 ==="
)

show = [
    "gene_symbol",
    "architecture_class_38D",
    "n_seasons_evaluable_recalc",
    "n_supported_seasons",
    "supported_fraction",
    "session38D_stability_class",
    "season_mean_logFC",
    "effect_range",
    "CV_abs",
]

print(
    priority[
        show
    ]
    .sort_values(
        "gene_symbol"
    )
    .to_string(
        index=False
    )
)

print()
print(
    "=== BRIDGE 7 ==="
)

print(
    bridge[
        show
    ]
    .sort_values(
        "gene_symbol"
    )
    .to_string(
        index=False
    )
)

print()
print(
    "=== MOST MAGNITUDE-STABLE GENES ==="
)

print(
    master.loc[
        master[
            "n_seasons_evaluable_recalc"
        ] >= 3,
        [
            "gene_symbol",
            "architecture_class_38D",
            "session38D_stability_class",
            "CV_abs",
            "effect_range",
        ]
    ]
    .sort_values(
        [
            "CV_abs",
            "effect_range",
        ]
    )
    .head(20)
    .to_string(
        index=False
    )
)

print()
print(
    "=== MOST MAGNITUDE-VARIABLE GENES ==="
)

print(
    master.loc[
        master[
            "n_seasons_evaluable_recalc"
        ] >= 3,
        [
            "gene_symbol",
            "architecture_class_38D",
            "session38D_stability_class",
            "CV_abs",
            "effect_range",
        ]
    ]
    .sort_values(
        [
            "CV_abs",
            "effect_range",
        ],
        ascending=[
            False,
            False
        ]
    )
    .head(20)
    .to_string(
        index=False
    )
)

print()
print(
    "Written outputs:"
)

for f in [
    "Session38D_CROSSSEASON_STABILITY_MASTER_v1.0.tsv",
    "Session38D_STABILITY_CLASS_SUMMARY_v1.0.tsv",
    "Session38D_ARCHITECTURE_CLASS_STABILITY_SUMMARY_v1.0.tsv",
    "Session38D_PRIORITY5_STABILITY_v1.0.tsv",
    "Session38D_BRIDGE7_STABILITY_v1.0.tsv",
]:
    print(
        OUT / f
    )
