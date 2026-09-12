import pandas as pd
from pathlib import Path

ROOT = Path.home() / "Influenza_RSV_Project"

platforms = {
    "GPL10558": {
        "expression":
            ROOT / "data/processed/DS001_GSE38900/preprocessing/"
                   "DS001_GSE38900_GPL10558_expression_vsn.tsv.gz",
        "annotation":
            ROOT / "data/interim/DS001_GSE38900/annotation/"
                   "GPL10558_annotation_clean.tsv",
    },
    "GPL6884": {
        "expression":
            ROOT / "data/processed/DS001_GSE38900/preprocessing/"
                   "DS001_GSE38900_GPL6884_expression_vsn.tsv.gz",
        "annotation":
            ROOT / "data/interim/DS001_GSE38900/annotation/"
                   "GPL6884_annotation_clean.tsv",
    },
}

ANALYSIS_READY = (
    ROOT / "data/processed/DS001_GSE38900/analysis_ready"
)

RESULTS = (
    ROOT / "results/DS001_GSE38900/gene_level/tables"
)

ANALYSIS_READY.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

summary_rows = []

for platform, paths in platforms.items():

    print("\n" + "=" * 80)
    print(platform)
    print("=" * 80)

    # ----------------------------------------------------------
    # Read expression
    # ----------------------------------------------------------
    expr = pd.read_csv(paths["expression"], sep="\t")

    expr = expr.rename(
        columns={expr.columns[0]: "probe_id"}
    )

    expr["probe_id"] = expr["probe_id"].astype(str)

    sample_cols = [
        c for c in expr.columns
        if c != "probe_id"
    ]

    # ----------------------------------------------------------
    # Read annotation
    # ----------------------------------------------------------
    ann = pd.read_csv(
        paths["annotation"],
        sep="\t",
        dtype=str
    )

    ann["probe_id"] = ann["probe_id"].astype(str)

    ann = ann.drop_duplicates(
        subset="probe_id",
        keep="first"
    )

    keep_ann = [
        "probe_id",
        "gene_symbol",
        "entrez_gene_id",
        "refseq_id",
        "chromosome",
        "illumina_gene",
    ]

    keep_ann = [
        c for c in keep_ann
        if c in ann.columns
    ]

    # ----------------------------------------------------------
    # Compute IQR
    # ----------------------------------------------------------
    values = expr[sample_cols]

    q25 = values.quantile(0.25, axis=1)
    q75 = values.quantile(0.75, axis=1)

    metrics = pd.DataFrame({
        "probe_id": expr["probe_id"],
        "iqr": q75 - q25,
        "mean_expression": values.mean(axis=1),
        "variance": values.var(axis=1, ddof=1),
    })

    # ----------------------------------------------------------
    # Merge annotation
    # ----------------------------------------------------------
    merged = metrics.merge(
        ann[keep_ann],
        on="probe_id",
        how="left",
        validate="one_to_one",
    )

    merged["gene_symbol"] = (
        merged["gene_symbol"]
        .astype("string")
        .str.strip()
    )

    invalid_symbols = {
        "", "nan", "NA", "N/A",
        "---", "None", "<NA>"
    }

    usable = merged[
        merged["gene_symbol"].notna()
        & ~merged["gene_symbol"].isin(invalid_symbols)
    ].copy()

    # ----------------------------------------------------------
    # Count probes per gene
    # ----------------------------------------------------------
    counts = (
        usable.groupby("gene_symbol")["probe_id"]
        .transform("count")
    )

    usable["n_probes_for_gene"] = counts

    # ----------------------------------------------------------
    # Highest-IQR selection
    # deterministic probe_id tie-break
    # ----------------------------------------------------------
    ranked = usable.sort_values(
        by=["gene_symbol", "iqr", "probe_id"],
        ascending=[True, False, True],
        kind="mergesort",
    )

    retained = (
        ranked
        .drop_duplicates(
            subset="gene_symbol",
            keep="first"
        )
        .copy()
    )

    retained["selection_strategy"] = (
        "highest_IQR_within_gene_symbol"
    )

    retained_ids = set(retained["probe_id"])

    discarded = usable[
        ~usable["probe_id"].isin(retained_ids)
    ].copy()

    discarded["discard_reason"] = (
        "duplicate_gene_symbol_lower_IQR"
    )

    # Unannotated probes separately
    unannotated = merged[
        ~merged["probe_id"].isin(
            usable["probe_id"]
        )
    ].copy()

    unannotated["discard_reason"] = (
        "no_valid_gene_symbol"
    )

    # ----------------------------------------------------------
    # Build gene-level expression matrix
    # ----------------------------------------------------------
    gene_expr = expr[
        expr["probe_id"].isin(retained_ids)
    ].copy()

    probe_to_gene = retained[
        [
            "probe_id",
            "gene_symbol",
            "entrez_gene_id",
            "refseq_id",
            "chromosome",
            "iqr",
            "mean_expression",
            "variance",
            "n_probes_for_gene",
            "selection_strategy",
        ]
    ].copy()

    gene_expr = gene_expr.merge(
        probe_to_gene[
            [
                "probe_id",
                "gene_symbol",
                "entrez_gene_id",
            ]
        ],
        on="probe_id",
        how="left",
        validate="one_to_one",
    )

    # Put gene identifiers first
    gene_expr = gene_expr[
        [
            "gene_symbol",
            "entrez_gene_id",
            "probe_id",
        ]
        + sample_cols
    ]

    # Sort consistently
    gene_expr = gene_expr.sort_values(
        "gene_symbol"
    ).reset_index(drop=True)

    retained = retained.sort_values(
        "gene_symbol"
    )

    discarded = discarded.sort_values(
        ["gene_symbol", "probe_id"]
    )

    # ----------------------------------------------------------
    # Save outputs
    # ----------------------------------------------------------
    gene_expr.to_csv(
        ANALYSIS_READY /
        f"{platform}_gene_expression.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    retained.to_csv(
        RESULTS /
        f"{platform}_retained_probe_list.tsv",
        sep="\t",
        index=False,
    )

    discarded.to_csv(
        RESULTS /
        f"{platform}_discarded_duplicate_probe_list.tsv",
        sep="\t",
        index=False,
    )

    unannotated.to_csv(
        RESULTS /
        f"{platform}_discarded_unannotated_probe_list.tsv",
        sep="\t",
        index=False,
    )

    probe_to_gene.to_csv(
        RESULTS /
        f"{platform}_probe_to_gene_mapping.tsv",
        sep="\t",
        index=False,
    )

    # ----------------------------------------------------------
    # Validation
    # ----------------------------------------------------------
    assert gene_expr["gene_symbol"].is_unique
    assert len(gene_expr) == retained["gene_symbol"].nunique()

    n_initial = len(expr)
    n_gene = len(gene_expr)
    n_duplicate_removed = len(discarded)
    n_unannotated_removed = len(unannotated)

    summary_rows.append({
        "platform": platform,
        "initial_probe_count": n_initial,
        "sample_count": len(sample_cols),
        "annotated_probe_count": len(usable),
        "final_gene_count": n_gene,
        "duplicate_probes_removed": n_duplicate_removed,
        "unannotated_probes_removed": n_unannotated_removed,
        "collapsing_strategy": "highest_IQR",
    })

    print(f"Initial probes: {n_initial:,}")
    print(f"Samples: {len(sample_cols):,}")
    print(f"Annotated probes: {len(usable):,}")
    print(f"Final genes: {n_gene:,}")
    print(
        f"Duplicate probes removed: "
        f"{n_duplicate_removed:,}"
    )
    print(
        f"Unannotated probes removed: "
        f"{n_unannotated_removed:,}"
    )

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    RESULTS /
    "DS001_GSE38900_gene_level_construction_summary.tsv",
    sep="\t",
    index=False,
)

print("\n" + "=" * 80)
print("FINAL GENE-LEVEL SUMMARY")
print("=" * 80)

print(summary.to_string(index=False))

print("\nGene-level matrix construction completed successfully.")
