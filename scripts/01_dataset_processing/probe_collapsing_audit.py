import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path.home() / "Influenza_RSV_Project"

platforms = {
    "GPL10558": {
        "expression": ROOT / "data/processed/DS001_GSE38900/preprocessing/DS001_GSE38900_GPL10558_expression_vsn.tsv.gz",
        "annotation": ROOT / "data/interim/DS001_GSE38900/annotation/GPL10558_annotation_clean.tsv",
    },
    "GPL6884": {
        "expression": ROOT / "data/processed/DS001_GSE38900/preprocessing/DS001_GSE38900_GPL6884_expression_vsn.tsv.gz",
        "annotation": ROOT / "data/interim/DS001_GSE38900/annotation/GPL6884_annotation_clean.tsv",
    },
}

OUT = ROOT / "results/DS001_GSE38900/gene_level/tables"
OUT.mkdir(parents=True, exist_ok=True)

summary_rows = []
agreement_rows = []

for platform, files in platforms.items():

    print(f"\n{'='*80}")
    print(f"Processing {platform}")
    print(f"{'='*80}")

    # ----------------------------------------------------------
    # Read expression matrix
    # ----------------------------------------------------------
    expr = pd.read_csv(files["expression"], sep="\t")

    # Session 14 stored probe IDs as dataframe index,
    # which becomes "Unnamed: 0" after export.
    first_col = expr.columns[0]
    expr = expr.rename(columns={first_col: "probe_id"})
    expr["probe_id"] = expr["probe_id"].astype(str)

    sample_cols = [c for c in expr.columns if c != "probe_id"]

    # Ensure numeric expression
    expr[sample_cols] = expr[sample_cols].apply(pd.to_numeric, errors="coerce")

    print(f"Expression probes: {len(expr):,}")
    print(f"Samples: {len(sample_cols):,}")

    # ----------------------------------------------------------
    # Read annotation
    # ----------------------------------------------------------
    ann = pd.read_csv(files["annotation"], sep="\t", dtype=str)
    ann["probe_id"] = ann["probe_id"].astype(str)

    # Keep one annotation record per probe
    ann = ann.drop_duplicates(subset="probe_id", keep="first")

    # ----------------------------------------------------------
    # Probe-level expression statistics
    # ----------------------------------------------------------
    values = expr[sample_cols]

    expr["mean_expression"] = values.mean(axis=1)
    expr["variance"] = values.var(axis=1, ddof=1)
    expr["sd"] = values.std(axis=1, ddof=1)
    expr["q25"] = values.quantile(0.25, axis=1)
    expr["q75"] = values.quantile(0.75, axis=1)
    expr["iqr"] = expr["q75"] - expr["q25"]

    # ----------------------------------------------------------
    # Join expression with annotation
    # ----------------------------------------------------------
    keep_ann = [
        "probe_id",
        "gene_symbol",
        "entrez_gene_id",
        "refseq_id",
        "illumina_gene",
        "chromosome",
        "has_gene_symbol",
        "has_entrez_gene_id",
        "probes_per_gene_symbol",
    ]

    keep_ann = [c for c in keep_ann if c in ann.columns]

    audit = expr[
        [
            "probe_id",
            "mean_expression",
            "variance",
            "sd",
            "q25",
            "q75",
            "iqr",
        ]
    ].merge(
        ann[keep_ann],
        on="probe_id",
        how="left",
        validate="one_to_one",
    )

    # ----------------------------------------------------------
    # Clean gene symbol / Entrez fields
    # ----------------------------------------------------------
    audit["gene_symbol"] = audit["gene_symbol"].astype("string").str.strip()
    audit["entrez_gene_id"] = audit["entrez_gene_id"].astype("string").str.strip()

    invalid_symbols = {"", "nan", "NA", "N/A", "---", "None", "<NA>"}
    audit["valid_gene_symbol"] = (
        audit["gene_symbol"].notna()
        & ~audit["gene_symbol"].isin(invalid_symbols)
    )

    invalid_entrez = {"", "nan", "NA", "N/A", "---", "None", "<NA>"}
    audit["valid_entrez_id"] = (
        audit["entrez_gene_id"].notna()
        & ~audit["entrez_gene_id"].isin(invalid_entrez)
    )

    # Primary collapsing universe = probes with usable gene symbols.
    annotated = audit.loc[audit["valid_gene_symbol"]].copy()

    annotated["n_probes_for_gene"] = (
        annotated.groupby("gene_symbol")["probe_id"].transform("count")
    )

    # ----------------------------------------------------------
    # Rank probes within each gene
    # Deterministic tie-breaking by probe_id
    # ----------------------------------------------------------
    def select_probe(df, metric):
        tmp = df.sort_values(
            by=["gene_symbol", metric, "probe_id"],
            ascending=[True, False, True],
            kind="mergesort",
        )

        selected = (
            tmp.groupby("gene_symbol", as_index=False)
            .first()[["gene_symbol", "probe_id", metric]]
            .rename(
                columns={
                    "probe_id": f"probe_{metric}",
                    metric: f"value_{metric}",
                }
            )
        )
        return selected

    sel_iqr = select_probe(annotated, "iqr")
    sel_var = select_probe(annotated, "variance")
    sel_mean = select_probe(annotated, "mean_expression")

    compare = (
        sel_iqr
        .merge(sel_var, on="gene_symbol", how="inner")
        .merge(sel_mean, on="gene_symbol", how="inner")
    )

    gene_counts = (
        annotated.groupby("gene_symbol")
        .agg(
            n_probes=("probe_id", "count"),
            any_valid_entrez=("valid_entrez_id", "max"),
        )
        .reset_index()
    )

    compare = compare.merge(gene_counts, on="gene_symbol", how="left")

    compare["iqr_equals_variance"] = (
        compare["probe_iqr"] == compare["probe_variance"]
    )

    compare["iqr_equals_mean"] = (
        compare["probe_iqr"] == compare["probe_mean_expression"]
    )

    compare["variance_equals_mean"] = (
        compare["probe_variance"] == compare["probe_mean_expression"]
    )

    compare["all_three_agree"] = (
        compare["iqr_equals_variance"]
        & compare["iqr_equals_mean"]
    )

    duplicated = compare.loc[compare["n_probes"] > 1].copy()

    # ----------------------------------------------------------
    # Output probe-level audit
    # ----------------------------------------------------------
    audit_out = OUT / f"{platform}_probe_level_collapsing_audit.tsv.gz"
    audit.to_csv(audit_out, sep="\t", index=False, compression="gzip")

    compare_out = OUT / f"{platform}_probe_selection_strategy_comparison.tsv"
    compare.to_csv(compare_out, sep="\t", index=False)

    duplicated_out = OUT / f"{platform}_duplicated_gene_strategy_comparison.tsv"
    duplicated.to_csv(duplicated_out, sep="\t", index=False)

    # ----------------------------------------------------------
    # Summary statistics
    # ----------------------------------------------------------
    n_expression = len(expr)
    n_matched_annotation = audit["gene_symbol"].notna().sum()
    n_symbol_probes = audit["valid_gene_symbol"].sum()
    n_entrez_probes = audit["valid_entrez_id"].sum()

    n_genes = annotated["gene_symbol"].nunique()
    n_single_probe_genes = (
        gene_counts["n_probes"].eq(1).sum()
    )
    n_multi_probe_genes = (
        gene_counts["n_probes"].gt(1).sum()
    )

    if len(duplicated) > 0:
        iqr_var_pct = 100 * duplicated["iqr_equals_variance"].mean()
        iqr_mean_pct = 100 * duplicated["iqr_equals_mean"].mean()
        var_mean_pct = 100 * duplicated["variance_equals_mean"].mean()
        all_three_pct = 100 * duplicated["all_three_agree"].mean()
    else:
        iqr_var_pct = np.nan
        iqr_mean_pct = np.nan
        var_mean_pct = np.nan
        all_three_pct = np.nan

    summary_rows.append({
        "platform": platform,
        "expression_probe_count": n_expression,
        "annotation_gene_symbol_nonmissing": n_matched_annotation,
        "probes_with_valid_gene_symbol": n_symbol_probes,
        "probes_with_valid_entrez_id": n_entrez_probes,
        "unique_gene_symbols": n_genes,
        "single_probe_genes": n_single_probe_genes,
        "multi_probe_genes": n_multi_probe_genes,
        "probes_removed_if_one_per_symbol": n_symbol_probes - n_genes,
        "unannotated_probes_excluded": n_expression - n_symbol_probes,
    })

    agreement_rows.append({
        "platform": platform,
        "multi_probe_genes": len(duplicated),
        "iqr_variance_same_probe_percent": iqr_var_pct,
        "iqr_mean_same_probe_percent": iqr_mean_pct,
        "variance_mean_same_probe_percent": var_mean_pct,
        "all_three_same_probe_percent": all_three_pct,
    })

    print(f"Valid symbol probes: {n_symbol_probes:,}")
    print(f"Unique gene symbols: {n_genes:,}")
    print(f"Single-probe genes: {n_single_probe_genes:,}")
    print(f"Multi-probe genes: {n_multi_probe_genes:,}")
    print()
    print("Agreement among duplicated genes:")
    print(f"IQR vs variance: {iqr_var_pct:.2f}%")
    print(f"IQR vs mean: {iqr_mean_pct:.2f}%")
    print(f"Variance vs mean: {var_mean_pct:.2f}%")
    print(f"All three: {all_three_pct:.2f}%")

summary = pd.DataFrame(summary_rows)
agreement = pd.DataFrame(agreement_rows)

summary.to_csv(
    OUT / "DS001_GSE38900_probe_collapsing_inventory.tsv",
    sep="\t",
    index=False,
)

agreement.to_csv(
    OUT / "DS001_GSE38900_probe_strategy_agreement.tsv",
    sep="\t",
    index=False,
)

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(summary.to_string(index=False))

print("\n" + "="*80)
print("STRATEGY AGREEMENT")
print("="*80)
print(agreement.to_string(index=False))

print("\nProbe-collapsing audit completed successfully.")
