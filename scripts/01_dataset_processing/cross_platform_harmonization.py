import pandas as pd
from pathlib import Path

ROOT = Path.home() / "Influenza_RSV_Project"

INPUT_DIR = ROOT / "data/processed/DS001_GSE38900/analysis_ready"
TABLE_DIR = ROOT / "results/DS001_GSE38900/gene_level/tables"

TABLE_DIR.mkdir(parents=True, exist_ok=True)

files = {
    "GPL10558": INPUT_DIR / "GPL10558_gene_expression.tsv.gz",
    "GPL6884": INPUT_DIR / "GPL6884_gene_expression.tsv.gz",
}

dfs = {}

for platform, path in files.items():
    df = pd.read_csv(path, sep="\t", dtype={"gene_symbol": str, "entrez_gene_id": str})
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.strip()
    df["entrez_gene_id"] = df["entrez_gene_id"].astype("string").str.strip()
    dfs[platform] = df

g1 = dfs["GPL10558"][["gene_symbol", "entrez_gene_id", "probe_id"]].copy()
g2 = dfs["GPL6884"][["gene_symbol", "entrez_gene_id", "probe_id"]].copy()

genes1 = set(g1["gene_symbol"])
genes2 = set(g2["gene_symbol"])

shared = sorted(genes1 & genes2)
only_10558 = sorted(genes1 - genes2)
only_6884 = sorted(genes2 - genes1)

print("=" * 80)
print("GENE SYMBOL OVERLAP")
print("=" * 80)
print(f"GPL10558 genes: {len(genes1):,}")
print(f"GPL6884 genes: {len(genes2):,}")
print(f"Shared genes: {len(shared):,}")
print(f"GPL10558-specific genes: {len(only_10558):,}")
print(f"GPL6884-specific genes: {len(only_6884):,}")

# ----------------------------------------------------------
# Shared-gene mapping
# ----------------------------------------------------------
shared_df = (
    g1[g1["gene_symbol"].isin(shared)]
    .rename(columns={
        "entrez_gene_id": "GPL10558_entrez_gene_id",
        "probe_id": "GPL10558_probe_id",
    })
    .merge(
        g2[g2["gene_symbol"].isin(shared)]
        .rename(columns={
            "entrez_gene_id": "GPL6884_entrez_gene_id",
            "probe_id": "GPL6884_probe_id",
        }),
        on="gene_symbol",
        how="inner",
        validate="one_to_one",
    )
)

def normalize_entrez(x):
    if pd.isna(x):
        return None
    x = str(x).strip()
    if x in {"", "nan", "NA", "N/A", "---", "<NA>", "None"}:
        return None
    return x

shared_df["GPL10558_entrez_gene_id"] = shared_df["GPL10558_entrez_gene_id"].map(normalize_entrez)
shared_df["GPL6884_entrez_gene_id"] = shared_df["GPL6884_entrez_gene_id"].map(normalize_entrez)

shared_df["entrez_both_present"] = (
    shared_df["GPL10558_entrez_gene_id"].notna()
    & shared_df["GPL6884_entrez_gene_id"].notna()
)

shared_df["entrez_match"] = (
    shared_df["entrez_both_present"]
    & (
        shared_df["GPL10558_entrez_gene_id"]
        == shared_df["GPL6884_entrez_gene_id"]
    )
)

shared_df["entrez_mismatch"] = (
    shared_df["entrez_both_present"]
    & ~shared_df["entrez_match"]
)

shared_df["entrez_missing_either_platform"] = (
    ~shared_df["entrez_both_present"]
)

shared_df = shared_df.sort_values("gene_symbol").reset_index(drop=True)

# ----------------------------------------------------------
# Shared gene list
# ----------------------------------------------------------
shared_gene_list = shared_df[
    [
        "gene_symbol",
        "GPL10558_entrez_gene_id",
        "GPL6884_entrez_gene_id",
        "entrez_match",
    ]
].copy()

shared_gene_list.to_csv(
    INPUT_DIR / "shared_gene_list.tsv",
    sep="\t",
    index=False,
)

# ----------------------------------------------------------
# Platform-specific lists
# ----------------------------------------------------------
pd.DataFrame({"gene_symbol": only_10558}).to_csv(
    TABLE_DIR / "GPL10558_platform_specific_gene_list.tsv",
    sep="\t",
    index=False,
)

pd.DataFrame({"gene_symbol": only_6884}).to_csv(
    TABLE_DIR / "GPL6884_platform_specific_gene_list.tsv",
    sep="\t",
    index=False,
)

shared_df.to_csv(
    TABLE_DIR / "DS001_GSE38900_shared_gene_mapping.tsv",
    sep="\t",
    index=False,
)

# ----------------------------------------------------------
# Entrez consistency
# ----------------------------------------------------------
n_shared = len(shared_df)
n_entrez_both = int(shared_df["entrez_both_present"].sum())
n_entrez_match = int(shared_df["entrez_match"].sum())
n_entrez_mismatch = int(shared_df["entrez_mismatch"].sum())
n_entrez_missing = int(shared_df["entrez_missing_either_platform"].sum())

if n_entrez_both > 0:
    entrez_match_pct = 100 * n_entrez_match / n_entrez_both
else:
    entrez_match_pct = float("nan")

# ----------------------------------------------------------
# Summary table
# ----------------------------------------------------------
summary = pd.DataFrame([
    {
        "GPL10558_gene_count": len(genes1),
        "GPL6884_gene_count": len(genes2),
        "shared_gene_count": len(shared),
        "GPL10558_specific_gene_count": len(only_10558),
        "GPL6884_specific_gene_count": len(only_6884),
        "shared_genes_with_entrez_on_both": n_entrez_both,
        "shared_genes_entrez_match": n_entrez_match,
        "shared_genes_entrez_mismatch": n_entrez_mismatch,
        "shared_genes_missing_entrez_either": n_entrez_missing,
        "entrez_match_percent_among_both_present": entrez_match_pct,
    }
])

summary.to_csv(
    TABLE_DIR / "DS001_GSE38900_cross_platform_harmonization_summary.tsv",
    sep="\t",
    index=False,
)

# ----------------------------------------------------------
# Overlap summary for figure / manuscript
# ----------------------------------------------------------
venn_summary = pd.DataFrame([
    {"category": "GPL10558_only", "gene_count": len(only_10558)},
    {"category": "shared", "gene_count": len(shared)},
    {"category": "GPL6884_only", "gene_count": len(only_6884)},
])

venn_summary.to_csv(
    TABLE_DIR / "DS001_GSE38900_gene_overlap_summary.tsv",
    sep="\t",
    index=False,
)

print()
print("=" * 80)
print("ENTREZ CONSISTENCY")
print("=" * 80)
print(f"Shared genes with Entrez on both platforms: {n_entrez_both:,}")
print(f"Entrez matches: {n_entrez_match:,}")
print(f"Entrez mismatches: {n_entrez_mismatch:,}")
print(f"Missing Entrez on either platform: {n_entrez_missing:,}")
print(f"Entrez match rate among both-present: {entrez_match_pct:.2f}%")

print()
print("=" * 80)
print("CROSS-PLATFORM HARMONIZATION SUMMARY")
print("=" * 80)
print(summary.to_string(index=False))

print("\nCross-platform harmonization completed successfully.")
