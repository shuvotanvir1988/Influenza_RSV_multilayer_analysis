import pandas as pd
from pathlib import Path

DOMAIN_FILE = Path(
    "results/DS001_GSE38900/pathway_analysis/tables/"
    "preregistered_domains/"
    "DS001_GSE38900_preregistered_domain_pathways.tsv.gz"
)

GENESET_FILES = [
    Path(
        "results/DS001_GSE38900/pathway_analysis/gene_sets/"
        "MSigDB_Hallmark_Homo_sapiens.tsv.gz"
    ),
    Path(
        "results/DS001_GSE38900/pathway_analysis/gene_sets/"
        "MSigDB_Reactome_Homo_sapiens.tsv.gz"
    ),
    Path(
        "results/DS001_GSE38900/pathway_analysis/gene_sets/"
        "MSigDB_GO_BP_Homo_sapiens.tsv.gz"
    ),
]

OUTDIR = Path(
    "results/rnaseq_validation/validation/pathway_validation"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Frozen pathway definitions
# ------------------------------------------------------------

d = pd.read_csv(DOMAIN_FILE, sep="\t")

keep = [
    "pathway",
    "preregistered_domain",
    "biological_domain",
    "NES_RSV_GPL10558",
    "padj_RSV_GPL10558",
    "NES_RSV_GPL6884",
    "padj_RSV_GPL6884",
    "RSV_replicated",
]

frozen = d[keep].drop_duplicates().copy()

print("=== FROZEN PATHWAY TABLE ===")
print("Rows:", len(frozen))
print("Unique pathways:", frozen["pathway"].nunique())

print("\nDomains:")
print(
    frozen.groupby("preregistered_domain")["pathway"]
    .nunique()
    .sort_values(ascending=False)
)

# ------------------------------------------------------------
# Load frozen MSigDB files
# ------------------------------------------------------------

sets = []

for f in GENESET_FILES:
    x = pd.read_csv(f, sep="\t")

    x = x[
        ["gs_name", "gene_symbol", "db_version"]
    ].copy()

    x = x.rename(columns={"gs_name": "pathway"})

    sets.append(x)

gs = pd.concat(sets, ignore_index=True)

print("\n=== MSigDB ===")
print("Membership rows:", len(gs))
print("Unique pathways:", gs["pathway"].nunique())
print("Versions:")
print(gs["db_version"].value_counts(dropna=False))

# ------------------------------------------------------------
# Join ONLY preregistered pathways
# ------------------------------------------------------------

membership = frozen.merge(
    gs,
    on="pathway",
    how="left"
)

# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

matched = (
    membership.groupby("pathway")["gene_symbol"]
    .apply(lambda x: x.notna().any())
)

missing = matched[~matched].index.tolist()

print("\n=== MATCH AUDIT ===")
print("Frozen unique pathways:", frozen["pathway"].nunique())
print("Matched pathways:", int(matched.sum()))
print("Missing pathways:", len(missing))

if missing:
    print("\nMISSING:")
    for x in missing:
        print(x)

coverage = (
    membership.dropna(subset=["gene_symbol"])
    .groupby(
        ["preregistered_domain", "pathway"],
        as_index=False
    )
    .agg(
        gene_set_size=("gene_symbol", "nunique"),
        db_version=("db_version", "first")
    )
)

print("\n=== DOMAIN COUNTS ===")
print(
    coverage.groupby("preregistered_domain")["pathway"]
    .nunique()
    .sort_values(ascending=False)
)

print("\n=== GENE-SET SIZE SUMMARY ===")
print(
    coverage["gene_set_size"]
    .describe()
    .to_string()
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

membership.to_csv(
    OUTDIR /
    "GSE155925_preregistered_pathway_membership_FROZEN_v1.0.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)

coverage.to_csv(
    OUTDIR /
    "GSE155925_preregistered_pathway_coverage_FROZEN_v1.0.tsv",
    sep="\t",
    index=False
)

frozen.to_csv(
    OUTDIR /
    "GSE155925_preregistered_pathways_FROZEN_v1.0.tsv",
    sep="\t",
    index=False
)

print("\nWritten:")
print(OUTDIR)
