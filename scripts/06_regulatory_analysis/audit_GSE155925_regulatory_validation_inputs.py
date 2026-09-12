from pathlib import Path
import pandas as pd

ROOT = Path.home() / "Influenza_RSV_Project"

expr_file = ROOT / "results/rnaseq_validation/GSE155925_DE/GSE155925_VST_expression.tsv"
meta_candidates = [
    ROOT / "results/rnaseq_validation/GSE155925_DE/GSE155925_DESeq2_primary_metadata.tsv",
    ROOT / "results/rnaseq_validation/GSE155925_DE/GSE155925_primary_metadata.tsv",
    ROOT / "results/rnaseq_validation/GSE155925_DE/GSE155925_metadata.tsv",
]

print("=== GSE155925 VST EXPRESSION ===")
expr = pd.read_csv(expr_file, sep="\t", nrows=5)
print("Columns:")
print(list(expr.columns))
print("\nFirst rows:")
print(expr.head().to_string(index=False))

print("\n=== POSSIBLE METADATA FILES ===")
for p in meta_candidates:
    print(p.relative_to(ROOT), p.exists())

print("\n=== GSE155925 RESULT FILES ===")
for p in sorted((ROOT / "results/rnaseq_validation/GSE155925_DE").glob("*")):
    print(p.name)

print("\n=== DRIVER TABLE ===")
drivers = pd.read_csv(
    ROOT / "results/regulatory_driver_analysis/tables/REGULATORY_HIGH_CONFIDENCE_DRIVERS_v1.0.tsv",
    sep="\t"
)
print("Frozen high-confidence drivers:", len(drivers))
print(drivers[["TF","regulatory_class","rsv_delta","RSV_replication_tier","frozen170_targets"]].head(50).to_string(index=False))
