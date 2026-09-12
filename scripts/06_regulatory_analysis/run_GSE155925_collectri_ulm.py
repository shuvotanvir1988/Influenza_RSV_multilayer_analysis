from pathlib import Path
import pandas as pd
import decoupler as dc

ROOT = Path.home() / "Influenza_RSV_Project"

EXPR = ROOT / "results/rnaseq_validation/GSE155925_DE/GSE155925_VST_expression.tsv"
NET = ROOT / "results/regulatory_driver_analysis/design/COLLECTRI_HUMAN_FROZEN_v1.0.tsv.gz"
META = ROOT / "results/rnaseq_validation/tables/GSE155925_sample_eligibility_FROZEN_v1.0.tsv"

OUT = ROOT / "results/regulatory_driver_analysis/rnaseq_validation"
QC = ROOT / "results/regulatory_driver_analysis/qc"
OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

expr = pd.read_csv(EXPR, sep="\t")
net = pd.read_csv(NET, sep="\t")
meta = pd.read_csv(META, sep="\t")

if "gene_id" not in expr.columns:
    raise ValueError("Missing gene_id column in GSE155925 VST matrix.")

# Parse SYMBOL:ENSEMBL, preserving the part before the first colon.
expr["gene_symbol"] = expr["gene_id"].astype(str).str.split(":", n=1).str[0].str.strip()

# Remove missing/empty symbols and collapse duplicate symbols deterministically by mean.
expr = expr[expr["gene_symbol"].notna() & (expr["gene_symbol"] != "")].copy()

sample_cols = list(expr.columns[1:-1])

num = expr[sample_cols].apply(pd.to_numeric, errors="raise")
num["gene_symbol"] = expr["gene_symbol"].values
collapsed = num.groupby("gene_symbol", sort=True)[sample_cols].mean()

eligible = meta.loc[meta["primary_eligible"] == True].copy()

eligible_set = set(eligible["count_matrix_label"].astype(str))
vst_set = set(sample_cols)

if not eligible_set.issubset(vst_set):
    missing = sorted(eligible_set - vst_set)
    raise ValueError(f"Eligible samples missing from VST matrix: {missing}")

ordered_samples = [s for s in sample_cols if s in eligible_set]

x = collapsed[ordered_samples].T

print("GSE155925 ULM matrix:", x.shape)

scores, padj = dc.mt.ulm(
    data=x,
    net=net,
    tmin=5,
    verbose=True
)

if list(scores.index) != ordered_samples:
    raise ValueError("ULM sample order mismatch.")

if scores.isna().any().any():
    raise ValueError("Missing activity values in GSE155925 ULM matrix.")

score_file = OUT / "GSE155925_COLLECTRI_ULM_ACTIVITY_MIN5_v1.0.tsv.gz"
padj_file = OUT / "GSE155925_COLLECTRI_ULM_PADJ_MIN5_v1.0.tsv.gz"
qc_file = QC / "GSE155925_COLLECTRI_ULM_INFERENCE_QC_v1.0.tsv"

scores.to_csv(
    score_file, sep="\t", compression="gzip",
    index=True, index_label="count_matrix_label"
)
padj.to_csv(
    padj_file, sep="\t", compression="gzip",
    index=True, index_label="count_matrix_label"
)

qc = pd.DataFrame([{
    "dataset": "GSE155925",
    "eligible_samples": len(ordered_samples),
    "unique_gene_symbols": collapsed.shape[0],
    "ulm_tmin": 5,
    "returned_TFs": scores.shape[1],
    "missing_activity_values": int(scores.isna().sum().sum()),
    "missing_padj_values": int(padj.isna().sum().sum()),
}])

qc.to_csv(qc_file, sep="\t", index=False)

print("\n=== GSE155925 ULM QC ===")
print(qc.to_string(index=False))
print("\nWritten:")
print(score_file)
print(padj_file)
print(qc_file)
