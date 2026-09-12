from pathlib import Path
import pandas as pd
import decoupler as dc

ROOT = Path.home() / "Influenza_RSV_Project"

EXPR = ROOT / "data/processed/DS001_GSE38900/analysis_ready/GPL10558_gene_expression.tsv.gz"
NET = ROOT / "results/regulatory_driver_analysis/design/COLLECTRI_HUMAN_FROZEN_v1.0.tsv.gz"
COV = ROOT / "results/regulatory_driver_analysis/qc/GPL10558_COLLECTRI_TARGET_COVERAGE_v1.0.tsv"

OUT = ROOT / "results/regulatory_driver_analysis/replication"
QC = ROOT / "results/regulatory_driver_analysis/qc"
OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

expr = pd.read_csv(EXPR, sep="\t")
net = pd.read_csv(NET, sep="\t")
cov = pd.read_csv(COV, sep="\t")

sample_cols = list(expr.columns[3:])

if expr["gene_symbol"].duplicated().any():
    raise ValueError("Duplicate gene symbols detected in GPL10558 matrix.")

x = (
    expr.set_index("gene_symbol")[sample_cols]
        .apply(pd.to_numeric, errors="raise")
        .T
)

expected_tfs = set(
    cov.loc[cov["measurable_targets"] >= 5, "source"].astype(str)
)

print("Expression matrix for ULM:", x.shape)
print("Expected >=5-target TFs:", len(expected_tfs))

scores, padj = dc.mt.ulm(
    data=x,
    net=net,
    tmin=5,
    verbose=True
)

if list(scores.index) != sample_cols:
    raise ValueError("Returned ULM sample order does not match expression matrix.")

returned_tfs = set(scores.columns.astype(str))
if returned_tfs != expected_tfs:
    missing = sorted(expected_tfs - returned_tfs)
    extra = sorted(returned_tfs - expected_tfs)
    raise ValueError(
        f"TF mismatch: missing={len(missing)}, extra={len(extra)}; "
        f"missing_examples={missing[:10]}, extra_examples={extra[:10]}"
    )

if scores.isna().any().any():
    raise ValueError("Missing values in ULM activity matrix.")
if padj.isna().any().any():
    raise ValueError("Missing values in ULM adjusted-p matrix.")

score_file = OUT / "GPL10558_COLLECTRI_ULM_ACTIVITY_MIN5_v1.0.tsv.gz"
padj_file = OUT / "GPL10558_COLLECTRI_ULM_PADJ_MIN5_v1.0.tsv.gz"
qc_file = QC / "GPL10558_COLLECTRI_ULM_INFERENCE_QC_v1.0.tsv"

scores.to_csv(
    score_file, sep="\t", compression="gzip",
    index=True, index_label="expression_sample_id"
)
padj.to_csv(
    padj_file, sep="\t", compression="gzip",
    index=True, index_label="expression_sample_id"
)

qc = pd.DataFrame([{
    "platform": "GPL10558",
    "input_genes": x.shape[1],
    "input_samples": x.shape[0],
    "ulm_tmin": 5,
    "expected_eligible_TFs": len(expected_tfs),
    "returned_TFs": scores.shape[1],
    "returned_samples": scores.shape[0],
    "missing_activity_values": int(scores.isna().sum().sum()),
    "missing_padj_values": int(padj.isna().sum().sum()),
}])

qc.to_csv(qc_file, sep="\t", index=False)

print("\n=== ULM QC ===")
print(qc.to_string(index=False))
print("\nWritten:")
print(score_file)
print(padj_file)
print(qc_file)
print("\nGPL10558 ULM INFERENCE COMPLETE.")
