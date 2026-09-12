from pathlib import Path
import pandas as pd
import decoupler as dc

ROOT = Path.home() / "Influenza_RSV_Project"

NET = ROOT / "results/regulatory_driver_analysis/design/COLLECTRI_HUMAN_FROZEN_v1.0.tsv.gz"
IN = ROOT / "results/regulatory_driver_analysis/scrnaseq_validation/logCPM"
OUT = ROOT / "results/regulatory_driver_analysis/scrnaseq_validation/ULM"
QC = ROOT / "results/regulatory_driver_analysis/qc"

OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

net = pd.read_csv(NET, sep="\t")

files = sorted(IN.glob("GSE283746_*_TMM_logCPM_for_ULM_v1.0.tsv.gz"))
if len(files) != 11:
    raise ValueError(f"Expected 11 logCPM files, found {len(files)}")

summary = []

for f in files:
    stem = f.name.replace("_TMM_logCPM_for_ULM_v1.0.tsv.gz", "")
    safe = stem.replace("GSE283746_", "")

    print("\n========================================")
    print(safe)
    print("========================================")

    d = pd.read_csv(f, sep="\t")

    if d["gene_symbol"].duplicated().any():
        raise ValueError(f"Duplicate gene symbols remain in {f.name}")

    sample_cols = list(d.columns[1:])
    x = d.set_index("gene_symbol")[sample_cols].T

    scores, padj = dc.mt.ulm(
        data=x,
        net=net,
        tmin=5,
        verbose=False
    )

    score_file = OUT / f"GSE283746_{safe}_COLLECTRI_ULM_ACTIVITY_MIN5_v1.0.tsv.gz"
    padj_file = OUT / f"GSE283746_{safe}_COLLECTRI_ULM_PADJ_MIN5_v1.0.tsv.gz"

    scores.to_csv(
        score_file, sep="\t", compression="gzip",
        index=True, index_label="Sample"
    )
    padj.to_csv(
        padj_file, sep="\t", compression="gzip",
        index=True, index_label="Sample"
    )

    summary.append({
        "cell_type_safe": safe,
        "samples": scores.shape[0],
        "TFs": scores.shape[1],
        "missing_activity_values": int(scores.isna().sum().sum()),
        "missing_padj_values": int(padj.isna().sum().sum()),
        "activity_file": str(score_file.relative_to(ROOT)),
    })

summary_df = pd.DataFrame(summary)

summary_df.to_csv(
    QC / "GSE283746_REGULATORY_ULM_QC_v1.0.tsv",
    sep="\t",
    index=False
)

print("\n=== ULM SUMMARY ===")
print(summary_df.to_string(index=False))
print("\nGSE283746 CELL-TYPE ULM COMPLETE.")
