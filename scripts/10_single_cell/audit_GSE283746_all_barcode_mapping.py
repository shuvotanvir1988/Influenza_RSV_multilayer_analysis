import anndata as ad
import pandas as pd
import gzip
import re
from pathlib import Path

H5AD = Path(
    "data/scrnaseq_validation/GSE283746/processed/"
    "GSE283746_COVIDRSV_scRNA.h5ad"
)

COUNTDIR = Path(
    "data/scrnaseq_validation/GSE283746/filtered_counts"
)

OUTDIR = Path(
    "results/scrnaseq_validation/GSE283746/design"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

a = ad.read_h5ad(H5AD, backed="r")
obs = a.obs.copy()

target = obs[
    obs["Combined Condition"].isin(["RSV", "Healthy"])
].copy()

samples = sorted(target["Sample"].astype(str).unique())

rows = []

print("=== ALL-PARTICIPANT BARCODE MAPPING AUDIT ===")
print("Participants:", len(samples))

for sample in samples:

    sub = target[
        target["Sample"].astype(str).eq(sample)
    ]

    h5_cells = sub.index.astype(str)

    # remove terminal -<sample>
    cleaned = {
        re.sub(rf"-{re.escape(sample)}$", "", x)
        for x in h5_cells
    }

    matches = list(
        COUNTDIR.glob(
            f"*_{sample}_scRNA_filtered_barcodes.tsv.gz"
        )
    )

    if len(matches) != 1:
        raise RuntimeError(
            f"{sample}: expected exactly one barcode file, "
            f"found {len(matches)}"
        )

    barcode_file = matches[0]

    with gzip.open(barcode_file, "rt") as fh:
        barcodes = {
            x.strip()
            for x in fh
            if x.strip()
        }

    mapped = cleaned & barcodes
    missing = cleaned - barcodes

    recovery = (
        len(mapped) / len(cleaned)
        if len(cleaned)
        else float("nan")
    )

    condition = (
        sub["Combined Condition"]
        .astype(str)
        .iloc[0]
    )

    rows.append({
        "Sample": sample,
        "Condition": condition,
        "H5AD_final_cells": len(cleaned),
        "filtered_barcodes": len(barcodes),
        "mapped_cells": len(mapped),
        "missing_cells": len(missing),
        "recovery_fraction": recovery
    })

    print(
        f"{sample:8s} "
        f"{condition:7s} "
        f"H5AD={len(cleaned):5d} "
        f"filtered={len(barcodes):5d} "
        f"mapped={len(mapped):5d} "
        f"missing={len(missing):3d} "
        f"recovery={recovery:.4f}"
    )

out = pd.DataFrame(rows)

print("\n=== SUMMARY ===")
print("Participants:", len(out))
print("Complete recovery:", int((out["recovery_fraction"] == 1).sum()))
print("Any missing cells:", int((out["missing_cells"] > 0).sum()))
print("Minimum recovery:", out["recovery_fraction"].min())

print("\nBy condition:")
print(
    out.groupby("Condition")
    .agg(
        participants=("Sample", "size"),
        H5AD_cells=("H5AD_final_cells", "sum"),
        mapped_cells=("mapped_cells", "sum"),
        missing_cells=("missing_cells", "sum")
    )
    .to_string()
)

out.to_csv(
    OUTDIR /
    "GSE283746_all_participant_barcode_mapping_FROZEN_v1.0.tsv",
    sep="\t",
    index=False
)

a.file.close()

print("\nAudit complete.")
