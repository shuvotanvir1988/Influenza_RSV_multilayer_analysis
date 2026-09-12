import anndata as ad
import pandas as pd
import numpy as np
import gzip
import re
from pathlib import Path
from scipy.io import mmread
from scipy import sparse

H5AD = Path(
    "data/scrnaseq_validation/GSE283746/processed/"
    "GSE283746_COVIDRSV_scRNA.h5ad"
)

COUNTDIR = Path(
    "data/scrnaseq_validation/GSE283746/filtered_counts"
)

OUTDIR = Path(
    "results/scrnaseq_validation/GSE283746/isg_state/pseudobulk"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

MIN_CELLS = 25

STATES = {
    "ISG_Naive_CD4_T": "ISG+ Naive CD4+ T",
    "Reference_Naive_CD4_T": "Naive CD4+ T",
}

a = ad.read_h5ad(H5AD, backed="r")
obs = a.obs.copy()

target = obs[
    obs["Combined Condition"].isin(["RSV", "Healthy"])
].copy()

samples = sorted(target["Sample"].astype(str).unique())

pb = {k: {} for k in STATES}
qc_rows = []

reference_feature_ids = None
reference_symbols = None

for sample in samples:

    sub = target[
        target["Sample"].astype(str).eq(sample)
    ].copy()

    condition = str(sub["Combined Condition"].iloc[0])

    matrix_file = list(
        COUNTDIR.glob(
            f"*_{sample}_scRNA_filtered_matrix.mtx.gz"
        )
    )[0]

    barcode_file = list(
        COUNTDIR.glob(
            f"*_{sample}_scRNA_filtered_barcodes.tsv.gz"
        )
    )[0]

    feature_file = list(
        COUNTDIR.glob(
            f"*_{sample}_scRNA_filtered_features.tsv.gz"
        )
    )[0]

    features = pd.read_csv(
        feature_file,
        sep="\t",
        header=None,
        compression="gzip",
        dtype=str
    )

    feature_ids = features.iloc[:,0].astype(str).to_numpy()
    gene_symbols = features.iloc[:,1].astype(str).to_numpy()

    if reference_feature_ids is None:
        reference_feature_ids = feature_ids.copy()
        reference_symbols = gene_symbols.copy()
    else:
        if not (
            np.array_equal(feature_ids, reference_feature_ids)
            and np.array_equal(gene_symbols, reference_symbols)
        ):
            raise RuntimeError(
                f"{sample}: feature definitions differ"
            )

    with gzip.open(barcode_file, "rt") as fh:
        barcodes = [
            x.strip() for x in fh if x.strip()
        ]

    bc_to_col = {
        bc:i for i,bc in enumerate(barcodes)
    }

    h5_names = sub.index.astype(str)

    clean = [
        re.sub(
            rf"-{re.escape(sample)}$",
            "",
            x
        )
        for x in h5_names
    ]

    if any(bc not in bc_to_col for bc in clean):
        raise RuntimeError(
            f"{sample}: missing H5AD cell barcode"
        )

    m = mmread(matrix_file)
    if not sparse.issparse(m):
        m = sparse.coo_matrix(m)
    m = m.tocsr()

    ann = sub["Annotation"].astype(str).to_numpy()

    for state_id, annotation in STATES.items():

        idx = np.where(ann == annotation)[0]
        n_cells = len(idx)
        eligible = n_cells >= MIN_CELLS

        if n_cells > 0:
            matrix_cols = np.array(
                [bc_to_col[clean[i]] for i in idx],
                dtype=int
            )

            summed = np.asarray(
                m[:, matrix_cols].sum(axis=1)
            ).ravel()

            total_umis = int(summed.sum())
            detected = int((summed > 0).sum())

        else:
            summed = None
            total_umis = 0
            detected = 0

        qc_rows.append({
            "Sample": sample,
            "Condition": condition,
            "state_id": state_id,
            "Annotation": annotation,
            "n_cells": n_cells,
            "eligible_25cells": eligible,
            "total_UMIs": total_umis,
            "detected_genes": detected
        })

        if eligible:
            pb[state_id][sample] = summed.astype(np.int64)

qc = pd.DataFrame(qc_rows)

qc.to_csv(
    OUTDIR /
    "GSE283746_NaiveCD4_ISG_state_pseudobulk_QC.tsv",
    sep="\t",
    index=False
)

summary_rows = []

for state_id, sample_vectors in pb.items():

    ordered = sorted(sample_vectors.keys())

    mat = np.column_stack(
        [sample_vectors[s] for s in ordered]
    )

    out = pd.DataFrame(
        mat,
        columns=ordered
    )

    out.insert(
        0,
        "gene_symbol",
        reference_symbols
    )

    out.insert(
        0,
        "feature_id",
        reference_feature_ids
    )

    outfile = (
        OUTDIR /
        f"GSE283746_{state_id}_counts.tsv.gz"
    )

    out.to_csv(
        outfile,
        sep="\t",
        index=False,
        compression="gzip"
    )

    q = qc[
        (qc["state_id"] == state_id) &
        (qc["eligible_25cells"])
    ]

    summary_rows.append({
        "state_id": state_id,
        "total_samples": len(ordered),
        "RSV_samples":
            int((q["Condition"] == "RSV").sum()),
        "Healthy_samples":
            int((q["Condition"] == "Healthy").sum()),
        "minimum_cells":
            int(q["n_cells"].min()),
        "median_cells":
            float(q["n_cells"].median()),
        "minimum_UMIs":
            int(q["total_UMIs"].min()),
        "median_UMIs":
            float(q["total_UMIs"].median()),
        "output_file": str(outfile)
    })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUTDIR /
    "GSE283746_NaiveCD4_ISG_state_pseudobulk_summary_FROZEN_v1.0.tsv",
    sep="\t",
    index=False
)

print("=== STATE PSEUDOBULK SUMMARY ===")
print(summary.to_string(index=False))

print("\n=== PARTICIPANT ELIGIBILITY ===")
print(
    qc.groupby(
        ["state_id","Condition"]
    )["eligible_25cells"]
    .sum()
    .to_string()
)

a.file.close()
