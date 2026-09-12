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

DESIGNDIR = Path(
    "results/scrnaseq_validation/GSE283746/design"
)

OUTDIR = Path(
    "results/scrnaseq_validation/GSE283746/pseudobulk"
)

COUNT_OUT = OUTDIR / "counts"
QC_OUT = OUTDIR / "qc"

COUNT_OUT.mkdir(parents=True, exist_ok=True)
QC_OUT.mkdir(parents=True, exist_ok=True)

ELIGIBILITY = DESIGNDIR / (
    "GSE283746_primary_celltype_eligibility_FROZEN_v1.0.tsv"
)

MIN_CELLS = 50

# ============================================================
# 1. Load frozen cell-type specification
# ============================================================

elig = pd.read_csv(ELIGIBILITY, sep="\t")

if "GroupedAnnotation" not in elig.columns:
    # handles index exported as first unnamed column
    first = elig.columns[0]
    elig = elig.rename(columns={first: "GroupedAnnotation"})

primary_types = elig.loc[
    elig["primary_eligible"].astype(str).str.lower().eq("true"),
    "GroupedAnnotation"
].tolist()

print("=== FROZEN PRIMARY CELL TYPES ===")
print("Count:", len(primary_types))
for x in primary_types:
    print(x)

if len(primary_types) != 11:
    raise RuntimeError(
        f"Expected 11 frozen primary cell types, found {len(primary_types)}"
    )

# ============================================================
# 2. Load H5AD annotations
# ============================================================

a = ad.read_h5ad(H5AD, backed="r")
obs = a.obs.copy()

target = obs[
    obs["Combined Condition"].isin(["RSV", "Healthy"])
].copy()

samples = sorted(target["Sample"].astype(str).unique())

print("\n=== TARGET COHORT ===")
print("Participants:", len(samples))
print(
    target[["Sample", "Combined Condition"]]
    .drop_duplicates()["Combined Condition"]
    .value_counts()
    .to_string()
)

if len(samples) != 36:
    raise RuntimeError(
        f"Expected 36 target participants, found {len(samples)}"
    )

# Participant-level metadata
meta_cols = [
    "Sample",
    "Combined Condition",
    "Condition",
    "Age",
    "Sex",
    "GA",
    "Race/Ethnicity",
    "LOS",
    "Steroids",
    "Oxygen",
    "Days Sx",
    "CDSS",
    "Prematurity"
]

meta_cols = [x for x in meta_cols if x in target.columns]

participant_meta = (
    target[meta_cols]
    .drop_duplicates()
    .sort_values("Sample")
    .reset_index(drop=True)
)

# ============================================================
# 3. Feature reference
# ============================================================

reference_symbols = None
reference_feature_ids = None
reference_feature_file = None

# Storage:
# one dictionary per cell type:
# sample -> gene-count vector
pb = {ct: {} for ct in primary_types}

qc_rows = []
mapping_rows = []
feature_audit_rows = []

# ============================================================
# 4. Process each participant sequentially
# ============================================================

for n, sample in enumerate(samples, start=1):

    print("\n========================================")
    print(f"{n}/{len(samples)} {sample}")
    print("========================================")

    sub = target[
        target["Sample"].astype(str).eq(sample)
    ].copy()

    condition = str(sub["Combined Condition"].iloc[0])

    matrix_files = list(
        COUNTDIR.glob(
            f"*_{sample}_scRNA_filtered_matrix.mtx.gz"
        )
    )

    barcode_files = list(
        COUNTDIR.glob(
            f"*_{sample}_scRNA_filtered_barcodes.tsv.gz"
        )
    )

    feature_files = list(
        COUNTDIR.glob(
            f"*_{sample}_scRNA_filtered_features.tsv.gz"
        )
    )

    if not (
        len(matrix_files) == 1 and
        len(barcode_files) == 1 and
        len(feature_files) == 1
    ):
        raise RuntimeError(
            f"{sample}: expected one matrix/barcode/features triplet; "
            f"found matrix={len(matrix_files)}, "
            f"barcodes={len(barcode_files)}, "
            f"features={len(feature_files)}"
        )

    matrix_file = matrix_files[0]
    barcode_file = barcode_files[0]
    feature_file = feature_files[0]

    # --------------------------------------------------------
    # Features
    # --------------------------------------------------------

    features = pd.read_csv(
        feature_file,
        sep="\t",
        header=None,
        compression="gzip",
        dtype=str
    )

    if features.shape[1] < 2:
        raise RuntimeError(
            f"{sample}: feature file has <2 columns"
        )

    feature_ids = features.iloc[:, 0].astype(str).to_numpy()
    gene_symbols = features.iloc[:, 1].astype(str).to_numpy()

    if reference_symbols is None:
        reference_symbols = gene_symbols.copy()
        reference_feature_ids = feature_ids.copy()
        reference_feature_file = str(feature_file)

        print("Reference features:", len(reference_symbols))

    identical = (
        len(gene_symbols) == len(reference_symbols)
        and np.array_equal(gene_symbols, reference_symbols)
        and np.array_equal(feature_ids, reference_feature_ids)
    )

    feature_audit_rows.append({
        "Sample": sample,
        "n_features": len(gene_symbols),
        "identical_to_reference": identical,
        "feature_file": str(feature_file)
    })

    if not identical:
        raise RuntimeError(
            f"{sample}: feature definitions differ from reference"
        )

    # --------------------------------------------------------
    # Barcodes
    # --------------------------------------------------------

    with gzip.open(barcode_file, "rt") as fh:
        barcodes = [
            x.strip()
            for x in fh
            if x.strip()
        ]

    barcode_to_col = {
        bc: i for i, bc in enumerate(barcodes)
    }

    # H5AD index -> original 10x barcode
    h5_names = sub.index.astype(str)

    clean_barcodes = [
        re.sub(
            rf"-{re.escape(sample)}$",
            "",
            x
        )
        for x in h5_names
    ]

    missing = [
        bc for bc in clean_barcodes
        if bc not in barcode_to_col
    ]

    mapping_rows.append({
        "Sample": sample,
        "Condition": condition,
        "H5AD_final_cells": len(clean_barcodes),
        "filtered_barcodes": len(barcodes),
        "mapped_cells": len(clean_barcodes) - len(missing),
        "missing_cells": len(missing),
        "recovery_fraction":
            (len(clean_barcodes) - len(missing)) /
            len(clean_barcodes)
    })

    if missing:
        raise RuntimeError(
            f"{sample}: {len(missing)} final H5AD cells "
            f"missing from filtered matrix"
        )

    # --------------------------------------------------------
    # Count matrix
    # --------------------------------------------------------

    m = mmread(matrix_file)

    if not sparse.issparse(m):
        m = sparse.coo_matrix(m)

    m = m.tocsr()

    expected_shape = (
        len(reference_symbols),
        len(barcodes)
    )

    if m.shape != expected_shape:
        raise RuntimeError(
            f"{sample}: matrix shape {m.shape}; "
            f"expected {expected_shape}"
        )

    # Cell annotations in exact same order as clean_barcodes
    ann = sub["GroupedAnnotation"].astype(str).to_numpy()

    # --------------------------------------------------------
    # Aggregate by frozen primary cell type
    # --------------------------------------------------------

    for ct in primary_types:

        idx_local = np.where(ann == ct)[0]
        n_cells = len(idx_local)

        eligible = n_cells >= MIN_CELLS

        if n_cells > 0:

            matrix_cols = np.array(
                [
                    barcode_to_col[clean_barcodes[i]]
                    for i in idx_local
                ],
                dtype=int
            )

            # genes x selected cells -> summed genes
            summed = np.asarray(
                m[:, matrix_cols].sum(axis=1)
            ).ravel()

            total_umis = int(summed.sum())
            detected_genes = int(np.sum(summed > 0))

        else:
            summed = None
            total_umis = 0
            detected_genes = 0

        qc_rows.append({
            "Sample": sample,
            "Condition": condition,
            "GroupedAnnotation": ct,
            "n_cells": n_cells,
            "eligible_50cells": eligible,
            "total_UMIs": total_umis,
            "detected_genes": detected_genes
        })

        if eligible:
            if not np.allclose(
                summed,
                np.round(summed),
                rtol=0,
                atol=0
            ):
                raise RuntimeError(
                    f"{sample} / {ct}: non-integer pseudobulk counts"
                )

            pb[ct][sample] = summed.astype(np.int64)

    del m

# ============================================================
# 5. Global mapping/feature QC
# ============================================================

mapping_df = pd.DataFrame(mapping_rows)
feature_df = pd.DataFrame(feature_audit_rows)
qc = pd.DataFrame(qc_rows)

mapping_df.to_csv(
    QC_OUT / "GSE283746_pseudobulk_barcode_mapping.tsv",
    sep="\t",
    index=False
)

feature_df.to_csv(
    QC_OUT / "GSE283746_feature_identity_audit.tsv",
    sep="\t",
    index=False
)

qc.to_csv(
    QC_OUT / "GSE283746_participant_celltype_pseudobulk_QC.tsv",
    sep="\t",
    index=False
)

print("\n========================================")
print("BARCODE MAPPING SUMMARY")
print("========================================")

print(
    mapping_df.groupby("Condition")
    .agg(
        participants=("Sample", "size"),
        final_cells=("H5AD_final_cells", "sum"),
        mapped=("mapped_cells", "sum"),
        missing=("missing_cells", "sum")
    )
    .to_string()
)

print("\nFeature-identical participants:")
print(
    int(feature_df["identical_to_reference"].sum()),
    "/",
    len(feature_df)
)

# ============================================================
# 6. Write feature reference
# ============================================================

feature_reference = pd.DataFrame({
    "feature_id": reference_feature_ids,
    "gene_symbol": reference_symbols
})

feature_reference.to_csv(
    OUTDIR / "GSE283746_feature_reference.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)

# ============================================================
# 7. Write one pseudobulk matrix per primary cell type
# ============================================================

matrix_summary = []

for ct in primary_types:

    sample_vectors = pb[ct]

    ordered_samples = sorted(sample_vectors.keys())

    if len(ordered_samples) == 0:
        raise RuntimeError(
            f"{ct}: no eligible pseudobulk samples"
        )

    mat = np.column_stack(
        [sample_vectors[s] for s in ordered_samples]
    )

    out = pd.DataFrame(
        mat,
        columns=ordered_samples
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

    safe = (
        ct.replace("+", "plus")
          .replace("/", "_")
          .replace(" ", "_")
    )

    outfile = (
        COUNT_OUT /
        f"GSE283746_pseudobulk_{safe}_counts.tsv.gz"
    )

    out.to_csv(
        outfile,
        sep="\t",
        index=False,
        compression="gzip"
    )

    subqc = qc[
        (qc["GroupedAnnotation"] == ct) &
        (qc["eligible_50cells"])
    ]

    n_rsv = int(
        (subqc["Condition"] == "RSV").sum()
    )

    n_healthy = int(
        (subqc["Condition"] == "Healthy").sum()
    )

    matrix_summary.append({
        "GroupedAnnotation": ct,
        "samples_total": len(ordered_samples),
        "RSV_samples": n_rsv,
        "Healthy_samples": n_healthy,
        "genes": mat.shape[0],
        "total_UMIs": int(mat.sum()),
        "minimum_cells_per_pseudobulk":
            int(subqc["n_cells"].min()),
        "median_cells_per_pseudobulk":
            float(subqc["n_cells"].median()),
        "minimum_UMIs_per_pseudobulk":
            int(subqc["total_UMIs"].min()),
        "median_UMIs_per_pseudobulk":
            float(subqc["total_UMIs"].median()),
        "output_file": str(outfile)
    })

summary_df = pd.DataFrame(matrix_summary)

summary_df.to_csv(
    OUTDIR /
    "GSE283746_primary_pseudobulk_matrix_summary_FROZEN_v1.0.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 8. Eligible sample metadata
# ============================================================

meta_rows = []

for ct in primary_types:

    eligible_samples = set(
        qc.loc[
            (qc["GroupedAnnotation"] == ct) &
            (qc["eligible_50cells"]),
            "Sample"
        ].astype(str)
    )

    mm = participant_meta[
        participant_meta["Sample"]
        .astype(str)
        .isin(eligible_samples)
    ].copy()

    mm["GroupedAnnotation"] = ct

    meta_rows.append(mm)

pb_meta = pd.concat(
    meta_rows,
    ignore_index=True
)

pb_meta.to_csv(
    OUTDIR /
    "GSE283746_primary_pseudobulk_sample_metadata_FROZEN_v1.0.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 9. Print final QC
# ============================================================

print("\n========================================")
print("PRIMARY PSEUDOBULK MATRICES")
print("========================================")

print(
    summary_df[
        [
            "GroupedAnnotation",
            "RSV_samples",
            "Healthy_samples",
            "genes",
            "minimum_cells_per_pseudobulk",
            "median_cells_per_pseudobulk",
            "minimum_UMIs_per_pseudobulk",
            "median_UMIs_per_pseudobulk"
        ]
    ].to_string(index=False)
)

print("\nTotal pseudobulk profiles:")
print(summary_df["samples_total"].sum())

print("\nMinimum cells among retained pseudobulks:")
print(
    qc.loc[
        qc["eligible_50cells"],
        "n_cells"
    ].min()
)

print("\nAny non-eligible profiles written:")
expected = int(
    qc["eligible_50cells"].sum()
)
actual = int(
    summary_df["samples_total"].sum()
)
print(actual != expected)

print("Expected retained profiles:", expected)
print("Actual retained profiles:", actual)

a.file.close()

print("\nPseudobulk construction complete.")
