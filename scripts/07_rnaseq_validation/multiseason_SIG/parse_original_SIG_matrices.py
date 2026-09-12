import gzip
import shlex
import pandas as pd
from pathlib import Path
import re

BASE = Path("data/rnaseq_validation/influenza_SIG")
OUT = Path("results/rnaseq_validation/influenza_SIG/qc")
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE158592": "GSE158592_norm_SIG_Mar18_070920.txt.gz",
    "GSE155635": "GSE155635_norm_SIG2019_191219.txt.gz",
    "GSE196350": "GSE196350_norm_SIG_2020_190122_2.txt.gz",
    "GSE213168": "GSE213168_sbst1_norm_SIG_2020_010822_2.txt.gz",
}

def parse_shlex(gse, f):
    rows = []

    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            parts = shlex.split(line.strip())

            if i == 0:
                header = parts
            else:
                rows.append(parts)

    row_lengths = sorted(set(len(x) for x in rows))

    print(f"\n{gse}: header={len(header)} row_lengths={row_lengths}")

    if len(row_lengths) != 1:
        raise RuntimeError(f"{gse}: inconsistent row lengths")

    nrow = row_lengths[0]

    if nrow == len(header) + 1:
        header = ["row_id"] + header

    elif nrow != len(header):
        raise RuntimeError(
            f"{gse}: cannot reconcile header {len(header)} vs row {nrow}"
        )

    d = pd.DataFrame(rows, columns=header)
    return d


for gse, filename in FILES.items():

    f = BASE/gse/"processed"/filename

    print("\n====================================")
    print(gse)
    print("====================================")

    if gse == "GSE213168":

        d = pd.read_csv(
            f,
            sep="\t",
            compression="gzip",
            quotechar='"'
        )

        if str(d.columns[0]).startswith("Unnamed"):
            d = d.rename(columns={d.columns[0]: "row_id"})

    else:
        d = parse_shlex(gse, f)

    # ------------------------------------------
    # Identify sample columns
    # ------------------------------------------

    sample_cols = [
        c for c in d.columns
        if re.fullmatch(r"SIG[_-]?\d+", str(c))
    ]

    print("Rows:", d.shape[0])
    print("Columns:", d.shape[1])
    print("Sample columns:", len(sample_cols))

    print("\nFirst columns:")
    print(list(d.columns[:10]))

    print("\nLast sample columns:")
    print(sample_cols[-5:])

    # ------------------------------------------
    # Numeric expression audit
    # ------------------------------------------

    expr = d[sample_cols].apply(
        pd.to_numeric,
        errors="coerce"
    )

    print(
        "Expression NA after numeric conversion:",
        int(expr.isna().sum().sum())
    )

    print(
        "Expression range:",
        float(expr.min().min()),
        "to",
        float(expr.max().max())
    )

    # ------------------------------------------
    # Gene audit
    # ------------------------------------------

    gene_col = "gene_symbol"

    if gene_col not in d.columns:
        raise RuntimeError(f"{gse}: gene_symbol missing")

    print(
        "Genes:",
        len(d),
        "unique gene symbols:",
        d[gene_col].nunique()
    )

    # ------------------------------------------
    # Save clean version
    # ------------------------------------------

    outfile = (
        BASE/gse/"processed"/
        f"{gse}_PARSED_normalized_expression.tsv.gz"
    )

    d.to_csv(
        outfile,
        sep="\t",
        index=False,
        compression="gzip"
    )

    # QC summary
    pd.DataFrame([{
        "GSE": gse,
        "genes": len(d),
        "columns_total": d.shape[1],
        "sample_columns": len(sample_cols),
        "unique_gene_symbols": d["gene_symbol"].nunique(),
        "expression_NA": int(expr.isna().sum().sum()),
        "expression_min": float(expr.min().min()),
        "expression_max": float(expr.max().max()),
    }]).to_csv(
        OUT/f"{gse}_parsed_matrix_QC.tsv",
        sep="\t",
        index=False
    )

    print("Saved:", outfile)

print("\nAll matrices parsed successfully.")
