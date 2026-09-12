from pathlib import Path
import csv
import re
import pandas as pd

RAW = Path(
    "data/proteomics/raw/Influenza-P01_Schughart2025/"
    "sbst3_norm_SIG_Somalogic_UTHSC_2021_291122.txt"
)

META = Path(
    "data/proteomics/metadata/Influenza-P01_Schughart2025/"
    "sbst3_target_SIG_proteome_Soma_251122a.xlsx"
)

OUT = Path(
    "results/proteomics_validation/Influenza-P01/analysis_inputs"
)

OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Parse SomaScan file
# ------------------------------------------------------------

rows = []

with open(RAW, "r", newline="", errors="replace") as fh:
    reader = csv.reader(
        fh,
        delimiter=" ",
        quotechar='"',
        skipinitialspace=True
    )

    for row in reader:
        row = [x for x in row if x != ""]
        if row:
            rows.append(row)

header = rows[0]
data = rows[1:]

if len(data[0]) == len(header) + 1:
    header = ["row_index"] + header

if any(len(r) != len(header) for r in data):
    raise SystemExit(
        "STOP: malformed source rows."
    )

df = pd.DataFrame(
    data,
    columns=header
)

sample_cols = [
    c for c in df.columns
    if re.fullmatch(r"SIG_\d+", c)
]

human = df.loc[
    df["Organism"].astype(str).str.strip() == "Human"
].copy()

assert len(human) == 7298
assert len(sample_cols) == 84

annotation_cols = [
    "row_index",
    "anlyt_ID",
    "TargetFullName",
    "Target",
    "UniProt",
    "EntrezGeneID",
    "EntrezGeneSymbol",
    "Organism",
    "SeqId",
    "SeqIdVersion",
]

annotation = human[
    annotation_cols
].copy()

expr = human[
    ["row_index"] + sample_cols
].copy()

for c in sample_cols:
    expr[c] = pd.to_numeric(
        expr[c],
        errors="raise"
    )

assert (
    expr[sample_cols]
    .isna()
    .sum()
    .sum()
    == 0
)

# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------

meta = pd.read_excel(
    META,
    sheet_name="Sheet 1"
).copy()

assert len(meta) == 84
assert meta["sample_ID"].nunique() == 84
assert set(meta["sample_ID"]) == set(sample_cols)

# Reorder metadata exactly to proteome matrix.
meta = (
    meta
    .set_index("sample_ID")
    .loc[sample_cols]
    .reset_index()
)

assert meta["infect_status"].value_counts().to_dict() == {
    "infected": 61,
    "healthy_control": 23,
}

assert meta["sex"].value_counts().to_dict() == {
    "f": 53,
    "m": 31,
}

assert meta["ICU"].value_counts().to_dict() == {
    "n": 58,
    "y": 26,
}

# ------------------------------------------------------------
# Write analysis-ready files
# ------------------------------------------------------------

annotation_file = OUT / (
    "Influenza-P01_human_assay_annotation_v1.0.tsv"
)

expr_file = OUT / (
    "Influenza-P01_human_expression_matrix_v1.0.tsv"
)

meta_file = OUT / (
    "Influenza-P01_sample_metadata_v1.0.tsv"
)

annotation.to_csv(
    annotation_file,
    sep="\t",
    index=False
)

expr.to_csv(
    expr_file,
    sep="\t",
    index=False
)

meta.to_csv(
    meta_file,
    sep="\t",
    index=False
)

print("=== ANALYSIS INPUT QC ===")
print("Human assays:", len(annotation))
print("Samples:", len(sample_cols))

print("\nInfection status:")
print(
    meta["infect_status"]
    .value_counts()
    .to_string()
)

print("\nSex:")
print(
    meta["sex"]
    .value_counts()
    .to_string()
)

print("\nICU:")
print(
    meta["ICU"]
    .value_counts()
    .to_string()
)

print("\nExpression missing values:",
      int(expr[sample_cols].isna().sum().sum()))

print("\nWritten:")
print(annotation_file)
print(expr_file)
print(meta_file)

print("\nSTATUS: ANALYSIS INPUTS PASS")
