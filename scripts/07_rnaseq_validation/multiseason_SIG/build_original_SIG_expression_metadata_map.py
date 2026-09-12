import gzip
import csv
import re
import pandas as pd
from pathlib import Path

BASE = Path("data/rnaseq_validation/influenza_SIG")
OUT = Path("results/rnaseq_validation/influenza_SIG/metadata")
OUT.mkdir(parents=True, exist_ok=True)

INFO = {
    "GSE158592": {
        "year": 2018,
        "matrix": "GSE158592_PARSED_normalized_expression.tsv.gz",
    },
    "GSE155635": {
        "year": 2019,
        "matrix": "GSE155635_PARSED_normalized_expression.tsv.gz",
    },
    "GSE196350": {
        "year": 2020,
        "matrix": "GSE196350_PARSED_normalized_expression.tsv.gz",
    },
    "GSE213168": {
        "year": 2022,
        "matrix": "GSE213168_PARSED_normalized_expression.tsv.gz",
    },
}

all_records = []

for gse, info in INFO.items():

    year = info["year"]

    meta_file = (
        BASE/gse/"metadata"/f"{gse}_series_matrix.txt.gz"
    )

    matrix_file = (
        BASE/gse/"processed"/info["matrix"]
    )

    # ------------------------------------------------
    # Parse GEO metadata
    # ------------------------------------------------

    titles = None
    gsms = None
    descriptions = None
    sources = None
    characteristics = []

    with gzip.open(
        meta_file,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as fh:

        for line in fh:

            parts = next(
                csv.reader(
                    [line.rstrip("\n")],
                    delimiter="\t"
                )
            )

            if line.startswith("!Sample_title"):
                titles = [
                    x.strip('"') for x in parts[1:]
                ]

            elif line.startswith("!Sample_geo_accession"):
                gsms = [
                    x.strip('"') for x in parts[1:]
                ]

            elif line.startswith("!Sample_description"):
                descriptions = [
                    x.strip('"') for x in parts[1:]
                ]

            elif line.startswith("!Sample_source_name_ch1"):
                sources = [
                    x.strip('"') for x in parts[1:]
                ]

            elif line.startswith(
                "!Sample_characteristics_ch1"
            ):
                characteristics.append(
                    [
                        x.strip('"')
                        for x in parts[1:]
                    ]
                )

    n = len(gsms)

    metadata_records = []

    for i in range(n):

        title = titles[i]

        m = re.search(
            r"SIG[\s_-]*(\d+)",
            title,
            flags=re.I
        )

        if not m:
            raise RuntimeError(
                f"{gse}: cannot derive SIG ID from {title}"
            )

        sig_num = m.group(1)
        expression_label = f"SIG_{sig_num}"

        chars = [
            row[i]
            for row in characteristics
            if i < len(row)
        ]

        chartext = " | ".join(chars)
        lower = chartext.lower()

        if "healthy_control" in lower:
            status = "HEALTHY_CONTROL"

        elif (
            "status: infected" in lower
            or "disease state: infected" in lower
            or "treatment: infected" in lower
        ):
            status = "INFLUENZA_INFECTED"

        else:
            status = "UNKNOWN"

        sex = "NA"

        sx = re.search(
            r"(?:sex|gender):\s*"
            r"(f|m|female|male|na)",
            chartext,
            flags=re.I
        )

        if sx:
            s = sx.group(1).lower()

            if s in ["f", "female"]:
                sex = "Female"

            elif s in ["m", "male"]:
                sex = "Male"

        metadata_records.append({
            "GSE": gse,
            "year": year,
            "SIG_ID": f"SIG_{sig_num}",
            "expression_label": expression_label,
            "GSM": gsms[i],
            "title": title,
            "description":
                descriptions[i]
                if descriptions else "",
            "source":
                sources[i]
                if sources else "",
            "status": status,
            "sex": sex,
            "characteristics": chartext,
        })

    meta = pd.DataFrame(metadata_records)

    # ------------------------------------------------
    # Read exact expression columns
    # ------------------------------------------------

    header = pd.read_csv(
        matrix_file,
        sep="\t",
        compression="gzip",
        nrows=0
    ).columns.tolist()

    expr_cols = [
        c for c in header
        if re.fullmatch(r"SIG_\d+", str(c))
    ]

    expr = pd.DataFrame({
        "expression_label": expr_cols
    })

    # ------------------------------------------------
    # Join
    # ------------------------------------------------

    joined = expr.merge(
        meta,
        on="expression_label",
        how="left",
        validate="one_to_one"
    )

    joined["expression_present"] = True

    outfile = (
        OUT /
        f"{gse}_expression_metadata_map.tsv"
    )

    joined.to_csv(
        outfile,
        sep="\t",
        index=False
    )

    print("\n================================")
    print(gse, year)
    print("================================")

    print("GEO metadata samples:", len(meta))
    print("Expression samples:", len(expr))

    print(
        "Mapped expression samples:",
        int(joined["GSM"].notna().sum())
    )

    print(
        "Unmapped expression samples:",
        int(joined["GSM"].isna().sum())
    )

    print("\nStatus:")
    print(
        joined["status"]
        .value_counts(dropna=False)
    )

    print("\nSex × status:")
    print(
        pd.crosstab(
            joined["status"],
            joined["sex"]
        )
    )

    # GEO samples absent from expression matrix
    absent = (
        meta.loc[
            ~meta["expression_label"]
            .isin(expr["expression_label"])
        ]
    )

    print(
        "\nGEO samples absent from expression matrix:",
        len(absent)
    )

    if len(absent):
        print(
            absent[
                [
                    "SIG_ID",
                    "GSM",
                    "status",
                    "sex"
                ]
            ].to_string(index=False)
        )

    all_records.append(joined)

# ------------------------------------------------
# Combined cross-season audit
# ------------------------------------------------

combined = pd.concat(
    all_records,
    ignore_index=True
)

combined["participant_id"] = (
    combined["SIG_ID"]
)

combined.to_csv(
    OUT /
    "SIG_original_four_season_manifest.tsv",
    sep="\t",
    index=False
)

print("\n================================")
print("CROSS-SEASON PARTICIPANT AUDIT")
print("================================")

print(
    "Total expression samples:",
    len(combined)
)

print(
    "Unique SIG participant IDs:",
    combined["participant_id"].nunique()
)

counts = (
    combined.groupby("participant_id")
    .size()
    .sort_values(ascending=False)
)

repeated = counts[counts > 1]

print(
    "Participant IDs occurring in >1 season:",
    len(repeated)
)

print("\nRepeated IDs:")
print(repeated.to_string())

print("\nStatus consistency for repeated IDs:")

for pid in repeated.index:

    z = combined[
        combined["participant_id"] == pid
    ]

    print(
        "\n",
        pid,
        z[
            [
                "year",
                "GSE",
                "GSM",
                "status",
                "sex"
            ]
        ].to_dict("records")
    )
