import gzip
import csv
import re
import pandas as pd
from pathlib import Path

GS = {
    "2018": "GSE158592",
    "2019": "GSE155635",
    "2020": "GSE196350",
    "2022": "GSE213168",
}

BASE = Path("data/rnaseq_validation/influenza_SIG")
OUTDIR = Path("results/rnaseq_validation/influenza_SIG/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)

records = []

for year, gse in GS.items():

    f = BASE/gse/"metadata"/f"{gse}_series_matrix.txt.gz"

    fields = {}
    characteristics = []

    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:

        for line in fh:

            if line.startswith("!Sample_title"):
                fields["title"] = next(
                    csv.reader([line.rstrip("\n")], delimiter="\t")
                )[1:]

            elif line.startswith("!Sample_geo_accession"):
                fields["gsm"] = next(
                    csv.reader([line.rstrip("\n")], delimiter="\t")
                )[1:]

            elif line.startswith("!Sample_description"):
                fields["description"] = next(
                    csv.reader([line.rstrip("\n")], delimiter="\t")
                )[1:]

            elif line.startswith("!Sample_source_name_ch1"):
                fields["source"] = next(
                    csv.reader([line.rstrip("\n")], delimiter="\t")
                )[1:]

            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(
                    next(
                        csv.reader([line.rstrip("\n")], delimiter="\t")
                    )[1:]
                )

    n = len(fields["title"])

    for i in range(n):

        title = fields["title"][i].strip('"')

        m = re.search(r"SIG[\s_-]*(\d+)", title, flags=re.I)

        if not m:
            continue

        num = m.group(1)

        combined_sample = f"SIG_{num}_{year}"

        chars = []
        for row in characteristics:
            if i < len(row):
                chars.append(row[i].strip('"'))

        chartext = " | ".join(chars)

        # --------------------------------------------
        # Disease status
        # --------------------------------------------
        lower = chartext.lower()

        if (
            "healthy_control" in lower
            or "healthy control" in lower
        ):
            status = "HEALTHY_CONTROL"

        elif (
            "status: infected" in lower
            or "disease state: infected" in lower
            or "treatment: infected" in lower
        ):
            status = "INFLUENZA_INFECTED"

        else:
            status = "UNKNOWN"

        # --------------------------------------------
        # Sex/gender
        # --------------------------------------------
        sex = None

        msex = re.search(
            r"(?:sex|gender):\s*([fm]|male|female|NA)",
            chartext,
            flags=re.I
        )

        if msex:
            sx = msex.group(1).lower()

            if sx in ["f", "female"]:
                sex = "Female"
            elif sx in ["m", "male"]:
                sex = "Male"
            else:
                sex = "NA"

        records.append({
            "combined_sample": combined_sample,
            "year": int(year),
            "GSE": gse,
            "GSM": fields["gsm"][i].strip('"'),
            "original_title": title,
            "status": status,
            "sex": sex,
            "characteristics": chartext,
        })

meta = pd.DataFrame(records)

# ------------------------------------------------------------
# Read exact 208 sample columns from parsed matrix
# ------------------------------------------------------------

matrix_file = (
    BASE/"GSE272879"/"processed"/
    "GSE272879_combined_matrix_PARSED.tsv.gz"
)

d = pd.read_csv(
    matrix_file,
    sep="\t",
    nrows=1
)

sample_cols = [
    c for c in d.columns
    if re.fullmatch(
        r"SIG_.+_(2018|2019|2020|2022)",
        str(c)
    )
]

manifest = pd.DataFrame({
    "combined_sample": sample_cols
})

joined = manifest.merge(
    meta,
    on="combined_sample",
    how="left",
    validate="one_to_one"
)

outfile = OUTDIR/"SIG208_sample_manifest_FROZEN_v1.0.tsv"
joined.to_csv(outfile, sep="\t", index=False)

print("=== GSE272879 / SIG208 METADATA AUDIT ===")

print("\nExpression samples:", len(manifest))
print("Metadata matched:", joined["GSM"].notna().sum())
print("Metadata missing:", joined["GSM"].isna().sum())

print("\n=== STATUS COUNTS ===")
print(joined["status"].value_counts(dropna=False))

print("\n=== STATUS BY YEAR ===")
print(pd.crosstab(joined["year"], joined["status"]))

print("\n=== SEX BY STATUS ===")
print(pd.crosstab(joined["status"], joined["sex"]))

print("\n=== SEX BY YEAR / STATUS ===")
print(
    joined.groupby(
        ["year", "status", "sex"],
        dropna=False
    ).size()
)

print("\n=== UNMATCHED SAMPLES ===")
bad = joined[joined["GSM"].isna()]

if len(bad):
    print(bad["combined_sample"].to_string(index=False))
else:
    print("NONE")

print("\n=== UNKNOWN STATUS ===")
unknown = joined[joined["status"] == "UNKNOWN"]

if len(unknown):
    print(
        unknown[
            [
                "combined_sample",
                "GSE",
                "GSM",
                "characteristics"
            ]
        ].to_string(index=False)
    )
else:
    print("NONE")

print("\nWritten:")
print(outfile)
