import gzip
import csv
import os
from collections import defaultdict

INPUT = (
    "data/rnaseq_validation/GSE155925/metadata/"
    "GSE155925_series_matrix.txt.gz"
)

OUTPUT = (
    "data/rnaseq_validation/GSE155925/metadata/"
    "GSE155925_sample_metadata_raw.tsv"
)

wanted = [
    "!Sample_geo_accession",
    "!Sample_title",
    "!Sample_source_name_ch1",
    "!Sample_characteristics_ch1",
    "!Sample_description",
    "!Sample_relation",
]

rows = defaultdict(list)

with gzip.open(INPUT, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if not line.startswith("!Sample_"):
            continue

        parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
        key = parts[0]

        if key in wanted:
            vals = [x.strip('"') for x in parts[1:]]
            rows[key].append(vals)

# Determine number of samples from GEO accession row
geo_rows = rows.get("!Sample_geo_accession", [])

if not geo_rows:
    raise RuntimeError("No !Sample_geo_accession line found.")

geo = geo_rows[0]
n = len(geo)

print("Number of GEO samples:", n)

records = []

for i in range(n):
    record = {
        "sample_index": i + 1,
        "count_matrix_label": f"Case {i+1}",
        "geo_accession": geo[i],
    }

    # Single-row fields
    for field, clean_name in [
        ("!Sample_title", "title"),
        ("!Sample_source_name_ch1", "source"),
        ("!Sample_description", "description"),
    ]:
        vals = rows.get(field, [])
        if vals and i < len(vals[0]):
            record[clean_name] = vals[0][i]
        else:
            record[clean_name] = ""

    # Multi-row characteristics
    chars = []
    for vals in rows.get("!Sample_characteristics_ch1", []):
        if i < len(vals):
            chars.append(vals[i])

    record["characteristics"] = " | ".join(chars)

    # Relations, often BioSample/SRA
    rels = []
    for vals in rows.get("!Sample_relation", []):
        if i < len(vals):
            rels.append(vals[i])

    record["relations"] = " | ".join(rels)

    records.append(record)

os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)

cols = [
    "sample_index",
    "count_matrix_label",
    "geo_accession",
    "title",
    "source",
    "characteristics",
    "description",
    "relations",
]

with open(OUTPUT, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
    writer.writeheader()
    writer.writerows(records)

print("Written:", OUTPUT)
print("\nFirst 10 samples:")
for r in records[:10]:
    print(
        r["count_matrix_label"],
        r["geo_accession"],
        r["title"],
        r["characteristics"],
        sep="\t"
    )
