#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET_ID = "DS001_GSE38900"

RAW_DIR = PROJECT_ROOT / "data" / "raw" / DATASET_ID
MATRIX_DIR = RAW_DIR / "matrix"
INTERIM_DIR = PROJECT_ROOT / "data" / "interim" / DATASET_ID
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / DATASET_ID
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
LOG_DIR = PROJECT_ROOT / "logs"

for directory in [INTERIM_DIR, PROCESSED_DIR, TABLE_DIR, LOG_DIR]:
    directory.mkdir(parents=True, exist_ok=True)


def clean_geo_value(value: str) -> str:
    value = value.strip()
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    return value.strip()


def split_geo_line(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def normalize_column_name(text: str) -> str:
    text = text.strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def parse_characteristic(value: str) -> tuple[str, str]:
    value = clean_geo_value(value)
    if ":" in value:
        key, val = value.split(":", 1)
        return normalize_column_name(key), val.strip()
    return "unparsed_characteristic", value


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    metadata_rows: dict[str, list[str]] = {}
    expression_lines: list[str] = []
    in_expression_table = False

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")

            if line == "!series_matrix_table_begin":
                in_expression_table = True
                continue
            if line == "!series_matrix_table_end":
                in_expression_table = False
                continue
            if in_expression_table:
                expression_lines.append(line)
                continue

            if line.startswith("!Sample_"):
                fields = split_geo_line(line)
                key = fields[0].replace("!Sample_", "")
                values = [clean_geo_value(v) for v in fields[1:]]

                if key in metadata_rows:
                    suffix = 2
                    new_key = f"{key}_{suffix}"
                    while new_key in metadata_rows:
                        suffix += 1
                        new_key = f"{key}_{suffix}"
                    metadata_rows[new_key] = values
                else:
                    metadata_rows[key] = values

    if not expression_lines:
        raise ValueError(f"No expression table found in {path}")

    reader = csv.reader(expression_lines, delimiter="\t", quotechar='"')
    records = list(reader)
    header = [clean_geo_value(x) for x in records[0]]
    data = [[clean_geo_value(x) for x in row] for row in records[1:] if row]

    expression_df = pd.DataFrame(data, columns=header)
    probe_column = expression_df.columns[0]

    for column in expression_df.columns[1:]:
        expression_df[column] = pd.to_numeric(expression_df[column], errors="coerce")

    sample_ids = list(expression_df.columns[1:])
    sample_metadata = pd.DataFrame({"geo_accession": sample_ids})

    for field, values in metadata_rows.items():
        if len(values) == len(sample_ids):
            sample_metadata[field] = values

    characteristics_columns = [
        c for c in sample_metadata.columns if c.startswith("characteristics_ch1")
    ]
    parsed: dict[str, list[str | None]] = defaultdict(
        lambda: [None] * len(sample_metadata)
    )

    unparsed_counter = 0
    for row_index, row in sample_metadata.iterrows():
        for column in characteristics_columns:
            raw_value = row.get(column)
            if raw_value is None or raw_value == "":
                continue

            key, value = parse_characteristic(str(raw_value))
            if key == "unparsed_characteristic":
                unparsed_counter += 1
                key = f"unparsed_characteristic_{unparsed_counter}"

            existing = parsed[key][row_index]
            parsed[key][row_index] = value if existing is None else f"{existing}; {value}"

    for key, values in parsed.items():
        sample_metadata[key] = values

    match = re.search(r"(GPL\d+)", path.name)
    platform_id = match.group(1) if match else "unknown"

    # Preserve GEO-provided columns if they conflict with our internal fields.
    rename_map = {}
    for column in ["dataset_id", "gse_accession", "platform_id", "matrix_file"]:
        if column in sample_metadata.columns:
            rename_map[column] = f"geo_{column}"

    if rename_map:
        sample_metadata = sample_metadata.rename(columns=rename_map)

    sample_metadata.insert(0, "dataset_id", "DS001")
    sample_metadata.insert(1, "gse_accession", "GSE38900")
    sample_metadata.insert(2, "platform_id", platform_id)
    sample_metadata.insert(3, "matrix_file", path.name)

    summary = {
        "dataset_id": "DS001",
        "gse_accession": "GSE38900",
        "platform_id": platform_id,
        "matrix_file": path.name,
        "probe_column": probe_column,
        "number_of_probes": int(expression_df.shape[0]),
        "number_of_samples": int(expression_df.shape[1] - 1),
        "metadata_columns_recovered": int(sample_metadata.shape[1]),
        "expression_missing_values": int(
            expression_df.iloc[:, 1:].isna().sum().sum()
        ),
    }
    return sample_metadata, expression_df, summary


def metadata_completeness(metadata: pd.DataFrame) -> pd.DataFrame:
    results = []
    for column in metadata.columns:
        series = metadata[column].astype("string")
        missing = (
            series.isna()
            | series.str.strip().eq("")
            | series.str.lower().isin(
                ["na", "n/a", "unknown", "not available", "none"]
            )
        )
        results.append({
            "column": column,
            "n_samples": len(metadata),
            "n_missing": int(missing.sum()),
            "percent_missing": round(100 * float(missing.mean()), 2),
            "n_unique_nonmissing": int(series[~missing].nunique(dropna=True)),
        })
    return pd.DataFrame(results).sort_values(
        ["percent_missing", "column"], ascending=[True, True]
    )


def main() -> None:
    matrix_files = sorted(MATRIX_DIR.glob("GSE38900-GPL*_series_matrix.txt.gz"))
    if not matrix_files:
        print(f"ERROR: No Series Matrix files found in {MATRIX_DIR}", file=sys.stderr)
        sys.exit(1)

    all_metadata = []
    summaries = []

    for matrix_file in matrix_files:
        print(f"Parsing: {matrix_file.name}")
        metadata, expression, summary = parse_series_matrix(matrix_file)
        platform_id = summary["platform_id"]

        metadata.to_csv(
            INTERIM_DIR / f"DS001_GSE38900_{platform_id}_sample_metadata_raw.tsv",
            sep="\t", index=False
        )
        expression.to_csv(
            INTERIM_DIR / f"DS001_GSE38900_{platform_id}_expression_matrix.tsv.gz",
            sep="\t", index=False, compression="gzip"
        )

        all_metadata.append(metadata)
        summaries.append(summary)

        print(f"  Platform: {platform_id}")
        print(f"  Samples: {summary['number_of_samples']}")
        print(f"  Probes: {summary['number_of_probes']}")
        print(f"  Missing expression values: {summary['expression_missing_values']}")

    combined = pd.concat(all_metadata, ignore_index=True, sort=False)
    duplicates = combined[
        combined.duplicated(subset=["geo_accession"], keep=False)
    ].copy()

    combined_output = (
        PROCESSED_DIR / "DS001_GSE38900_master_sample_metadata_initial.tsv"
    )
    combined.to_csv(combined_output, sep="\t", index=False)

    summary_df = pd.DataFrame(summaries)
    summary_output = TABLE_DIR / "DS001_GSE38900_matrix_summary.tsv"
    summary_df.to_csv(summary_output, sep="\t", index=False)

    completeness_output = TABLE_DIR / "DS001_GSE38900_metadata_completeness.tsv"
    metadata_completeness(combined).to_csv(
        completeness_output, sep="\t", index=False
    )

    duplicate_output = TABLE_DIR / "DS001_GSE38900_duplicate_sample_check.tsv"
    duplicates.to_csv(duplicate_output, sep="\t", index=False)

    run_summary = {
        "dataset_id": "DS001",
        "gse_accession": "GSE38900",
        "matrix_files_parsed": len(matrix_files),
        "platforms": summary_df["platform_id"].tolist(),
        "total_samples": int(len(combined)),
        "unique_geo_samples": int(combined["geo_accession"].nunique()),
        "duplicate_geo_sample_rows": int(len(duplicates)),
        "combined_metadata_columns": int(combined.shape[1]),
    }

    with (LOG_DIR / "DS001_GSE38900_metadata_extraction_summary.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(run_summary, handle, indent=2)

    print("\nMetadata extraction complete.")
    print(f"Total samples: {len(combined)}")
    print(f"Unique GEO samples: {combined['geo_accession'].nunique()}")
    print(f"Duplicate GEO sample rows: {len(duplicates)}")
    print(f"Master metadata: {combined_output}")
    print(f"Matrix summary: {summary_output}")
    print(f"Completeness report: {completeness_output}")


if __name__ == "__main__":
    main()
