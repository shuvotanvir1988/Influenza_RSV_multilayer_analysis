#!/usr/bin/env python3

from pathlib import Path
import json
import sys

import pandas as pd


ROOT = Path.home() / "Influenza_RSV_Project"
DATASET = "DS001_GSE38900"

INTERIM = ROOT / "data" / "interim" / DATASET
PROCESSED = ROOT / "data" / "processed" / DATASET
RESULTS = ROOT / "results" / DATASET / "preprocessing" / "tables"
LOGS = ROOT / "logs" / "session14"

PLATFORMS = {
    "GPL10558": {
        "expression": INTERIM / f"{DATASET}_GPL10558_expression_matrix.tsv.gz",
        "metadata": INTERIM / f"{DATASET}_GPL10558_sample_metadata_raw.tsv",
        "expected_features": 47323,
        "expected_samples": 36,
    },
    "GPL6884": {
        "expression": INTERIM / f"{DATASET}_GPL6884_expression_matrix.tsv.gz",
        "metadata": INTERIM / f"{DATASET}_GPL6884_sample_metadata_raw.tsv",
        "expected_features": 48803,
        "expected_samples": 205,
    },
}

HARMONIZED = (
    PROCESSED / f"{DATASET}_master_sample_metadata_harmonized.tsv"
)


def identify_sample_column(df: pd.DataFrame) -> str:
    preferred = [
        "sample_id",
        "gsm",
        "gsm_id",
        "geo_accession",
        "sample_accession",
        "accession",
    ]

    lowered = {
        str(column).strip().lower(): column
        for column in df.columns
    }

    for name in preferred:
        if name in lowered:
            return lowered[name]

    for column in df.columns:
        values = df[column].dropna().astype(str).str.strip()
        if len(values) > 0 and values.str.match(r"^GSM\d+$").mean() > 0.8:
            return column

    raise RuntimeError(
        "Could not identify a GSM/sample-ID column. "
        f"Columns found: {list(df.columns)}"
    )


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    if not HARMONIZED.exists():
        raise FileNotFoundError(f"Missing harmonized metadata: {HARMONIZED}")

    harmonized = pd.read_csv(
        HARMONIZED,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    harmonized_id_col = identify_sample_column(harmonized)
    harmonized_ids = set(
        harmonized[harmonized_id_col]
        .dropna()
        .astype(str)
        .str.strip()
    )

    rows = []
    platform_sample_sets = {}

    print(f"Project root: {ROOT}")
    print(f"Harmonized metadata rows: {len(harmonized)}")
    print(f"Harmonized sample column: {harmonized_id_col}")
    print()

    for platform, files in PLATFORMS.items():
        expression_path = files["expression"]
        metadata_path = files["metadata"]

        if not expression_path.exists():
            raise FileNotFoundError(expression_path)

        if not metadata_path.exists():
            raise FileNotFoundError(metadata_path)

        expression = pd.read_csv(
            expression_path,
            sep="\t",
            compression="gzip",
            index_col=0,
            low_memory=False,
        )

        metadata = pd.read_csv(
            metadata_path,
            sep="\t",
            dtype=str,
            low_memory=False,
        )

        metadata_id_col = identify_sample_column(metadata)
        metadata_ids = set(
            metadata[metadata_id_col]
            .dropna()
            .astype(str)
            .str.strip()
        )

        column_ids = {
            str(value).strip()
            for value in expression.columns
        }

        index_ids = {
            str(value).strip()
            for value in expression.index
        }

        column_overlap = column_ids & metadata_ids
        index_overlap = index_ids & metadata_ids

        if len(column_overlap) >= len(index_overlap):
            orientation = "features_by_samples"
            sample_ids = [str(x).strip() for x in expression.columns]
            sample_count = expression.shape[1]
            feature_count = expression.shape[0]
        else:
            orientation = "samples_by_features"
            sample_ids = [str(x).strip() for x in expression.index]
            sample_count = expression.shape[0]
            feature_count = expression.shape[1]

        expression_ids = set(sample_ids)
        platform_sample_sets[platform] = expression_ids

        missing_values = int(expression.isna().sum().sum())
        duplicate_expression_ids = len(sample_ids) - len(expression_ids)

        duplicate_metadata_ids = int(
            metadata[metadata_id_col]
            .dropna()
            .astype(str)
            .str.strip()
            .duplicated()
            .sum()
        )

        metadata_missing_from_expression = metadata_ids - expression_ids
        expression_missing_from_metadata = expression_ids - metadata_ids
        expression_missing_from_harmonized = expression_ids - harmonized_ids

        checks = {
            "sample_count_expected":
                sample_count == files["expected_samples"],
            "feature_count_expected":
                feature_count == files["expected_features"],
            "metadata_expression_alignment":
                len(metadata_missing_from_expression) == 0
                and len(expression_missing_from_metadata) == 0,
            "harmonized_metadata_alignment":
                len(expression_missing_from_harmonized) == 0,
            "no_duplicate_expression_ids":
                duplicate_expression_ids == 0,
            "no_duplicate_metadata_ids":
                duplicate_metadata_ids == 0,
            "no_missing_expression_values":
                missing_values == 0,
        }

        status = "PASS" if all(checks.values()) else "FAIL"

        rows.append({
            "dataset": DATASET,
            "platform": platform,
            "status": status,
            "orientation": orientation,
            "feature_count": feature_count,
            "expected_feature_count": files["expected_features"],
            "sample_count": sample_count,
            "expected_sample_count": files["expected_samples"],
            "metadata_rows": len(metadata),
            "metadata_sample_column": metadata_id_col,
            "metadata_missing_from_expression":
                len(metadata_missing_from_expression),
            "expression_missing_from_metadata":
                len(expression_missing_from_metadata),
            "expression_missing_from_harmonized":
                len(expression_missing_from_harmonized),
            "duplicate_expression_ids":
                duplicate_expression_ids,
            "duplicate_metadata_ids":
                duplicate_metadata_ids,
            "missing_expression_values":
                missing_values,
            **checks,
        })

        print(f"=== {platform} ===")
        print(f"Expression shape: {expression.shape}")
        print(f"Orientation:      {orientation}")
        print(f"Features:         {feature_count}")
        print(f"Samples:          {sample_count}")
        print(f"Metadata rows:    {len(metadata)}")
        print(f"Missing values:   {missing_values}")
        print(f"Status:           {status}")

        for check, passed in checks.items():
            print(f"  {'PASS' if passed else 'FAIL'}  {check}")

        print()

    platform_names = list(platform_sample_sets)

    cross_platform_overlap = (
        platform_sample_sets[platform_names[0]]
        & platform_sample_sets[platform_names[1]]
    )

    total_unique_samples = len(
        set().union(*platform_sample_sets.values())
    )

    summary = pd.DataFrame(rows)

    summary_path = (
        RESULTS / f"{DATASET}_session14_input_validation.tsv"
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    overall_pass = (
        (summary["status"] == "PASS").all()
        and total_unique_samples == 241
        and len(cross_platform_overlap) == 0
    )

    overall = {
        "dataset": DATASET,
        "platforms": platform_names,
        "total_unique_expression_samples": total_unique_samples,
        "expected_total_unique_samples": 241,
        "cross_platform_duplicate_samples":
            len(cross_platform_overlap),
        "harmonized_metadata_rows": len(harmonized),
        "harmonized_unique_sample_ids": len(harmonized_ids),
        "overall_status": "PASS" if overall_pass else "FAIL",
    }

    json_path = (
        LOGS / f"{DATASET}_session14_input_validation.json"
    )
    json_path.write_text(
        json.dumps(overall, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=== OVERALL ===")
    for key, value in overall.items():
        print(f"{key}: {value}")

    print()
    print(f"Validation table: {summary_path}")
    print(f"Validation JSON:  {json_path}")

    if not overall_pass:
        print(
            "\nValidation failed. Do not begin annotation yet.",
            file=sys.stderr,
        )
        return 1

    print(
        "\nValidation passed. Inputs are ready for platform annotation."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
