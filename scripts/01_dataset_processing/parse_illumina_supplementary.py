#!/usr/bin/env python3
"""
Script: parse_illumina_supplementary.py
Project: Influenza A and RSV Systems Immunology
Session: 14, revised Step 2
Version: 1.0

Purpose
-------
Parse GEO Illumina supplementary exports containing alternating signal and
detection P-value columns. Produce separate signal and detection matrices for
GPL10558 and GPL6884.

Inputs
------
data/raw/DS001_GSE38900/supplementary/
    GSE38900_non-normalized_GSM1226237-GSM1226272.txt.gz
    GSE38900_non-normalized.txt.gz

Outputs
-------
data/interim/DS001_GSE38900/illumina_raw/
    DS001_GSE38900_GPL10558_signal_matrix.tsv.gz
    DS001_GSE38900_GPL10558_detection_pvalues.tsv.gz
    DS001_GSE38900_GPL6884_signal_matrix.tsv.gz
    DS001_GSE38900_GPL6884_detection_pvalues.tsv.gz

results/DS001_GSE38900/preprocessing/tables/
    DS001_GSE38900_illumina_parsing_qc.tsv

logs/session14/
    DS001_GSE38900_illumina_parsing_summary.json
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import json
import sys
import time

import numpy as np
import pandas as pd


ROOT = Path.home() / "Influenza_RSV_Project"
DATASET = "DS001_GSE38900"

SUPPLEMENTARY_DIR = (
    ROOT / "data" / "raw" / DATASET / "supplementary"
)

OUTPUT_DIR = (
    ROOT / "data" / "interim" / DATASET / "illumina_raw"
)

RESULTS_DIR = (
    ROOT
    / "results"
    / DATASET
    / "preprocessing"
    / "tables"
)

LOG_DIR = ROOT / "logs" / "session14"


PLATFORMS = {
    "GPL10558": {
        "input_file":
            SUPPLEMENTARY_DIR
            / "GSE38900_non-normalized_GSM1226237-GSM1226272.txt.gz",
        "expected_probes": 47323,
        "expected_samples": 36,
    },
    "GPL6884": {
        "input_file":
            SUPPLEMENTARY_DIR
            / "GSE38900_non-normalized.txt.gz",
        "expected_probes": 48803,
        "expected_samples": 205,
    },
}


def normalize_header(value: object) -> str:
    """Return a stripped string representation of a header."""
    return str(value).strip()


def is_detection_column(header: str) -> bool:
    """Identify detection P-value columns."""
    lowered = header.strip().lower()

    return (
        "detection" in lowered
        and (
            "pval" in lowered
            or "p-value" in lowered
            or "p value" in lowered
        )
    )


def make_unique_sample_names(
    raw_names: list[str],
    platform: str,
) -> list[str]:
    """
    Standardize and verify sample names.

    This function intentionally preserves the original sample labels except for
    surrounding whitespace. Duplicate labels are not silently renamed.
    """
    cleaned = [name.strip() for name in raw_names]

    duplicates = pd.Series(cleaned).duplicated(keep=False)

    if duplicates.any():
        duplicate_names = sorted(
            pd.Series(cleaned)[duplicates].unique().tolist()
        )

        raise RuntimeError(
            f"{platform}: duplicate signal sample labels detected: "
            f"{duplicate_names[:20]}"
        )

    return cleaned


def parse_platform(
    platform: str,
    settings: dict[str, object],
) -> dict[str, object]:
    input_path = Path(settings["input_file"])
    expected_probes = int(settings["expected_probes"])
    expected_samples = int(settings["expected_samples"])

    if not input_path.exists():
        raise FileNotFoundError(
            f"{platform}: input file not found: {input_path}"
        )

    print(f"\n=== {platform} ===")
    print(f"Input: {input_path}")

    table = pd.read_csv(
        input_path,
        sep="\t",
        compression="gzip",
        dtype=str,
        low_memory=False,
    )

    if table.shape[1] < 3:
        raise RuntimeError(
            f"{platform}: fewer than three columns were found."
        )

    probe_column = table.columns[0]

    if normalize_header(probe_column).upper() != "ID_REF":
        raise RuntimeError(
            f"{platform}: expected first column ID_REF, "
            f"found {probe_column!r}"
        )

    remaining_columns = list(table.columns[1:])

    expected_total_columns = 1 + (2 * expected_samples)

    if table.shape[1] != expected_total_columns:
        raise RuntimeError(
            f"{platform}: expected {expected_total_columns} columns, "
            f"found {table.shape[1]}"
        )

    signal_column_positions: list[int] = []
    detection_column_positions: list[int] = []
    signal_headers: list[str] = []

    for position, header in enumerate(
        remaining_columns,
        start=1,
    ):
        normalized_header = normalize_header(header)

        if position % 2 == 1:
            if is_detection_column(normalized_header):
                raise RuntimeError(
                    f"{platform}: column {position + 1} was expected "
                    f"to be a signal column but was detection-like: "
                    f"{normalized_header!r}"
                )

            signal_column_positions.append(position)
            signal_headers.append(normalized_header)

        else:
            if not is_detection_column(normalized_header):
                raise RuntimeError(
                    f"{platform}: column {position + 1} was expected "
                    f"to be a detection P-value column but was: "
                    f"{normalized_header!r}"
                )

            detection_column_positions.append(position)

    if len(signal_column_positions) != expected_samples:
        raise RuntimeError(
            f"{platform}: expected {expected_samples} signal columns, "
            f"found {len(signal_column_positions)}"
        )

    if len(detection_column_positions) != expected_samples:
        raise RuntimeError(
            f"{platform}: expected {expected_samples} detection columns, "
            f"found {len(detection_column_positions)}"
        )

    sample_names = make_unique_sample_names(
        signal_headers,
        platform,
    )

    probe_ids = (
        table.iloc[:, 0]
        .astype(str)
        .str.strip()
    )

    missing_probe_ids = int(
        probe_ids.isna().sum()
        + (probe_ids == "").sum()
    )

    duplicate_probe_ids = int(
        probe_ids.duplicated().sum()
    )

    if missing_probe_ids > 0:
        raise RuntimeError(
            f"{platform}: {missing_probe_ids} missing probe IDs found."
        )

    if duplicate_probe_ids > 0:
        raise RuntimeError(
            f"{platform}: {duplicate_probe_ids} duplicate probe IDs found."
        )

    signal_raw = table.iloc[
        :,
        signal_column_positions,
    ].copy()

    detection_raw = table.iloc[
        :,
        detection_column_positions,
    ].copy()

    signal_raw.columns = sample_names
    detection_raw.columns = sample_names

    signal = signal_raw.apply(
        pd.to_numeric,
        errors="coerce",
    )

    detection = detection_raw.apply(
        pd.to_numeric,
        errors="coerce",
    )

    signal.index = probe_ids
    detection.index = probe_ids

    signal.index.name = "ID_REF"
    detection.index.name = "ID_REF"

    signal_missing = int(
        signal.isna().sum().sum()
    )

    detection_missing = int(
        detection.isna().sum().sum()
    )

    nonnumeric_signal_values = signal_missing
    nonnumeric_detection_values = detection_missing

    detection_below_zero = int(
        (detection < 0).sum().sum()
    )

    detection_above_one = int(
        (detection > 1).sum().sum()
    )

    if detection_below_zero > 0 or detection_above_one > 0:
        raise RuntimeError(
            f"{platform}: detection P-values outside [0,1]. "
            f"Below 0: {detection_below_zero}; "
            f"above 1: {detection_above_one}"
        )

    if signal.shape != detection.shape:
        raise RuntimeError(
            f"{platform}: signal and detection matrix shapes differ. "
            f"Signal: {signal.shape}; detection: {detection.shape}"
        )

    if list(signal.columns) != list(detection.columns):
        raise RuntimeError(
            f"{platform}: signal and detection sample orders differ."
        )

    if list(signal.index) != list(detection.index):
        raise RuntimeError(
            f"{platform}: signal and detection probe orders differ."
        )

    probe_count_pass = signal.shape[0] == expected_probes
    sample_count_pass = signal.shape[1] == expected_samples

    if not probe_count_pass:
        raise RuntimeError(
            f"{platform}: expected {expected_probes} probes, "
            f"found {signal.shape[0]}"
        )

    if not sample_count_pass:
        raise RuntimeError(
            f"{platform}: expected {expected_samples} samples, "
            f"found {signal.shape[1]}"
        )

    signal_path = (
        OUTPUT_DIR
        / f"{DATASET}_{platform}_signal_matrix.tsv.gz"
    )

    detection_path = (
        OUTPUT_DIR
        / f"{DATASET}_{platform}_detection_pvalues.tsv.gz"
    )

    signal.to_csv(
        signal_path,
        sep="\t",
        compression="gzip",
        index=True,
        na_rep="",
    )

    detection.to_csv(
        detection_path,
        sep="\t",
        compression="gzip",
        index=True,
        na_rep="",
    )

    finite_signal = signal.to_numpy(dtype=float)
    finite_detection = detection.to_numpy(dtype=float)

    signal_values = finite_signal[
        np.isfinite(finite_signal)
    ]

    detection_values = finite_detection[
        np.isfinite(finite_detection)
    ]

    qc = {
        "dataset": DATASET,
        "platform": platform,
        "input_file": str(input_path),
        "input_rows": int(table.shape[0]),
        "input_columns": int(table.shape[1]),
        "probe_count": int(signal.shape[0]),
        "expected_probe_count": expected_probes,
        "sample_count": int(signal.shape[1]),
        "expected_sample_count": expected_samples,
        "signal_value_count": int(signal.size),
        "detection_value_count": int(detection.size),
        "duplicate_probe_ids": duplicate_probe_ids,
        "duplicate_sample_names": 0,
        "missing_probe_ids": missing_probe_ids,
        "missing_signal_values": signal_missing,
        "missing_detection_values": detection_missing,
        "nonnumeric_signal_values": nonnumeric_signal_values,
        "nonnumeric_detection_values":
            nonnumeric_detection_values,
        "signal_min": float(np.min(signal_values)),
        "signal_q01": float(np.quantile(signal_values, 0.01)),
        "signal_median": float(np.median(signal_values)),
        "signal_q99": float(np.quantile(signal_values, 0.99)),
        "signal_max": float(np.max(signal_values)),
        "negative_signal_values": int(
            np.sum(signal_values < 0)
        ),
        "detection_min": float(np.min(detection_values)),
        "detection_median": float(np.median(detection_values)),
        "detection_max": float(np.max(detection_values)),
        "detection_values_below_zero":
            detection_below_zero,
        "detection_values_above_one":
            detection_above_one,
        "signal_output_file": str(signal_path),
        "detection_output_file": str(detection_path),
        "status": "PASS",
    }

    print(f"Probes:                 {signal.shape[0]:,}")
    print(f"Samples:                {signal.shape[1]:,}")
    print(f"Missing signal values:  {signal_missing:,}")
    print(f"Missing detection:      {detection_missing:,}")
    print(
        f"Signal range:           "
        f"{qc['signal_min']:.6g} to {qc['signal_max']:.6g}"
    )
    print(
        f"Detection range:        "
        f"{qc['detection_min']:.6g} to "
        f"{qc['detection_max']:.6g}"
    )
    print("Status:                 PASS")
    print(f"Signal output:          {signal_path}")
    print(f"Detection output:       {detection_path}")

    return qc


def main() -> int:
    started = time.time()
    started_at = datetime.now().astimezone()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    print(
        f"Session 14 revised Step 2 started: "
        f"{started_at.isoformat()}"
    )

    qc_rows: list[dict[str, object]] = []

    try:
        for platform, settings in PLATFORMS.items():
            qc_rows.append(
                parse_platform(platform, settings)
            )

        qc_table = pd.DataFrame(qc_rows)

        qc_path = (
            RESULTS_DIR
            / f"{DATASET}_illumina_parsing_qc.tsv"
        )

        qc_table.to_csv(
            qc_path,
            sep="\t",
            index=False,
        )

        finished_at = datetime.now().astimezone()
        elapsed_seconds = time.time() - started

        summary = {
            "dataset": DATASET,
            "session": "14",
            "step": "revised Step 2",
            "started_at": started_at.isoformat(),
            "finished_at": finished_at.isoformat(),
            "elapsed_seconds": round(elapsed_seconds, 3),
            "platform_count": len(qc_rows),
            "platforms": [
                row["platform"]
                for row in qc_rows
            ],
            "all_platforms_pass": bool(
                all(row["status"] == "PASS" for row in qc_rows)
            ),
            "overall_status": "PASS",
            "qc_table": str(qc_path),
        }

        json_path = (
            LOG_DIR
            / f"{DATASET}_illumina_parsing_summary.json"
        )

        json_path.write_text(
            json.dumps(summary, indent=2) + "\n",
            encoding="utf-8",
        )

        print("\n=== OVERALL ===")
        print("Status: PASS")
        print(f"QC table: {qc_path}")
        print(f"JSON log: {json_path}")
        print(f"Elapsed seconds: {elapsed_seconds:.2f}")

        return 0

    except Exception as error:
        finished_at = datetime.now().astimezone()
        elapsed_seconds = time.time() - started

        failure_summary = {
            "dataset": DATASET,
            "session": "14",
            "step": "revised Step 2",
            "started_at": started_at.isoformat(),
            "finished_at": finished_at.isoformat(),
            "elapsed_seconds": round(elapsed_seconds, 3),
            "overall_status": "FAIL",
            "error_type": type(error).__name__,
            "error_message": str(error),
        }

        failure_path = (
            LOG_DIR
            / f"{DATASET}_illumina_parsing_failure.json"
        )

        failure_path.write_text(
            json.dumps(failure_summary, indent=2) + "\n",
            encoding="utf-8",
        )

        print(
            f"\nERROR: {type(error).__name__}: {error}",
            file=sys.stderr,
        )
        print(
            f"Failure log: {failure_path}",
            file=sys.stderr,
        )

        return 1


if __name__ == "__main__":
    raise SystemExit(main())
