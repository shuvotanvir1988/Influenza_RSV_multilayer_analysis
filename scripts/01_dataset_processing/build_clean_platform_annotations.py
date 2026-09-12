#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import json
import sys

import pandas as pd


ROOT = Path.home() / "Influenza_RSV_Project"
DATASET = "DS001_GSE38900"

PLATFORM_DIR = (
    ROOT / "data" / "raw" / DATASET / "platform"
)

EXPRESSION_DIR = (
    ROOT / "data" / "interim" / DATASET
)

OUTPUT_DIR = (
    ROOT / "data" / "interim" / DATASET / "annotation"
)

RESULTS_DIR = (
    ROOT
    / "results"
    / DATASET
    / "preprocessing"
    / "tables"
)

LOG_DIR = ROOT / "logs" / "session14"

PLATFORMS = ["GPL10558", "GPL6884"]

RETAIN_COLUMNS = [
    "ID",
    "Species",
    "Source",
    "Search_Key",
    "Transcript",
    "ILMN_Gene",
    "Source_Reference_ID",
    "RefSeq_ID",
    "Unigene_ID",
    "Entrez_Gene_ID",
    "GI",
    "Accession",
    "Symbol",
    "Protein_Product",
    "Probe_Id",
    "Array_Address_Id",
    "Probe_Type",
    "Probe_Start",
    "SEQUENCE",
    "Chromosome",
    "Probe_Chr_Orientation",
    "Probe_Coordinates",
    "Cytoband",
    "Definition",
    "Synonyms",
    "Obsolete_Probe_Id",
    "GB_ACC",
]


def read_geo_platform_table(path: Path) -> pd.DataFrame:
    """Read the annotation table embedded in a GEO SOFT platform file."""

    rows: list[list[str]] = []
    header: list[str] | None = None
    inside = False

    with path.open(
        mode="r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")

            if line == "!platform_table_begin":
                inside = True
                continue

            if line == "!platform_table_end":
                break

            if not inside:
                continue

            fields = line.split("\t")

            if header is None:
                header = fields
                continue

            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))

            if len(fields) > len(header):
                fields = fields[:len(header)]

            rows.append(fields)

    if header is None:
        raise RuntimeError(
            f"No platform annotation table found in {path}"
        )

    table = pd.DataFrame(rows, columns=header)

    if "ID" not in table.columns:
        raise RuntimeError(
            f"Required ID column missing from {path}. "
            f"Columns: {list(table.columns)}"
        )

    return table


def clean_text_column(series: pd.Series) -> pd.Series:
    cleaned = (
        series
        .fillna("")
        .astype(str)
        .str.strip()
    )

    missing_tokens = {
        "",
        "---",
        "NA",
        "N/A",
        "na",
        "n/a",
        "null",
        "None",
    }

    return cleaned.mask(cleaned.isin(missing_tokens), pd.NA)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    qc_rows: list[dict[str, object]] = []
    all_pass = True

    for platform in PLATFORMS:
        platform_path = PLATFORM_DIR / f"{platform}.txt"

        expression_path = (
            EXPRESSION_DIR
            / f"{DATASET}_{platform}_expression_matrix.tsv.gz"
        )

        if not platform_path.exists():
            raise FileNotFoundError(platform_path)

        if not expression_path.exists():
            raise FileNotFoundError(expression_path)

        annotation = read_geo_platform_table(platform_path)

        expression_probe_table = pd.read_csv(
            expression_path,
            sep="\t",
            compression="gzip",
            usecols=[0],
            dtype=str,
        )

        expression_probe_column = expression_probe_table.columns[0]

        expression_probes = (
            expression_probe_table[expression_probe_column]
            .astype(str)
            .str.strip()
        )

        expression_probe_set = set(expression_probes)

        for column in annotation.columns:
            annotation[column] = clean_text_column(
                annotation[column]
            )

        annotation["ID"] = annotation["ID"].astype("string").str.strip()

        duplicate_platform_probe_ids = int(
            annotation["ID"].duplicated().sum()
        )

        if duplicate_platform_probe_ids > 0:
            annotation = annotation.drop_duplicates(
                subset=["ID"],
                keep="first",
            ).copy()

        platform_probe_set = set(annotation["ID"].dropna())

        expression_annotation = annotation[
            annotation["ID"].isin(expression_probe_set)
        ].copy()

        missing_from_annotation = (
            expression_probe_set - platform_probe_set
        )

        extra_platform_probes = (
            platform_probe_set - expression_probe_set
        )

        available_retain_columns = [
            column
            for column in RETAIN_COLUMNS
            if column in expression_annotation.columns
        ]

        clean = expression_annotation[
            available_retain_columns
        ].copy()

        rename_map = {
            "ID": "probe_id",
            "Species": "species",
            "Source": "annotation_source",
            "Search_Key": "search_key",
            "Transcript": "transcript",
            "ILMN_Gene": "illumina_gene",
            "Source_Reference_ID": "source_reference_id",
            "RefSeq_ID": "refseq_id",
            "Unigene_ID": "unigene_id",
            "Entrez_Gene_ID": "entrez_gene_id",
            "GI": "gi",
            "Accession": "accession",
            "Symbol": "gene_symbol",
            "Protein_Product": "protein_product",
            "Probe_Id": "manufacturer_probe_id",
            "Array_Address_Id": "array_address_id",
            "Probe_Type": "probe_type",
            "Probe_Start": "probe_start",
            "SEQUENCE": "probe_sequence",
            "Chromosome": "chromosome",
            "Probe_Chr_Orientation": "probe_chr_orientation",
            "Probe_Coordinates": "probe_coordinates",
            "Cytoband": "cytoband",
            "Definition": "gene_definition",
            "Synonyms": "synonyms",
            "Obsolete_Probe_Id": "obsolete_probe_id",
            "GB_ACC": "genbank_accession",
        }

        clean = clean.rename(columns=rename_map)

        clean.insert(0, "platform", platform)

        expression_order = {
            probe: position
            for position, probe in enumerate(expression_probes)
        }

        clean["_expression_order"] = (
            clean["probe_id"]
            .map(expression_order)
        )

        clean = (
            clean
            .sort_values("_expression_order")
            .drop(columns="_expression_order")
            .reset_index(drop=True)
        )

        clean["has_gene_symbol"] = (
            clean.get(
                "gene_symbol",
                pd.Series(pd.NA, index=clean.index),
            )
            .notna()
        )

        clean["has_entrez_gene_id"] = (
            clean.get(
                "entrez_gene_id",
                pd.Series(pd.NA, index=clean.index),
            )
            .notna()
        )

        clean["has_refseq_id"] = (
            clean.get(
                "refseq_id",
                pd.Series(pd.NA, index=clean.index),
            )
            .notna()
        )

        if "gene_symbol" in clean.columns:
            symbol_counts = (
                clean["gene_symbol"]
                .dropna()
                .value_counts()
            )

            clean["probes_per_gene_symbol"] = (
                clean["gene_symbol"]
                .map(symbol_counts)
                .astype("Int64")
            )

            duplicated_gene_symbol_probes = int(
                clean["gene_symbol"]
                .notna()
                .mul(
                    clean["probes_per_gene_symbol"]
                    .fillna(0)
                    .gt(1)
                )
                .sum()
            )

            unique_gene_symbols = int(
                clean["gene_symbol"]
                .dropna()
                .nunique()
            )

            duplicated_gene_symbols = int(
                (symbol_counts > 1).sum()
            )

        else:
            clean["probes_per_gene_symbol"] = pd.Series(
                pd.NA,
                index=clean.index,
                dtype="Int64",
            )
            duplicated_gene_symbol_probes = 0
            unique_gene_symbols = 0
            duplicated_gene_symbols = 0

        clean_path = (
            OUTPUT_DIR
            / f"{platform}_annotation_clean.tsv"
        )

        clean.to_csv(
            clean_path,
            sep="\t",
            index=False,
            na_rep="",
        )

        platform_qc = {
            "dataset": DATASET,
            "platform": platform,
            "platform_table_rows_original":
                len(annotation),
            "expression_probe_count":
                len(expression_probe_set),
            "matched_expression_probes":
                len(clean),
            "missing_expression_probes":
                len(missing_from_annotation),
            "probe_coverage_percent":
                round(
                    100
                    * len(clean)
                    / len(expression_probe_set),
                    4,
                ),
            "extra_platform_probes":
                len(extra_platform_probes),
            "duplicate_platform_probe_ids":
                duplicate_platform_probe_ids,
            "probes_with_gene_symbol":
                int(clean["has_gene_symbol"].sum()),
            "probes_without_gene_symbol":
                int((~clean["has_gene_symbol"]).sum()),
            "probes_with_entrez_gene_id":
                int(clean["has_entrez_gene_id"].sum()),
            "probes_without_entrez_gene_id":
                int((~clean["has_entrez_gene_id"]).sum()),
            "probes_with_refseq_id":
                int(clean["has_refseq_id"].sum()),
            "unique_gene_symbols":
                unique_gene_symbols,
            "duplicated_gene_symbols":
                duplicated_gene_symbols,
            "probes_assigned_to_duplicated_gene_symbols":
                duplicated_gene_symbol_probes,
            "output_file":
                str(clean_path),
        }

        platform_qc["status"] = (
            "PASS"
            if platform_qc["missing_expression_probes"] == 0
            and platform_qc["matched_expression_probes"]
            == platform_qc["expression_probe_count"]
            else "FAIL"
        )

        if platform_qc["status"] != "PASS":
            all_pass = False

        qc_rows.append(platform_qc)

        individual_qc_path = (
            OUTPUT_DIR / f"{platform}_annotation_qc.tsv"
        )

        pd.DataFrame(
            [platform_qc]
        ).to_csv(
            individual_qc_path,
            sep="\t",
            index=False,
        )

        print(f"\n=== {platform} ===")
        print(
            f"Expression probes: "
            f"{platform_qc['expression_probe_count']:,}"
        )
        print(
            f"Matched probes:    "
            f"{platform_qc['matched_expression_probes']:,}"
        )
        print(
            f"Coverage:          "
            f"{platform_qc['probe_coverage_percent']:.2f}%"
        )
        print(
            f"Gene symbols:      "
            f"{platform_qc['unique_gene_symbols']:,}"
        )
        print(
            f"Probes with symbol:"
            f" {platform_qc['probes_with_gene_symbol']:,}"
        )
        print(
            f"Probes with Entrez:"
            f" {platform_qc['probes_with_entrez_gene_id']:,}"
        )
        print(
            f"Duplicated symbols:"
            f" {platform_qc['duplicated_gene_symbols']:,}"
        )
        print(f"Status:            {platform_qc['status']}")
        print(f"Output:            {clean_path}")

    combined_qc = pd.DataFrame(qc_rows)

    combined_qc_path = (
        RESULTS_DIR
        / f"{DATASET}_platform_annotation_qc_summary.tsv"
    )

    combined_qc.to_csv(
        combined_qc_path,
        sep="\t",
        index=False,
    )

    summary_json = {
        "dataset": DATASET,
        "platform_count": len(PLATFORMS),
        "all_platforms_pass": all_pass,
        "overall_status": "PASS" if all_pass else "FAIL",
        "annotation_qc_summary":
            str(combined_qc_path),
    }

    summary_json_path = (
        LOG_DIR
        / f"{DATASET}_platform_annotation_summary.json"
    )

    summary_json_path.write_text(
        json.dumps(summary_json, indent=2) + "\n",
        encoding="utf-8",
    )

    print("\n=== OVERALL ===")
    print(f"Status: {summary_json['overall_status']}")
    print(f"QC summary: {combined_qc_path}")
    print(f"JSON log:   {summary_json_path}")

    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
