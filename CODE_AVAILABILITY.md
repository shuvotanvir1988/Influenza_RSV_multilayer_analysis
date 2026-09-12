# Code availability

The computational code and reproducibility materials supporting this study are organized in this repository and are intended for public archival release with the associated manuscript.

## Repository contents

The repository includes:

- canonical analysis scripts associated with the final manuscript results;
- scripts used to generate the eight main and eleven supplementary figures;
- supplementary-table generation and provenance materials;
- workflow maps linking scripts to manuscript outputs;
- Conda, Python, and R software-environment records;
- SHA256 manifests for integrity verification;
- documentation describing the intended analysis execution order and reproducibility procedure; and
- selected frozen specifications, provenance records, and compact result snapshots required to document analytical steps for which the original execution script could not be recovered.

The principal analysis modules cover discovery-dataset processing, differential expression, cross-platform replication, construction of the frozen 170-gene response architecture, pathway analysis, cellular-composition analysis, regulatory inference, independent RNA-seq and multiseason validation, proteomic validation, CRISPR functional-evidence integration, single-cell analysis, multilayer integration, figure generation, and supplementary-table generation.

## Reproducibility documentation

The recommended analysis order is documented in `workflow/ANALYSIS_EXECUTION_ORDER.md`.

Script-to-output relationships are documented in:

- `workflow/SCRIPT_TO_OUTPUT_MAP.tsv`
- `workflow/FIGURE_TO_SCRIPT_MAP.tsv`
- `workflow/TABLE_TO_SCRIPT_MAP.tsv`

Detailed reproduction guidance is provided in `docs/reproducibility_guide.md`.

## Software environments

The main computational environment and the separate single-cell environment are documented under `requirements/` using Conda environment exports, Python package inventories, R package inventories, and a curated software-version table.

SHA256 manifests are supplied to permit verification of the archived reproducibility materials.

## CRISPR provenance limitation

The original execution script that generated the frozen Session 26 CRISPR result tables was not recovered.

For this analytical layer, the repository therefore preserves the frozen analysis specifications, provenance documentation, and compact result snapshots available from the completed analysis. Any future reconstructed reproduction script must be explicitly identified as reconstructed code and must not be represented as the original historical execution script.

This limitation does not apply to the other analytical modules for which canonical or final frozen scripts were recovered and included in the repository.

## Public archive

The final public repository URL and permanent archival DOI will be added here after repository release and DOI registration.

Repository: [TO BE ADDED AFTER PUBLIC RELEASE]

Archived release DOI: [TO BE ADDED AFTER DOI REGISTRATION]

## Manuscript Code Availability statement

Upon public release, the manuscript Code Availability statement can be finalized as follows:

> Custom code and reproducibility materials supporting this study are publicly available at [REPOSITORY URL] and have been archived at [ARCHIVE/DOI]. The repository includes canonical analysis and figure-generation scripts, workflow provenance maps, software-environment records, and SHA256 integrity manifests. The original execution script for the historical CRISPR integration step was not recovered; frozen analysis specifications, provenance records, and result snapshots for this step are provided in the repository. No other custom code restrictions apply.

The bracketed repository and archival identifiers should be replaced only after the corresponding public records exist.
