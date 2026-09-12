# Reproducibility guide

This guide describes how to reproduce the computational analyses associated with the influenza/RSV multilayer host-response study.

## 1. General principle

The project was developed iteratively and contains multiple analytical layers. Reproduction should therefore proceed module by module rather than through a single monolithic command.

The authoritative execution order is documented in `workflow/ANALYSIS_EXECUTION_ORDER.md`.

The repository preserves final canonical or frozen scripts wherever recoverable. Historical development copies and superseded versions are intentionally excluded from the public reproducibility layer.

## 2. Repository integrity

Before executing analyses, verify the frozen repository materials.

### Script layer

Use the SHA256 manifest under `provenance/manifests/` to verify the repository `scripts/` reproducibility layer.

### Environment layer

Verify the environment files with:

`cd requirements && sha256sum -c REQUIREMENTS_SHA256SUMS_v1.0.txt`

Any unexpected hash mismatch should be investigated before proceeding.

## 3. Main computational environment

The main analysis environment is `influenza_rsv`.

- Python: 3.12.13
- R: 4.4.3

The portable Conda specification is `requirements/environment.yml`.

The file was derived from a successful `conda env export --no-builds` export with only the machine-specific `prefix:` line removed.

The original full export is retained as `requirements/environment_full.yml`.

A `conda env export --from-history` export was not used because that command failed for the environment.

Additional package inventories are available in:

- `requirements/requirements.txt`
- `requirements/R_sessionInfo.txt`
- `requirements/R_packages.tsv`
- `requirements/software_versions.tsv`

## 4. Single-cell computational environment

Single-cell census and localization workflows use the separate environment `influenza_scrna_census`.

- Python: 3.11.15

Environment records are:

- `requirements/environment_scrna.yml`
- `requirements/environment_scrna_full.yml`
- `requirements/requirements_scrna.txt`

This environment should be treated separately from the main analysis environment.

## 5. Environment reconstruction

A user may attempt to reconstruct the main environment with:

`conda env create -f requirements/environment.yml`

and the single-cell environment with:

`conda env create -f requirements/environment_scrna.yml`

Because software repositories change over time, exact historical reconstruction may depend on package availability. The full exports, Python package inventories, R package inventory, and curated software-version table are therefore retained as additional provenance.

## 6. Software-version caveats

The curated software table is provided in `requirements/software_versions.tsv`.

The current main Python runtime reports NumPy 2.5.2, whereas the captured package-manager record contains NumPy 2.5.1. The runtime value is recorded in the curated software table, while the original environment export is preserved unchanged.

The regulatory-analysis environment provenance records decoupler 2.2.0. decoupler is not currently installed in the captured main environment, so its version is reported from frozen regulatory provenance rather than from the present runtime.

statsmodels 0.14.6 is present in the captured single-cell environment and is also documented in historical regulatory provenance. Earlier Session 38 execution logs include failed attempts in which statsmodels was not available; these failed historical logs are retained rather than rewritten.

## 7. Analysis execution order

The intended high-level order is:

1. dataset processing;
2. differential expression;
3. cross-platform replication and frozen 170-gene architecture;
4. pathway analysis;
5. cellular-composition analysis;
6. regulatory analysis;
7. independent RNA-seq and multiseason validation;
8. proteomic validation;
9. CRISPR functional-evidence integration;
10. single-cell analysis;
11. multilayer integration;
12. figure generation;
13. table generation.

Refer to `workflow/ANALYSIS_EXECUTION_ORDER.md` for the detailed module-level ordering.

## 8. Script-to-output provenance

Three workflow maps are supplied:

- `workflow/SCRIPT_TO_OUTPUT_MAP.tsv`
- `workflow/FIGURE_TO_SCRIPT_MAP.tsv`
- `workflow/TABLE_TO_SCRIPT_MAP.tsv`

These maps should be used to trace repository scripts and provenance files to manuscript-facing outputs.

## 9. Discovery dataset processing

The discovery dataset is GEO GSE38900.

The corresponding repository module is `scripts/01_dataset_processing/`.

The processing workflow includes platform annotation, phenotype harmonization, Illumina normalization, probe parsing, gene-level matrix construction, cross-platform harmonization, and quality-control checks.

The final gene-level matrices use the frozen highest-IQR probe-collapse strategy documented by the canonical processing scripts.

## 10. Differential expression

Primary and supporting differential-expression workflows are stored under `scripts/02_differential_expression/`.

The repository contains the final analytical scripts together with selected validation and sensitivity scripts needed to document the manuscript results.

## 11. Cross-platform replication and 170-gene architecture

The relevant module is `scripts/03_cross_platform_replication/`.

The `gene_architecture/` subdirectory contains the canonical scripts used to construct and audit the frozen 170-gene response architecture.

The authoritative gene-level evidence table is the frozen DS001/GSE38900 gene-evidence master produced by this workflow.

## 12. Pathway, cellular, and regulatory analyses

Pathway scripts are in `scripts/04_pathway_analysis/`.

Cellular-composition scripts are in `scripts/05_cellular_composition/`.

Regulatory-analysis scripts are in `scripts/06_regulatory_analysis/`.

The regulatory workflow uses frozen CollecTRI network provenance and decoupler-based transcription-factor activity inference.

## 13. Independent RNA-seq and multiseason validation

Independent GSE155925 validation and multiseason validation workflows are stored under `scripts/07_rnaseq_validation/`.

Primary and sensitivity analyses should be distinguished according to the frozen workflow documentation and manuscript methods.

## 14. Proteomics

Proteomic validation scripts are stored under `scripts/08_proteomics/`.

The analysis integrates the influenza SomaScan dataset with transcriptomic gene-level, pathway-level, and regulatory-program results.

Raw third-party proteomic data should be obtained from the original public source rather than redistributed from this repository unless redistribution is explicitly permitted.

## 15. CRISPR provenance limitation

CRISPR materials are stored under `scripts/09_crispr/`.

The original execution script that generated the frozen Session 26 CRISPR result tables was not recovered.

The repository therefore includes frozen analysis specifications, provenance records, and compact result snapshots that document the functional-evidence layer.

Any future script reconstructed to reproduce these results must be labeled `RECONSTRUCTED_REPRODUCTION_SCRIPT` and must not be represented as the original historical execution script.

## 16. Single-cell analysis

Single-cell workflows are stored under `scripts/10_single_cell/`.

These analyses use the separate `influenza_scrna_census` environment described above.

## 17. Multilayer integration

Integrated evidence-analysis scripts are stored under `scripts/11_multilayer_integration/`.

These scripts combine frozen upstream evidence from transcriptomic, proteomic, regulatory, CRISPR, and single-cell layers.

Upstream scientific results should not be silently regenerated or altered when reproducing downstream integration outputs.

## 18. Figure generation

Canonical figure generators are stored under `scripts/12_figures/`.

The repository contains generators for:

- 8 main manuscript figures;
- 11 supplementary figures.

Use `workflow/FIGURE_TO_SCRIPT_MAP.tsv` to identify the canonical generator corresponding to each manuscript figure.

## 19. Table generation

Supplementary-table scripts and provenance are stored under `scripts/13_tables/`.

Use `workflow/TABLE_TO_SCRIPT_MAP.tsv` to trace table outputs to their corresponding workflows.

## 20. Frozen results and sensitivity analyses

Primary, sensitivity, exploratory, and provenance-only scripts are not interchangeable.

When reproducing manuscript results, use the script classification and output mapping supplied by the workflow and provenance files rather than selecting scripts solely by filename or apparent version number.

## 21. Third-party data

Most source data are publicly available through GEO, Figshare, or the original publications.

Raw third-party datasets should be downloaded from their original repositories. This repository is intended to provide the computational workflow and reproducibility materials rather than redistribute all source data.

## 22. Recommended reproducibility procedure

A reviewer or external researcher should:

1. verify SHA256 manifests;
2. reconstruct the required environment;
3. obtain the original public datasets;
4. follow `workflow/ANALYSIS_EXECUTION_ORDER.md`;
5. use the workflow maps to identify canonical scripts;
6. compare regenerated outputs against frozen manuscript-facing outputs or reported summary statistics;
7. document any package or platform differences encountered during reproduction.

## 23. Known limitations

Exact reproduction can be affected by changes in third-party repositories, package availability, operating-system libraries, and remote resources.

The CRISPR layer has the additional limitation that the original historical execution script was not recovered.

The repository preserves available provenance and frozen outputs so that these limitations are explicit rather than hidden.
