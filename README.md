# Influenza-centered multilayer host-response analysis with comparative RSV analyses

This repository contains the computational code, workflow provenance, software-environment records, and reproducibility materials associated with an influenza-centered multilayer analysis with comparative respiratory syncytial virus (RSV) host-response analyses.

The study integrates discovery transcriptomics, cross-platform replication, pathway analysis, immune-cell composition, transcription-factor activity inference, independent RNA-seq validation, multiseason validation, plasma proteomics, functional CRISPR evidence, single-cell localization, and multilayer evidence integration.

## Repository scope

The repository is organized around the final analyses and manuscript outputs. Where multiple historical versions of a script existed during development, the repository contains the canonical or final frozen version associated with the manuscript wherever recoverable.

The `scripts/` directory is a reproducibility layer containing analytical scripts together with selected frozen specifications, provenance records, and compact result snapshots required to document workflows for which the original execution script could not be recovered.

## Repository structure

The principal directories are:

- `scripts/` — analysis, figure, and table reproducibility materials
- `workflow/` — execution order and script-to-output provenance maps
- `requirements/` — Conda, Python, R, and software-version records
- `provenance/` — SHA256 manifests and frozen decisions
- `docs/` — detailed reproducibility documentation

The `scripts/` reproducibility layer is organized into 13 modules, from dataset processing through manuscript tables.

## Major analysis modules

### 1. Dataset processing

Processing and harmonization of the discovery microarray dataset GSE38900, including platform annotation, phenotype harmonization, probe-level processing, gene-level matrix construction, cross-platform harmonization, and quality control.

### 2. Differential expression

Differential-expression workflows for the discovery cohorts, including primary models, sensitivity analyses, sample inventories, alignment checks, and input validation.

### 3. Cross-platform replication and gene architecture

Cross-platform RSV replication and construction of the frozen 170-gene response architecture, including shared-core, influenza-amplified, and RSV-amplified response classes.

### 4. Pathway analysis

Pathway-level analyses used to characterize influenza-associated and shared host-response programs.

### 5. Cellular composition

Immune-cell composition inference and associated statistical analyses.

### 6. Regulatory analysis

Transcription-factor activity inference using CollecTRI and decoupler, followed by identification and prioritization of regulatory drivers.

### 7. Independent and multiseason validation

Independent RNA-seq validation using GSE155925 and additional multiseason validation workflows.

### 8. Proteomics

Influenza plasma-proteomic validation, including gene-level, pathway-level, and regulatory-program comparisons with discovery transcriptomic results.

### 9. CRISPR functional evidence

Integration of influenza and RSV functional CRISPR-screen evidence.

The original execution script that generated the frozen Session 26 CRISPR result tables was not recovered. This module therefore includes the frozen analysis specifications, provenance documentation, and compact result snapshots needed to document the analysis. Any reconstructed reproduction code must be regarded as reconstructed code rather than the original historical execution script.

### 10. Single-cell analysis

Single-cell localization and supporting analyses using dedicated single-cell workflows and a separate Python environment.

### 11. Multilayer integration

Evidence integration across transcriptomic, proteomic, functional, regulatory, and single-cell layers.

### 12. Figures

Canonical generators for the eight main manuscript figures and eleven supplementary figures.

### 13. Tables

Scripts and provenance associated with manuscript supplementary tables.

## Software environments

Two computational environments are documented separately.

### Main environment

- Conda environment: `influenza_rsv`
- Python: 3.12.13
- R: 4.4.3

The portable file `requirements/environment.yml` was derived from a successful `conda env export --no-builds` export with only the machine-specific `prefix:` line removed.

A `conda env export --from-history` export was not used because that command failed for this environment.

The unmodified full export is retained as `requirements/environment_full.yml`.

Additional package inventories are provided in `requirements/requirements.txt`, `requirements/R_sessionInfo.txt`, `requirements/R_packages.tsv`, and `requirements/software_versions.tsv`.

### Single-cell environment

- Conda environment: `influenza_scrna_census`
- Python: 3.11.15

The corresponding environment files are `requirements/environment_scrna.yml`, `requirements/environment_scrna_full.yml`, and `requirements/requirements_scrna.txt`.

This environment was used for the dedicated single-cell census and localization workflows.

## Reproducibility and provenance

The intended analytical execution order is documented in `workflow/ANALYSIS_EXECUTION_ORDER.md`.

The provenance maps `workflow/SCRIPT_TO_OUTPUT_MAP.tsv`, `workflow/FIGURE_TO_SCRIPT_MAP.tsv`, and `workflow/TABLE_TO_SCRIPT_MAP.tsv` connect repository materials to manuscript outputs.

The repository script reproducibility layer is protected by SHA256 manifests under `provenance/manifests/`.

The software-environment files have a separate SHA256 manifest at `requirements/REQUIREMENTS_SHA256SUMS_v1.0.txt`.

These manifests allow later changes or accidental substitutions of frozen files to be detected.

## Integrity checks

Environment files can be verified with `cd requirements && sha256sum -c REQUIREMENTS_SHA256SUMS_v1.0.txt`.

The repository script layer can be verified using the corresponding SHA256 manifest under `provenance/manifests/`.

## Version-selection principle

The project was developed iteratively, and some analyses underwent multiple versions before freezing.

Final canonical versions were selected using, in order of priority:

1. explicit association with final frozen outputs;
2. output-hash agreement with manuscript-facing files;
3. documented frozen-version provenance;
4. later session-level canonical records;
5. filename versioning only as supporting evidence.

Historical development copies and clearly superseded scripts are intentionally excluded from the public reproducibility layer.

## Data sources

This study primarily reanalyzes publicly available datasets. Dataset identifiers, source repositories, and access information are described in the manuscript and in the repository Data Availability documentation.

Major data sources include:

- GEO GSE38900 — discovery transcriptomic dataset
- GEO GSE155925 — independent RNA-seq validation dataset
- GEO GSE149689 — single-cell localization dataset
- GEO GSE283746 — additional single-cell/regulatory dataset
- Figshare record 27826857 — influenza SomaScan proteomics dataset
- published influenza and RSV CRISPR-screen datasets used for functional evidence integration

Raw third-party datasets are not redistributed in this repository unless redistribution is explicitly permitted by the original source. Users should obtain source datasets from their original public repositories.

## Reproducing the workflow

Because the study integrates multiple datasets and analytical layers, complete reproduction is intended to proceed module by module rather than through a single monolithic command.

Start with `workflow/ANALYSIS_EXECUTION_ORDER.md` and then consult `docs/reproducibility_guide.md` for environment setup, execution order, provenance checks, and known limitations.

The workflow maps in `workflow/` should be used to trace analytical scripts to manuscript outputs.

## Known reproducibility limitation

The original execution script that generated the frozen Session 26 CRISPR result tables was not recovered. The public repository therefore preserves the frozen analysis specification, provenance documentation, and compact result snapshots for that module. Any future reconstructed script must be clearly identified as reconstructed reproduction code and not as the original historical execution script.

## Code and data availability

Detailed code-availability and data-availability information is provided in `CODE_AVAILABILITY.md` and `DATA_AVAILABILITY.md`.

The public GitHub repository is available at https://github.com/shuvotanvir1988/Influenza_RSV_multilayer_analysis. A permanent archival DOI will be added after the archived release is created.

## License

Original repository code and documentation are released under the MIT License; see `LICENSE`. Third-party materials remain subject to their original terms and licenses.

## Citation

Citation metadata will be provided in `CITATION.cff`. Once a permanent archival DOI is assigned, that DOI should be used when citing the code release.

## Contact

Questions about the computational workflow can be submitted through the public repository issue tracker after release or directed to the repository authors using the contact information associated with the published article.
