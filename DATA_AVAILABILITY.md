# Data availability

This study reanalyzes publicly available transcriptomic, proteomic, single-cell, and functional-screen datasets. Raw third-party datasets are generally not redistributed in this repository; users should obtain them from the original public repositories or associated publications.

## Discovery transcriptomic dataset

The principal discovery dataset is GEO accession `GSE38900`.

The study uses data from the two relevant Illumina platforms represented in this accession, followed by phenotype harmonization, probe processing, gene-level matrix construction, cross-platform replication, and downstream multilayer analyses.

The repository contains the computational workflows and provenance required to reproduce the processing steps but does not replace the original GEO record as the authoritative source of the raw data.

## Independent RNA-seq validation

Independent RSV RNA-seq validation uses GEO accession `GSE155925`.

The original sequencing data and associated metadata should be obtained from GEO and the corresponding source publication. Repository scripts document cohort eligibility, covariate handling, differential-expression modeling, and sensitivity analyses used in this study.

## Single-cell datasets

Single-cell analyses use publicly available datasets including:

- GEO `GSE149689`
- GEO `GSE283746`

The original single-cell data should be retrieved from the corresponding public repositories. The repository contains the analysis scripts, environment records, and provenance used for localization and supporting single-cell analyses.

## Influenza proteomics

Influenza proteomic validation uses the publicly available SomaScan dataset associated with Figshare record `27826857`.

The repository contains code for the study-specific proteomic analyses and cross-layer comparisons. Raw third-party proteomic data should be obtained directly from the original Figshare record unless redistribution is explicitly permitted by its original terms.

## CRISPR functional screens

Functional evidence was derived from published influenza and RSV CRISPR-screen datasets described in the manuscript and associated source publications.

The repository includes the frozen analysis specifications, provenance materials, and compact result snapshots used for the manuscript integration. The original source data remain governed by the availability terms of the corresponding publications and repositories.

## Derived study outputs

Study-derived intermediate and summary outputs required to document the manuscript analyses may be included in the public repository where redistribution is appropriate and file size permits.

The repository workflow maps and provenance files identify the relationship between analytical scripts and manuscript-facing outputs.

## Third-party data redistribution

Public availability of a source dataset does not automatically imply that every downloaded copy should be redistributed in a secondary code repository. Accordingly, this repository prioritizes reproducible retrieval from the authoritative public source together with preservation of analysis code, provenance, software environments, and compact derived outputs.

Users are responsible for complying with the access, citation, licensing, and reuse conditions specified by the original data providers.

## Permanent links and identifiers

The public code-repository URL is provided below. The permanent archival DOI will be added after archival release.

Repository: https://github.com/shuvotanvir1988/Influenza_RSV_multilayer_analysis

Archived release DOI: [TO BE ADDED AFTER DOI REGISTRATION]

## Manuscript Data Availability statement

A concise manuscript-facing Data Availability statement can be finalized after permanent archival as follows:

> The datasets analyzed in this study are publicly available from their original repositories. Discovery transcriptomic data are available from GEO under accession GSE38900; independent RNA-seq validation data are available under GSE155925; single-cell datasets include GSE149689 and GSE283746; and influenza SomaScan proteomic data are available through Figshare record 27826857. Published influenza and RSV CRISPR-screen datasets were obtained from the sources cited in the manuscript. Study-specific code, provenance records, and reproducibility materials are available at https://github.com/shuvotanvir1988/Influenza_RSV_multilayer_analysis and archived at [ARCHIVE/DOI].

The archival DOI placeholders should be replaced only after the permanent archived release exists.
