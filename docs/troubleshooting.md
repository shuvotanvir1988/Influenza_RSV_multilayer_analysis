# Troubleshooting

This document summarizes known reproducibility and environment issues for the Influenza/RSV multilayer analysis repository.

## 1. Run commands from the repository root

Many scripts use repository-relative paths or the `INFLUENZA_RSV_PROJECT_ROOT` environment variable.

Recommended:

```bash
cd Influenza_RSV_multilayer_analysis
export INFLUENZA_RSV_PROJECT_ROOT="$PWD"
```

## 2. Main and single-cell environments are separate

The main analysis environment and the single-cell environment were captured separately.

- Main environment: `influenza_rsv`
- Single-cell environment: `influenza_scrna_census`

Use the environment files in `requirements/` and consult `docs/reproducibility_guide.md` for details.

## 3. Conda environment recreation

`environment.yml` and `environment_scrna.yml` are the portable environment records. The corresponding `*_full.yml` files are retained as historical environment captures and may contain machine-specific installation prefixes.

Machine-specific prefixes in full environment records are provenance metadata and are not intended to be reused directly.

## 4. pip requirements contain conda-build file paths

Some entries in `requirements.txt` and `requirements_scrna.txt` contain `file:///home/conda/feedstock_root/...` references inherited from the original Conda installations.

These are environment provenance records. Prefer the Conda YAML files when recreating the environments.

## 5. NumPy version discrepancy in the main environment

Runtime inspection recorded NumPy 2.5.2, whereas Conda metadata reported 2.5.1. This discrepancy is documented rather than altered after analysis completion.

## 6. decoupler and statsmodels regulatory-analysis provenance

The frozen regulatory-analysis provenance records decoupler 2.2.0 and statsmodels 0.14.6. These packages were not importable in the final captured main environment at repository assembly time.

The original frozen regulatory outputs and scripts were therefore preserved without modifying the scientific environment retrospectively.

## 7. CRISPR historical execution script

The original execution script that generated the historical Session 26 CRISPR integration result tables could not be recovered.

The repository instead provides frozen analysis specifications, dataset and input manifests, provenance records, and compact snapshots of the final CRISPR results.

No reconstructed script is represented as the original historical execution script.

## 8. Public CRISPR provenance trace

The original producer-trace artifact contained machine-specific paths and shell-history material. It was retained privately and replaced in the public repository with a privacy-safe provenance summary.

The SHA256 of the original historical artifact is recorded in `provenance/CRISPR_TRACE_PUBLIC_SANITIZATION_v1.tsv`.

## 9. Figure generation

Canonical figure-generation scripts are listed in `workflow/FIGURE_TO_SCRIPT_MAP.tsv`.

Final manuscript figures were generated from frozen analysis outputs. Scientific results should not be changed merely to reproduce visual formatting.

## 10. Supplementary table generation

Supplementary-table scientific producers and packaging scripts are listed in `workflow/TABLE_TO_SCRIPT_MAP.tsv`.

The Session 39 staging shell scripts were portability-patched in the public repository; original pre-patch hashes are preserved in the provenance records.

## 11. Integrity verification

To verify the reproducibility-layer files:

```bash
sha256sum -c provenance/manifests/SCRIPT_SHA256SUMS_v1.0.txt
```

A repository-wide SHA256 manifest is generated only after the public release package is finalized.
