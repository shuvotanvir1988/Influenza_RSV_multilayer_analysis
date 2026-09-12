# Analysis Execution Order

This document describes the intended scientific execution order for the reproducibility repository accompanying the influenza/RSV multilayer analysis.

The repository contains canonical analysis scripts, quality-control and sensitivity scripts, figure generators, table-generation utilities, and selected frozen provenance artifacts. Not every file is intended to be executed sequentially.

---

## 1. Discovery dataset processing

Directory:

`scripts/01_dataset_processing/`

Purpose: prepare the GSE38900 discovery microarray datasets, harmonize phenotypes and platform annotations, normalize expression data, collapse probes to genes, and perform quality control.

Recommended sequence:

1. `validate_ds001_inputs.py`
2. `extract_gse38900_metadata.py`
3. `harmonize_gse38900_phenotypes.py`
4. `download_platform_annotations.py`
5. `build_clean_platform_annotations.py`
6. `parse_illumina_supplementary.py`
7. `normalize_illumina_vsn.R`
8. `build_gene_level_matrices.py`
9. `postcollapse_gene_qc.py`
10. `cross_platform_harmonization.py`

Supporting audit:

- `probe_collapsing_audit.py`

The resulting harmonized gene-level matrices and phenotype information provide the inputs for discovery differential-expression analysis.

---

## 2. Differential-expression analysis

Directory:

`scripts/02_differential_expression/`

Purpose: perform primary limma differential-expression analysis and associated sample, contrast, alignment, and sensitivity checks.

Recommended sequence:

1. `validate_de_inputs.py`
2. `build_de_sample_inventory.py`
3. `verify_de_sample_alignment.py`
4. `run_primary_limma_DE.R`
5. `run_DE_sensitivity_analyses.R`
6. `audit_de_contrasts.py`
7. `audit_primary_DE_results.py`

Primary discovery contrasts include influenza versus control, RSV versus control, and influenza versus RSV where applicable.

---

## 3. Cross-platform replication and frozen 170-gene architecture

Directory:

`scripts/03_cross_platform_replication/`

First run:

`analyze_cross_platform_RSV_replication.py`

This evaluates replication of the RSV transcriptional response across the two GSE38900 platforms.

Gene-architecture scripts are located in:

`scripts/03_cross_platform_replication/gene_architecture/`

Recommended sequence:

1. `build_gene_evidence_master.py`
2. `classify_pathogen_response_patterns.py`
3. `run_signature_sensitivity.py`
4. `audit_signature_concordance.py`

The authoritative gene-level evidence master is:

`DS001_GSE38900_gene_evidence_master.tsv`

The frozen architecture contains 170 genes:

- 130 shared-core genes
- 37 influenza-amplified genes
- 3 RSV-amplified genes

Downstream validation analyses should consume this frozen architecture rather than independently reselecting genes.

---

## 4. Pathway analysis

Directory:

`scripts/04_pathway_analysis/`

Purpose: generate ranked gene lists, freeze pathway resources, perform Hallmark and Reactome enrichment analysis, assess pathway robustness, and summarize preregistered biological domains.

Representative sequence:

1. `validate_pathway_inputs.py`
2. `freeze_msigdb_gene_sets.R`
3. `create_ranked_gene_lists.py`
4. `run_hallmark_fgsea.R`
5. `run_reactome_fgsea.R`
6. `compare_hallmark_contrasts.R`
7. `compare_reactome_contrasts.R`
8. `summarize_preregistered_domains.R`

Supporting analyses:

- `audit_gene_set_coverage.R`
- `audit_pathway_robustness.R`
- `fix_leading_edge_audit.R`
- `plot_session17_pathway_results.R`

---

## 5. Cellular-composition analysis

Directory:

`scripts/05_cellular_composition/`

Purpose: estimate immune-cell composition using MCP-counter, compare disease groups, replicate RSV composition findings, and evaluate whether pathway signals remain after accounting for estimated cellular composition.

Representative sequence:

1. `prepare_cell_composition_input.R`
2. `run_mcp_counter.R`
3. `qc_mcp_counter_scores.R`
4. `build_cell_composition_cohorts.R`
5. `test_mcp_counter_disease_groups.R`
6. `build_rsv_composition_replication.R`
7. `score_representative_pathways.R`
8. `test_composition_adjusted_pathway_effects.R`
9. `test_interferon_gene_composition_robustness.R`
10. `classify_pathway_composition_sensitivity.R`

---

## 6. Regulatory analysis

Directory:

`scripts/06_regulatory_analysis/`

Purpose: infer transcription-factor activity using the frozen CollecTRI network and ULM, test differential TF activity, map regulatory signals onto the frozen gene architecture, and evaluate regulatory evidence in external datasets.

Core discovery sequence:

1. `freeze_collectri_network.py`
2. `audit_collectri_coverage.py`
3. `run_GPL10558_collectri_ulm.py`
4. `run_GPL10558_TF_activity_limma.R`
5. `run_GPL6884_collectri_ulm.py`
6. `run_GPL6884_TF_activity_limma.R`
7. `map_regulators_to_frozen_gene_architecture.py`
8. `build_regulatory_architecture_classification.py`

Additional integration:

- `build_regulatory_cellular_support_matrix.py`

GSE155925 regulatory validation:

- `audit_GSE155925_regulatory_validation_inputs.py`
- `run_GSE155925_collectri_ulm.py`
- `run_GSE155925_TF_activity_limma.R`
- `run_GSE155925_TF_activity_unadjusted_sensitivity.R`
- `diagnose_GSE155925_regulatory_discrepancy.py`

GSE283746 regulatory analysis:

- `prepare_GSE283746_regulatory_logCPM.R`
- `run_GSE283746_collectri_ulm.py`
- `run_GSE283746_TF_activity_limma.R`

---

## 7. Independent RNA-seq and multiseason validation

Directory:

`scripts/07_rnaseq_validation/`

### 7a. GSE155925 RNA-seq validation

Purpose: construct the prespecified RSV-only versus virus-negative cohort, perform adjusted DESeq2 analysis, evaluate the frozen 170-gene architecture, assess model sensitivity, and test preregistered pathways.

Representative sequence:

1. `parse_GSE155925_metadata.py`
2. `build_GSE155925_eligibility.py`
3. `qc_GSE155925_counts.py`
4. `pca_GSE155925_qc.py`
5. `run_GSE155925_DESeq2_primary.R`
6. `build_GSE155925_frozen170_validation.py`
7. `build_GSE155925_frozen_pathway_membership.py`
8. `run_GSE155925_preregistered_pathway_GSEA.R`

Sensitivity and audit scripts:

- `run_GSE155925_age_sensitivity.R`
- `run_GSE155925_model_sensitivity.R`
- `compare_GSE155925_age_sensitivity_frozen170.py`
- `compare_GSE155925_models_frozen170.py`
- `reconcile_GSE155925_M0_direction.py`
- `audit_GSE155925_direction_reversal.py`
- `audit_GSE155925_pathway_validation.py`

### 7b. Multiseason influenza validation

Directory:

`scripts/07_rnaseq_validation/multiseason_SIG/`

Recommended sequence:

1. `build_SIG208_metadata.py`
2. `build_original_SIG_expression_metadata_map.py`
3. `parse_original_SIG_matrices.py`
4. `run_SIG_seasonal_limma.R`
5. `validate_SIG_seasons_frozen_architecture.py`

---

## 8. Influenza proteomic validation

Directory:

`scripts/08_proteomics/`

Purpose: test whether the frozen transcriptomic architecture is supported by independent influenza SomaScan proteomic evidence.

Recommended sequence:

1. `build_Influenza_P01_analysis_inputs.py`
2. `run_Influenza_P01_limma_v1.0.R`
3. `validate_Influenza_P01_frozen170_v1.0.py`
4. `run_Influenza_P01_frozen_pathway_validation_v1.0.R`
5. `validate_Influenza_P01_frozen37_regulatory_programs_v1.0.py`

This module evaluates gene-level transcript-protein concordance, pathway-level convergence, and regulatory-program support.

---

## 9. CRISPR functional evidence

Directory:

`scripts/09_crispr/`

Important provenance limitation:

The original Session26 executable script that generated the final frozen CRISPR downstream result tables was not recovered.

The public reproducibility layer therefore preserves:

- frozen dataset and input manifests
- frozen analysis specifications
- score rules
- producer/provenance documentation
- final compact result snapshots

These artifacts preserve the manuscript-level CRISPR evidence and analysis definition but must not be represented as the original executable CRISPR workflow.

Any future reconstructed implementation should be explicitly labeled:

`RECONSTRUCTED_REPRODUCTION_SCRIPT — NOT ORIGINAL EXECUTION SCRIPT`

---

## 10. Single-cell analyses

Directory:

`scripts/10_single_cell/`

### 10a. GSE283746

Representative sequence:

1. `audit_GSE283746_all_barcode_mapping.py`
2. `build_GSE283746_primary_pseudobulk.py`
3. `run_GSE283746_primary_edgeR.R`
4. `validate_GSE283746_frozen170.py`
5. `audit_GSE283746_ISG_state_coverage.py`
6. `build_GSE283746_NaiveCD4_ISG_pseudobulk.py`
7. `run_GSE283746_NaiveCD4_ISG_edgeR.R`
8. `compare_GSE283746_NaiveCD4_ISG_frozen170.py`
9. `run_GSE283746_ISG_state_abundance.R`

### 10b. GSE149689 supportive localization

Directory:

`scripts/10_single_cell/GSE149689_supportive_localization/`

Recommended sequence:

1. `run_Session38C_full_library_sizes_v1.py`
2. `run_Session38C_extract_pseudobulk_v1.1.py`
3. `run_Session38C_detection_robustness_v1.py`
4. `run_Session38C_localization_inference_v1.py`

---

## 11. Multilayer integration

Directory:

`scripts/11_multilayer_integration/`

Purpose: integrate discovery transcriptional evidence with proteomic, CRISPR, multiseason, regulatory, and single-cell evidence.

Recommended sequence:

1. `run_Session38A_evidence_weighted_prioritization_v1.1.py`
2. `run_Session38A_weight_sensitivity_v1.py`
3. `run_Session38B_molecular_functional_convergence_v1.1.py`
4. `run_Session38D_crossseason_stability_v1.1.py`

These analyses depend on frozen outputs from the relevant upstream evidence layers.

---

## 12. Manuscript figures

Directory:

`scripts/12_figures/`

Canonical main-figure generators are located in:

`scripts/12_figures/main/`

and generate Figures 1-8.

Canonical supplementary-figure generators are located in:

`scripts/12_figures/supplementary/`

and generate Supplementary Figures S1-S11.

The authoritative mapping is:

`workflow/FIGURE_TO_SCRIPT_MAP.tsv`

Historical Session39 Supplementary Figures S11-S13 were subsequently renumbered as final manuscript Supplementary Figures S9-S11.

The final manuscript supplementary-figure series therefore ends at S11.

---

## 13. Supplementary tables

Directory:

`scripts/13_tables/`

Purpose: generate, quality-control, register, and stage supplementary-table material.

The authoritative table provenance mapping is:

`workflow/TABLE_TO_SCRIPT_MAP.tsv`

---

## Reproducibility metadata

The repository script-layer integrity manifest is:

`provenance/manifests/SCRIPT_SHA256SUMS_v1.0.txt`

The current `scripts/` reproducibility layer contains 142 files.

This total includes executable Python, R, and shell scripts together with the frozen CRISPR provenance/specification and result-snapshot artifacts described above.

The one-to-one repository artifact registry is:

`workflow/SCRIPT_TO_OUTPUT_MAP.tsv`

The canonical manuscript figure-generator registry is:

`workflow/FIGURE_TO_SCRIPT_MAP.tsv`

The supplementary-table provenance registry is:

`workflow/TABLE_TO_SCRIPT_MAP.tsv`

---

## General execution principle

Analytical selection and biological classification were frozen upstream before downstream validation where required by the study design.

Downstream validation modules should therefore consume the frozen upstream definitions rather than recomputing or reselecting candidate genes, pathways, or regulators using validation datasets themselves.

This preserves the distinction among discovery, validation, sensitivity analysis, and supportive evidence.
