#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(limma)
})

# ============================================================
# SESSION 16 — PRIMARY DIFFERENTIAL EXPRESSION
# Dataset: DS001_GSE38900
# Method: limma empirical-Bayes linear models
# ============================================================

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

meta_file <- file.path(
    root,
    "data/processed/DS001_GSE38900",
    "DS001_GSE38900_master_sample_metadata_harmonized.tsv"
)

mapping_file <- file.path(
    root,
    "results/DS001_GSE38900/differential_expression/tables",
    "DS001_GSE38900_expression_metadata_sample_mapping.tsv"
)

expr_files <- list(
    GPL10558 = file.path(
        root,
        "data/processed/DS001_GSE38900/analysis_ready",
        "GPL10558_gene_expression.tsv.gz"
    ),
    GPL6884 = file.path(
        root,
        "data/processed/DS001_GSE38900/analysis_ready",
        "GPL6884_gene_expression.tsv.gz"
    )
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/differential_expression"
)

table_dir <- file.path(outdir, "tables")
design_dir <- file.path(outdir, "design")

dir.create(table_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(design_dir, recursive=TRUE, showWarnings=FALSE)

cat("============================================================\n")
cat("SESSION 16 — PRIMARY LIMMA DIFFERENTIAL EXPRESSION\n")
cat("============================================================\n\n")

cat("limma version:", as.character(packageVersion("limma")), "\n")
cat("R version:", R.version.string, "\n\n")

# ------------------------------------------------------------
# Read metadata
# ------------------------------------------------------------

meta <- read.delim(
    meta_file,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

mapping <- read.delim(
    mapping_file,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

# ------------------------------------------------------------
# Utility function
# ------------------------------------------------------------

save_full_results <- function(
    fit,
    coef_name,
    annotation,
    filename,
    contrast_label
) {

    tt <- topTable(
        fit,
        coef=coef_name,
        number=Inf,
        adjust.method="BH",
        sort.by="P"
    )

    tt$gene_symbol <- rownames(tt)

    tt <- merge(
        tt,
        annotation,
        by="gene_symbol",
        all.x=TRUE,
        sort=FALSE
    )

    tt$contrast <- contrast_label

    # Reorder key columns
    preferred <- c(
        "gene_symbol",
        "entrez_gene_id",
        "probe_id",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B",
        "contrast"
    )

    remaining <- setdiff(
        colnames(tt),
        preferred
    )

    tt <- tt[, c(
        preferred[preferred %in% colnames(tt)],
        remaining
    )]

    write.table(
        tt,
        file=filename,
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    return(tt)
}

summarize_results <- function(
    table,
    contrast,
    platform
) {

    data.frame(
        platform=platform,
        contrast=contrast,

        genes_tested=nrow(table),

        FDR_0_05=sum(
            table$adj.P.Val < 0.05,
            na.rm=TRUE
        ),

        FDR_0_05_logFC_ge_0_5=sum(
            table$adj.P.Val < 0.05 &
            abs(table$logFC) >= 0.5,
            na.rm=TRUE
        ),

        FDR_0_05_logFC_ge_1=sum(
            table$adj.P.Val < 0.05 &
            abs(table$logFC) >= 1,
            na.rm=TRUE
        ),

        up_FDR_0_05=sum(
            table$adj.P.Val < 0.05 &
            table$logFC > 0,
            na.rm=TRUE
        ),

        down_FDR_0_05=sum(
            table$adj.P.Val < 0.05 &
            table$logFC < 0,
            na.rm=TRUE
        ),

        stringsAsFactors=FALSE
    )
}

summary_list <- list()

# ============================================================
# GPL10558
# RSV acute vs healthy control
#
# Primary model:
# expression ~ group
#
# Reason:
# age units are not harmonized on GPL10558.
# Sex adjustment will be treated separately as sensitivity.
# ============================================================

cat("\n============================================================\n")
cat("GPL10558 PRIMARY ANALYSIS\n")
cat("============================================================\n")

expr10558 <- read.delim(
    gzfile(expr_files$GPL10558),
    stringsAsFactors=FALSE,
    check.names=FALSE
)

annotation10558 <- expr10558[, c(
    "gene_symbol",
    "entrez_gene_id",
    "probe_id"
)]

rownames(expr10558) <- expr10558$gene_symbol

expr10558 <- expr10558[, !colnames(expr10558) %in% c(
    "gene_symbol",
    "entrez_gene_id",
    "probe_id"
)]

map10558 <- mapping[
    mapping$platform_id == "GPL10558",
]

# Ensure expression-column order
map10558 <- map10558[
    match(
        colnames(expr10558),
        map10558$expression_sample_id
    ),
]

stopifnot(
    all(
        colnames(expr10558) ==
        map10558$expression_sample_id
    )
)

keep10558 <- map10558$harmonized_group %in% c(
    "Healthy control",
    "RSV acute"
)

expr10558_sub <- as.matrix(
    expr10558[, keep10558]
)

mode(expr10558_sub) <- "numeric"

meta10558 <- map10558[keep10558,]

group10558 <- factor(
    meta10558$harmonized_group,
    levels=c(
        "Healthy control",
        "RSV acute"
    )
)

design10558 <- model.matrix(
    ~ 0 + group10558
)

colnames(design10558) <- c(
    "Healthy_control",
    "RSV_acute"
)

cat("\nGPL10558 sample counts:\n")
print(table(group10558))

cat("\nGPL10558 design matrix dimensions:\n")
print(dim(design10558))

cat("\nGPL10558 design rank:\n")
print(qr(design10558)$rank)

if (
    qr(design10558)$rank != ncol(design10558)
) {
    stop("GPL10558 design matrix is not full rank.")
}

write.table(
    design10558,
    file=file.path(
        design_dir,
        "GPL10558_primary_design_matrix.tsv"
    ),
    sep="\t",
    quote=FALSE,
    col.names=NA
)

contrast10558 <- makeContrasts(
    RSV_acute_vs_control =
        RSV_acute - Healthy_control,

    levels=design10558
)

fit10558 <- lmFit(
    expr10558_sub,
    design10558
)

fit10558 <- contrasts.fit(
    fit10558,
    contrast10558
)

fit10558 <- eBayes(
    fit10558
)

result10558 <- save_full_results(
    fit=fit10558,
    coef_name="RSV_acute_vs_control",
    annotation=annotation10558,
    filename=file.path(
        table_dir,
        "GPL10558_RSV_acute_vs_control_full_DE.tsv"
    ),
    contrast_label="RSV acute vs healthy control"
)

summary_list[[length(summary_list)+1]] <-
    summarize_results(
        result10558,
        "RSV acute vs healthy control",
        "GPL10558"
    )

# ============================================================
# GPL6884
#
# Primary model:
# expression ~ group + age_months + sex
#
# Groups retained:
# Healthy control
# RSV acute
# Influenza A acute
# HRV acute
# RSV recovery
#
# Primary contrasts:
# 1. RSV acute vs healthy control
# 2. Influenza A acute vs healthy control
# 3. Influenza A acute vs RSV acute
# ============================================================

cat("\n============================================================\n")
cat("GPL6884 PRIMARY ANALYSIS\n")
cat("============================================================\n")

expr6884 <- read.delim(
    gzfile(expr_files$GPL6884),
    stringsAsFactors=FALSE,
    check.names=FALSE
)

annotation6884 <- expr6884[, c(
    "gene_symbol",
    "entrez_gene_id",
    "probe_id"
)]

rownames(expr6884) <- expr6884$gene_symbol

expr6884 <- expr6884[, !colnames(expr6884) %in% c(
    "gene_symbol",
    "entrez_gene_id",
    "probe_id"
)]

map6884 <- mapping[
    mapping$platform_id == "GPL6884",
]

map6884 <- map6884[
    match(
        colnames(expr6884),
        map6884$expression_sample_id
    ),
]

stopifnot(
    all(
        colnames(expr6884) ==
        map6884$expression_sample_id
    )
)

expr6884_mat <- as.matrix(expr6884)
mode(expr6884_mat) <- "numeric"

map6884$age_months_harmonized <- as.numeric(
    map6884$age_months_harmonized
)

if (
    any(is.na(map6884$age_months_harmonized))
) {
    stop(
        "GPL6884 contains missing harmonized ages."
    )
}

if (
    any(is.na(map6884$sex_standardized))
) {
    stop(
        "GPL6884 contains missing standardized sex."
    )
}

group6884 <- factor(
    map6884$harmonized_group,
    levels=c(
        "Healthy control",
        "RSV acute",
        "Influenza A acute",
        "HRV acute",
        "RSV recovery"
    )
)

sex6884 <- factor(
    map6884$sex_standardized,
    levels=c(
        "Female",
        "Male"
    )
)

age6884 <- map6884$age_months_harmonized

design6884 <- model.matrix(
    ~ 0 + group6884 + age6884 + sex6884
)

colnames(design6884) <- c(
    "Healthy_control",
    "RSV_acute",
    "Influenza_A_acute",
    "HRV_acute",
    "RSV_recovery",
    "age_months",
    "sex_Male"
)

cat("\nGPL6884 sample counts:\n")
print(table(group6884))

cat("\nGPL6884 sex counts:\n")
print(table(sex6884))

cat("\nGPL6884 age summary:\n")
print(summary(age6884))

cat("\nGPL6884 design dimensions:\n")
print(dim(design6884))

cat("\nGPL6884 design rank:\n")
print(qr(design6884)$rank)

if (
    qr(design6884)$rank != ncol(design6884)
) {
    stop("GPL6884 design matrix is not full rank.")
}

write.table(
    design6884,
    file=file.path(
        design_dir,
        "GPL6884_primary_design_matrix.tsv"
    ),
    sep="\t",
    quote=FALSE,
    col.names=NA
)

contrasts6884 <- makeContrasts(

    RSV_acute_vs_control =
        RSV_acute - Healthy_control,

    Influenza_A_vs_control =
        Influenza_A_acute - Healthy_control,

    Influenza_A_vs_RSV_acute =
        Influenza_A_acute - RSV_acute,

    levels=design6884
)

fit6884 <- lmFit(
    expr6884_mat,
    design6884
)

fit6884 <- contrasts.fit(
    fit6884,
    contrasts6884
)

fit6884 <- eBayes(
    fit6884
)

# ------------------------------------------------------------
# Contrast 1
# ------------------------------------------------------------

rsv6884 <- save_full_results(
    fit=fit6884,
    coef_name="RSV_acute_vs_control",
    annotation=annotation6884,
    filename=file.path(
        table_dir,
        "GPL6884_RSV_acute_vs_control_full_DE.tsv"
    ),
    contrast_label="RSV acute vs healthy control"
)

summary_list[[length(summary_list)+1]] <-
    summarize_results(
        rsv6884,
        "RSV acute vs healthy control",
        "GPL6884"
    )

# ------------------------------------------------------------
# Contrast 2
# ------------------------------------------------------------

flu6884 <- save_full_results(
    fit=fit6884,
    coef_name="Influenza_A_vs_control",
    annotation=annotation6884,
    filename=file.path(
        table_dir,
        "GPL6884_InfluenzaA_vs_control_full_DE.tsv"
    ),
    contrast_label="Influenza A acute vs healthy control"
)

summary_list[[length(summary_list)+1]] <-
    summarize_results(
        flu6884,
        "Influenza A acute vs healthy control",
        "GPL6884"
    )

# ------------------------------------------------------------
# Contrast 3
# ------------------------------------------------------------

flu_rsv6884 <- save_full_results(
    fit=fit6884,
    coef_name="Influenza_A_vs_RSV_acute",
    annotation=annotation6884,
    filename=file.path(
        table_dir,
        "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv"
    ),
    contrast_label="Influenza A acute vs RSV acute"
)

summary_list[[length(summary_list)+1]] <-
    summarize_results(
        flu_rsv6884,
        "Influenza A acute vs RSV acute",
        "GPL6884"
    )

# ============================================================
# Save master summary
# ============================================================

summary_df <- do.call(
    rbind,
    summary_list
)

summary_file <- file.path(
    table_dir,
    "DS001_GSE38900_primary_DE_summary.tsv"
)

write.table(
    summary_df,
    file=summary_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n============================================================\n")
cat("PRIMARY DE SUMMARY\n")
cat("============================================================\n\n")

print(summary_df)

cat("\nSaved summary:\n")
cat(summary_file, "\n")

cat("\nPrimary DE analysis completed successfully.\n")
