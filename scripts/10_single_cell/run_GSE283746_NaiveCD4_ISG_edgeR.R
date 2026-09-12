suppressPackageStartupMessages({
    library(edgeR)
    library(readr)
    library(dplyr)
})

count_dir <- paste0(
    "results/scrnaseq_validation/GSE283746/",
    "isg_state/pseudobulk/"
)

meta_file <- paste0(
    "results/scrnaseq_validation/GSE283746/pseudobulk/",
    "GSE283746_primary_pseudobulk_sample_metadata_FROZEN_v1.0.tsv"
)

outdir <- paste0(
    "results/scrnaseq_validation/GSE283746/",
    "isg_state/DE"
)

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# Frozen participant metadata
# ------------------------------------------------------------

meta0 <- read_tsv(
    meta_file,
    show_col_types = FALSE
) %>%
    distinct(Sample, .keep_all = TRUE)

cat("=== PARTICIPANT METADATA ===\n")
cat("Unique participants:", nrow(meta0), "\n")
cat("Columns:\n")
print(names(meta0))

required_meta <- c(
    "Sample",
    "Combined Condition",
    "Age",
    "Sex"
)

missing_meta <- setdiff(
    required_meta,
    names(meta0)
)

if (length(missing_meta) > 0) {
    stop(
        "Missing metadata columns: ",
        paste(missing_meta, collapse = ", ")
    )
}

# ------------------------------------------------------------
# Frozen states
# ------------------------------------------------------------

states <- list(

    ISG_Naive_CD4_T =
        "GSE283746_ISG_Naive_CD4_T_counts.tsv.gz",

    Reference_Naive_CD4_T =
        "GSE283746_Reference_Naive_CD4_T_counts.tsv.gz"
)

summary_rows <- list()

# ------------------------------------------------------------
# Analyze each state independently
# ------------------------------------------------------------

for (state in names(states)) {

    cat("\n========================================\n")
    cat(state, "\n")
    cat("========================================\n")

    count_file <- file.path(
        count_dir,
        states[[state]]
    )

    if (!file.exists(count_file)) {
        stop(
            "Missing count file: ",
            count_file
        )
    }

    x <- read_tsv(
        count_file,
        show_col_types = FALSE
    )

    required_count_cols <- c(
        "feature_id",
        "gene_symbol"
    )

    if (!all(required_count_cols %in% names(x))) {
        stop(
            state,
            ": feature_id/gene_symbol columns missing"
        )
    }

    sample_cols <- setdiff(
        names(x),
        required_count_cols
    )

    cat(
        "Count-matrix samples:",
        length(sample_cols),
        "\n"
    )

    # --------------------------------------------------------
    # Numeric count matrix
    # --------------------------------------------------------

    counts <- as.matrix(
        x[, sample_cols]
    )

    storage.mode(counts) <- "integer"

    gene_key <- make.unique(
        paste0(
            x$gene_symbol,
            ":",
            x$feature_id
        )
    )

    rownames(counts) <- gene_key

    # --------------------------------------------------------
    # Match frozen participant metadata
    # --------------------------------------------------------

    idx <- match(
        sample_cols,
        meta0$Sample
    )

    if (any(is.na(idx))) {

        cat("\nSamples missing metadata:\n")
        print(
            sample_cols[is.na(idx)]
        )

        stop(
            state,
            ": one or more samples missing frozen metadata"
        )
    }

    meta <- as.data.frame(
        meta0[idx, , drop = FALSE]
    )

    stopifnot(
        identical(
            sample_cols,
            as.character(meta$Sample)
        )
    )

    # --------------------------------------------------------
    # Disease / covariates
    # --------------------------------------------------------

    meta$RSV_status <- factor(
        meta$`Combined Condition`,
        levels = c(
            "Healthy",
            "RSV"
        )
    )

    meta$sex <- factor(
        meta$Sex,
        levels = c(
            "F",
            "M"
        )
    )

    meta$Age <- as.numeric(
        meta$Age
    )

    if (
        anyNA(meta$RSV_status) ||
        anyNA(meta$sex) ||
        anyNA(meta$Age)
    ) {
        stop(
            state,
            ": missing disease, age, or sex covariate"
        )
    }

    # Standardize age within state-specific analysis set.
    meta$age_z <- as.numeric(
        scale(meta$Age)
    )

    cat("\nGroup counts:\n")
    print(
        table(meta$RSV_status)
    )

    cat("\nSex counts:\n")
    print(
        table(meta$sex)
    )

    cat("\nAge by group:\n")
    print(
        tapply(
            meta$Age,
            meta$RSV_status,
            summary
        )
    )

    # --------------------------------------------------------
    # Primary model
    # --------------------------------------------------------

    design <- model.matrix(
        ~ age_z + sex + RSV_status,
        data = meta
    )

    cat(
        "\nDesign parameters:",
        ncol(design),
        " rank:",
        qr(design)$rank,
        " full rank:",
        qr(design)$rank == ncol(design),
        "\n"
    )

    if (
        qr(design)$rank !=
        ncol(design)
    ) {
        stop(
            state,
            ": primary design rank deficient"
        )
    }

    # --------------------------------------------------------
    # edgeR
    # --------------------------------------------------------

    y <- DGEList(
        counts = counts
    )

    genes_before <- nrow(y)

    keep <- filterByExpr(
        y,
        design = design
    )

    genes_after <- sum(keep)

    cat(
        "Genes before filter:",
        genes_before,
        "\n"
    )

    cat(
        "Genes after filter:",
        genes_after,
        "\n"
    )

    y <- y[
        keep,
        ,
        keep.lib.sizes = FALSE
    ]

    y <- calcNormFactors(y)

    y <- estimateDisp(
        y,
        design,
        robust = TRUE
    )

    fit <- glmQLFit(
        y,
        design,
        robust = TRUE
    )

    coef_name <- "RSV_statusRSV"

    if (
        !(coef_name %in%
          colnames(design))
    ) {
        stop(
            state,
            ": RSV coefficient not found. Design columns: ",
            paste(
                colnames(design),
                collapse = ", "
            )
        )
    }

    qlf <- glmQLFTest(
        fit,
        coef = which(
            colnames(design) ==
            coef_name
        )
    )

    tt <- topTags(
        qlf,
        n = Inf,
        sort.by = "PValue"
    )$table

    tt$gene_key <- rownames(tt)

    # Re-attach annotation using original key.
    annotation <- data.frame(
        gene_key = gene_key,
        feature_id = x$feature_id,
        gene_symbol = x$gene_symbol,
        stringsAsFactors = FALSE
    )

    tt <- tt %>%
        tibble::rownames_to_column(
            "tmp_rowname"
        ) %>%
        mutate(
            gene_key =
                tmp_rowname
        ) %>%
        select(
            -tmp_rowname
        ) %>%
        left_join(
            annotation,
            by = "gene_key"
        ) %>%
        select(
            feature_id,
            gene_symbol,
            logFC,
            logCPM,
            F,
            PValue,
            FDR
        )

    outfile <- file.path(
        outdir,
        paste0(
            "GSE283746_",
            state,
            "_RSV_vs_Healthy_edgeR.tsv"
        )
    )

    write_tsv(
        tt,
        outfile
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    n_fdr <- sum(
        tt$FDR < 0.05,
        na.rm = TRUE
    )

    n_fdr_lfc05 <- sum(
        tt$FDR < 0.05 &
        abs(tt$logFC) >= 0.5,
        na.rm = TRUE
    )

    n_fdr_lfc1 <- sum(
        tt$FDR < 0.05 &
        abs(tt$logFC) >= 1,
        na.rm = TRUE
    )

    cat(
        "\nFDR<0.05:",
        n_fdr,
        "\n"
    )

    cat(
        "FDR<0.05 & |logFC|>=0.5:",
        n_fdr_lfc05,
        "\n"
    )

    cat(
        "FDR<0.05 & |logFC|>=1:",
        n_fdr_lfc1,
        "\n"
    )

    summary_rows[[state]] <- data.frame(

        state_id =
            state,

        total_samples =
            nrow(meta),

        RSV_samples =
            sum(
                meta$RSV_status ==
                "RSV"
            ),

        Healthy_samples =
            sum(
                meta$RSV_status ==
                "Healthy"
            ),

        genes_tested =
            genes_after,

        FDR05 =
            n_fdr,

        FDR05_absLFC05 =
            n_fdr_lfc05,

        FDR05_absLFC1 =
            n_fdr_lfc1,

        design_parameters =
            ncol(design),

        design_rank =
            qr(design)$rank,

        stringsAsFactors = FALSE
    )
}

summary_df <- bind_rows(
    summary_rows
)

summary_file <- file.path(
    outdir,
    "GSE283746_NaiveCD4_ISG_state_edgeR_summary.tsv"
)

write_tsv(
    summary_df,
    summary_file
)

cat("\n========================================\n")
cat("FINAL SUMMARY\n")
cat("========================================\n")

print(summary_df)

cat("\nWritten:\n")
cat(summary_file, "\n")

cat(
    "\nState-specific edgeR analyses complete.\n"
)
