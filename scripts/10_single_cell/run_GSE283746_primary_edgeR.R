suppressPackageStartupMessages({
    library(edgeR)
    library(readr)
    library(dplyr)
    library(stringr)
})

count_dir <- "results/scrnaseq_validation/GSE283746/pseudobulk/counts"

meta_file <- paste0(
    "results/scrnaseq_validation/GSE283746/pseudobulk/",
    "GSE283746_primary_pseudobulk_sample_metadata_FROZEN_v1.0.tsv"
)

outdir <- "results/scrnaseq_validation/GSE283746/DE"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

meta_all <- read_tsv(
    meta_file,
    show_col_types=FALSE
)

# ------------------------------------------------------------
# Frozen cell-type -> filename mapping
# ------------------------------------------------------------

celltypes <- c(
    "Naive CD4+ T",
    "Naive CD8+ T",
    "TrB",
    "B",
    "NK",
    "Memory CD4+ T",
    "CD14+ Monocyte",
    "Cytotoxic T",
    "Treg",
    "Memory B",
    "CD16+ Monocyte"
)

safe_name <- function(x) {
    x |>
        str_replace_all("\\+", "plus") |>
        str_replace_all("/", "_") |>
        str_replace_all(" ", "_")
}

summary_rows <- list()

for (ct in celltypes) {

    cat("\n========================================\n")
    cat(ct, "\n")
    cat("========================================\n")

    safe <- safe_name(ct)

    count_file <- file.path(
        count_dir,
        paste0(
            "GSE283746_pseudobulk_",
            safe,
            "_counts.tsv.gz"
        )
    )

    x <- read_tsv(
        count_file,
        show_col_types=FALSE
    )

    feature_id <- x$feature_id
    gene_symbol <- x$gene_symbol

    count_cols <- setdiff(
        colnames(x),
        c("feature_id", "gene_symbol")
    )

    counts <- as.matrix(
        x[, count_cols]
    )

    storage.mode(counts) <- "integer"

    rownames(counts) <- make.unique(
        paste0(
            gene_symbol,
            ":",
            feature_id
        )
    )

    # --------------------------------------------------------
    # Metadata matching
    # --------------------------------------------------------

    meta <- meta_all |>
        filter(GroupedAnnotation == ct) |>
        distinct(Sample, .keep_all=TRUE)

    meta <- as.data.frame(meta)

    idx <- match(
        colnames(counts),
        meta$Sample
    )

    if (any(is.na(idx))) {
        stop(
            paste(
                ct,
                "count-matrix samples missing metadata"
            )
        )
    }

    meta <- meta[idx, , drop=FALSE]

    stopifnot(
        identical(
            colnames(counts),
            meta$Sample
        )
    )

    meta$RSV_status <- factor(
        meta$`Combined Condition`,
        levels=c("Healthy", "RSV")
    )

    meta$sex <- factor(
        meta$Sex,
        levels=c("F", "M")
    )

    meta$age_z <- as.numeric(
        scale(meta$Age)
    )

    cat("Samples:", nrow(meta), "\n")
    print(table(meta$RSV_status))
    print(table(meta$sex))

    # --------------------------------------------------------
    # edgeR object and expression filtering
    # --------------------------------------------------------

    y <- DGEList(
        counts=counts
    )

    design_primary <- model.matrix(
        ~ age_z + sex + RSV_status,
        data=meta
    )

    keep <- filterByExpr(
        y,
        design=design_primary
    )

    cat(
        "Genes before filter:",
        nrow(y),
        "\n"
    )

    cat(
        "Genes after filter:",
        sum(keep),
        "\n"
    )

    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)

    # --------------------------------------------------------
    # Primary model
    # --------------------------------------------------------

    qr_primary <- qr(design_primary)

    cat(
        "Primary parameters:",
        ncol(design_primary),
        " rank:",
        qr_primary$rank,
        " full rank:",
        qr_primary$rank == ncol(design_primary),
        "\n"
    )

    if (qr_primary$rank != ncol(design_primary)) {
        stop(
            paste(
                ct,
                "primary design is rank deficient"
            )
        )
    }

    y1 <- estimateDisp(
        y,
        design_primary,
        robust=TRUE
    )

    fit1 <- glmQLFit(
        y1,
        design_primary,
        robust=TRUE
    )

    coef_name <- grep(
        "^RSV_statusRSV$",
        colnames(design_primary),
        value=TRUE
    )

    if (length(coef_name) != 1) {
        stop(
            paste(
                ct,
                "could not identify RSV coefficient"
            )
        )
    }

    qlf1 <- glmQLFTest(
        fit1,
        coef=which(
            colnames(design_primary) == coef_name
        )
    )

    tt1 <- topTags(
        qlf1,
        n=Inf,
        sort.by="PValue"
    )$table

    tt1$gene_key <- rownames(tt1)

    split_primary <- str_split_fixed(
        tt1$gene_key,
        ":",
        2
    )

    tt1$gene_symbol <- split_primary[,1]
    tt1$feature_id <- split_primary[,2]

    tt1 <- tt1 |>
        select(
            feature_id,
            gene_symbol,
            logFC,
            logCPM,
            F,
            PValue,
            FDR
        )

    primary_file <- file.path(
        outdir,
        paste0(
            "GSE283746_",
            safe,
            "_RSV_vs_Healthy_edgeR_primary.tsv"
        )
    )

    write_tsv(
        tt1,
        primary_file
    )

    # --------------------------------------------------------
    # Sensitivity model: disease only
    # --------------------------------------------------------

    design_sens <- model.matrix(
        ~ RSV_status,
        data=meta
    )

    qr_sens <- qr(design_sens)

    if (qr_sens$rank != ncol(design_sens)) {
        stop(
            paste(
                ct,
                "sensitivity design rank deficient"
            )
        )
    }

    ys <- DGEList(
        counts=counts[keep, , drop=FALSE]
    )

    ys <- calcNormFactors(ys)

    ys <- estimateDisp(
        ys,
        design_sens,
        robust=TRUE
    )

    fits <- glmQLFit(
        ys,
        design_sens,
        robust=TRUE
    )

    qlfs <- glmQLFTest(
        fits,
        coef=2
    )

    tts <- topTags(
        qlfs,
        n=Inf,
        sort.by="PValue"
    )$table

    tts$gene_key <- rownames(tts)

    split_sens <- str_split_fixed(
        tts$gene_key,
        ":",
        2
    )

    tts$gene_symbol <- split_sens[,1]
    tts$feature_id <- split_sens[,2]

    tts <- tts |>
        select(
            feature_id,
            gene_symbol,
            logFC,
            logCPM,
            F,
            PValue,
            FDR
        )

    sens_file <- file.path(
        outdir,
        paste0(
            "GSE283746_",
            safe,
            "_RSV_vs_Healthy_edgeR_unadjusted.tsv"
        )
    )

    write_tsv(
        tts,
        sens_file
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    n_fdr <- sum(
        tt1$FDR < 0.05,
        na.rm=TRUE
    )

    n_fdr_lfc05 <- sum(
        tt1$FDR < 0.05 &
        abs(tt1$logFC) >= 0.5,
        na.rm=TRUE
    )

    n_fdr_lfc1 <- sum(
        tt1$FDR < 0.05 &
        abs(tt1$logFC) >= 1,
        na.rm=TRUE
    )

    summary_rows[[ct]] <- data.frame(
        GroupedAnnotation=ct,
        total_samples=nrow(meta),
        RSV_samples=sum(
            meta$RSV_status == "RSV"
        ),
        Healthy_samples=sum(
            meta$RSV_status == "Healthy"
        ),
        genes_tested=nrow(tt1),
        FDR05=n_fdr,
        FDR05_absLFC05=n_fdr_lfc05,
        FDR05_absLFC1=n_fdr_lfc1,
        primary_design_rank=qr_primary$rank,
        primary_design_parameters=
            ncol(design_primary),
        stringsAsFactors=FALSE
    )

    cat(
        "Primary FDR<0.05:",
        n_fdr,
        "\n"
    )

    cat(
        "Primary FDR<0.05 & |logFC|>=0.5:",
        n_fdr_lfc05,
        "\n"
    )

    cat(
        "Primary FDR<0.05 & |logFC|>=1:",
        n_fdr_lfc1,
        "\n"
    )
}

summary_df <- bind_rows(
    summary_rows
)

write_tsv(
    summary_df,
    file.path(
        outdir,
        "GSE283746_primary_edgeR_summary.tsv"
    )
)

cat("\n========================================\n")
cat("FINAL SUMMARY\n")
cat("========================================\n")

print(summary_df)

cat("\nAll primary pseudobulk analyses complete.\n")
