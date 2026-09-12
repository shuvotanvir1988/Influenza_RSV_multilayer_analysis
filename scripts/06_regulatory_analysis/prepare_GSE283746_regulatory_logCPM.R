suppressPackageStartupMessages({
    library(edgeR)
    library(readr)
    library(dplyr)
    library(stringr)
})

ROOT <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

count_dir <- file.path(
    ROOT,
    "results/scrnaseq_validation/GSE283746/pseudobulk/counts"
)

meta_file <- file.path(
    ROOT,
    "results/scrnaseq_validation/GSE283746/pseudobulk",
    "GSE283746_primary_pseudobulk_sample_metadata_FROZEN_v1.0.tsv"
)

outdir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/scrnaseq_validation/logCPM"
)

qcdir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/qc"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(qcdir, recursive=TRUE, showWarnings=FALSE)

meta_all <- read_tsv(meta_file, show_col_types=FALSE)

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

    if (!file.exists(count_file)) {
        stop(paste("Missing count file:", count_file))
    }

    x <- read_tsv(count_file, show_col_types=FALSE)

    count_cols <- setdiff(
        colnames(x),
        c("feature_id", "gene_symbol")
    )

    counts <- as.matrix(x[, count_cols])
    storage.mode(counts) <- "integer"

    rownames(counts) <- make.unique(
        paste0(x$gene_symbol, ":", x$feature_id)
    )

    meta <- meta_all |>
        filter(GroupedAnnotation == ct) |>
        distinct(Sample, .keep_all=TRUE)

    meta <- as.data.frame(meta)

    idx <- match(colnames(counts), meta$Sample)

    if (any(is.na(idx))) {
        stop(paste(ct, "count-matrix samples missing metadata"))
    }

    meta <- meta[idx, , drop=FALSE]

    stopifnot(identical(colnames(counts), meta$Sample))

    meta$RSV_status <- factor(
        meta$`Combined Condition`,
        levels=c("Healthy", "RSV")
    )

    meta$sex <- factor(
        meta$Sex,
        levels=c("F", "M")
    )

    meta$age_z <- as.numeric(scale(meta$Age))

    design_primary <- model.matrix(
        ~ age_z + sex + RSV_status,
        data=meta
    )

    if (qr(design_primary)$rank != ncol(design_primary)) {
        stop(paste(ct, "design matrix is not full rank"))
    }

    y <- DGEList(counts=counts)

    keep <- filterByExpr(
        y,
        design=design_primary
    )

    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)

    logcpm <- cpm(
        y,
        log=TRUE,
        prior.count=2
    )

    ann <- data.frame(
        feature_key=rownames(logcpm),
        stringsAsFactors=FALSE
    )

    ann$gene_symbol <- sub(":ENSG.*$", "", ann$feature_key)

    out <- data.frame(
        gene_symbol=ann$gene_symbol,
        logcpm,
        check.names=FALSE
    )

    # Deterministic duplicate-symbol collapse by mean.
    out <- aggregate(
        . ~ gene_symbol,
        data=out,
        FUN=mean
    )

    out_file <- file.path(
        outdir,
        paste0(
            "GSE283746_",
            safe,
            "_TMM_logCPM_for_ULM_v1.0.tsv.gz"
        )
    )

    write.table(
        out,
        gzfile(out_file),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    summary_rows[[length(summary_rows)+1]] <- data.frame(
        cell_type=ct,
        safe_name=safe,
        samples=ncol(counts),
        RSV_samples=sum(meta$RSV_status == "RSV"),
        Healthy_samples=sum(meta$RSV_status == "Healthy"),
        genes_before_filter=nrow(counts),
        genes_after_filter=sum(keep),
        unique_symbols_after_collapse=nrow(out),
        design_columns=ncol(design_primary),
        design_rank=qr(design_primary)$rank,
        output_file=out_file,
        stringsAsFactors=FALSE
    )
}

summary_df <- do.call(rbind, summary_rows)

write.table(
    summary_df,
    file=file.path(
        qcdir,
        "GSE283746_REGULATORY_LOGCPM_PREPROCESSING_QC_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n=== PREPROCESSING SUMMARY ===\n")
print(summary_df)
cat("\nGSE283746 REGULATORY LOGCPM PREPROCESSING COMPLETE.\n")
