suppressPackageStartupMessages({
    library(limma)
    library(readr)
    library(dplyr)
})

base <- "data/rnaseq_validation/influenza_SIG"
metabase <- "results/rnaseq_validation/influenza_SIG/metadata"
outbase <- "results/rnaseq_validation/influenza_SIG/DE"

dir.create(outbase, recursive=TRUE, showWarnings=FALSE)

info <- list(
    GSE158592=list(year=2018),
    GSE155635=list(year=2019),
    GSE196350=list(year=2020),
    GSE213168=list(year=2022)
)

for (gse in names(info)) {

    year <- info[[gse]]$year

    cat("\n========================================\n")
    cat(gse, year, "\n")
    cat("========================================\n")

    expr_file <- file.path(
        base, gse, "processed",
        paste0(gse, "_PARSED_normalized_expression.tsv.gz")
    )

    meta_file <- file.path(
        metabase,
        paste0(gse, "_expression_metadata_map.tsv")
    )

    d <- read_tsv(expr_file, show_col_types=FALSE)
    m <- read_tsv(meta_file, show_col_types=FALSE)

    sample_cols <- m$expression_label

    stopifnot(all(sample_cols %in% colnames(d)))

    # Collapse any duplicated gene symbols by maximum variance
    x <- d %>%
        select(gene_symbol, all_of(sample_cols))

    x$gene_symbol <- as.character(x$gene_symbol)

    mat <- as.matrix(x[, sample_cols])
    mode(mat) <- "numeric"

    variances <- apply(mat, 1, var, na.rm=TRUE)

    tmp <- data.frame(
        gene_symbol=x$gene_symbol,
        row_index=seq_len(nrow(x)),
        variance=variances
    )

    keep_idx <- tmp %>%
        arrange(gene_symbol, desc(variance)) %>%
        group_by(gene_symbol) %>%
        slice(1) %>%
        ungroup() %>%
        pull(row_index)

    mat <- mat[keep_idx, , drop=FALSE]
    rownames(mat) <- x$gene_symbol[keep_idx]

    # Metadata order
    m <- m[match(colnames(mat), m$expression_label), ]

    m$status <- factor(
        m$status,
        levels=c("HEALTHY_CONTROL", "INFLUENZA_INFECTED")
    )

    # --------------------------------------------
    # Primary unadjusted model
    # --------------------------------------------

    design0 <- model.matrix(~ status, data=m)

    fit0 <- lmFit(mat, design0)
    fit0 <- eBayes(fit0)

    coef0 <- grep(
        "statusINFLUENZA_INFECTED",
        colnames(design0),
        value=TRUE
    )

    tab0 <- topTable(
        fit0,
        coef=coef0,
        number=Inf,
        sort.by="P"
    )

    tab0$gene_symbol <- rownames(tab0)

    write_tsv(
        tab0,
        file.path(
            outbase,
            paste0(
                gse,
                "_influenza_vs_control_limma_unadjusted.tsv"
            )
        )
    )

    cat(
        "Primary FDR<0.05:",
        sum(tab0$adj.P.Val < 0.05, na.rm=TRUE),
        "\n"
    )

    # --------------------------------------------
    # Sex-adjusted sensitivity
    # --------------------------------------------

    ms <- m %>%
        filter(sex %in% c("Female", "Male"))

    mats <- mat[, ms$expression_label, drop=FALSE]

    ms$sex <- factor(ms$sex)
    ms$status <- factor(
        ms$status,
        levels=c("HEALTHY_CONTROL", "INFLUENZA_INFECTED")
    )

    design1 <- model.matrix(~ sex + status, data=ms)

    cat(
        "Sex model samples:",
        nrow(ms),
        " parameters:",
        ncol(design1),
        " rank:",
        qr(design1)$rank,
        "\n"
    )

    if (qr(design1)$rank == ncol(design1)) {

        fit1 <- lmFit(mats, design1)
        fit1 <- eBayes(fit1)

        coef1 <- grep(
            "statusINFLUENZA_INFECTED",
            colnames(design1),
            value=TRUE
        )

        tab1 <- topTable(
            fit1,
            coef=coef1,
            number=Inf,
            sort.by="P"
        )

        tab1$gene_symbol <- rownames(tab1)

        write_tsv(
            tab1,
            file.path(
                outbase,
                paste0(
                    gse,
                    "_influenza_vs_control_limma_sex_adjusted.tsv"
                )
            )
        )

        cat(
            "Sex-adjusted FDR<0.05:",
            sum(tab1$adj.P.Val < 0.05, na.rm=TRUE),
            "\n"
        )
    }
}

cat("\nAll seasonal limma analyses complete.\n")
