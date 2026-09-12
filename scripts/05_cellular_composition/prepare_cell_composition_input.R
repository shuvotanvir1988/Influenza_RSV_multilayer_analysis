#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

files <- c(
    GPL10558 = file.path(
        root,
        "data/processed/DS001_GSE38900/analysis_ready/GPL10558_gene_expression.tsv.gz"
    ),
    GPL6884 = file.path(
        root,
        "data/processed/DS001_GSE38900/analysis_ready/GPL6884_gene_expression.tsv.gz"
    )
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/audit"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()
duplicate_rows <- list()

for (platform in names(files)) {

    message("Processing ", platform)

    dat <- read.delim(
        gzfile(files[[platform]]),
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )

    required <- c("gene_symbol", "entrez_gene_id", "probe_id")

    if (!all(required %in% colnames(dat))) {
        stop(platform, ": expected identifier columns are missing")
    }

    sample_cols <- setdiff(colnames(dat), required)

    dat$gene_symbol <- trimws(as.character(dat$gene_symbol))

    valid <- !is.na(dat$gene_symbol) &
             dat$gene_symbol != "" &
             dat$gene_symbol != "NA"

    dat <- dat[valid, , drop = FALSE]

    dup_symbols <- unique(
        dat$gene_symbol[duplicated(dat$gene_symbol)]
    )

    if (length(dup_symbols) > 0) {

        dups <- dat[
            dat$gene_symbol %in% dup_symbols,
            c(required, sample_cols),
            drop = FALSE
        ]

        dups$platform <- platform
        duplicate_rows[[platform]] <- dups
    }

    x <- dat[, sample_cols, drop = FALSE]

    for (nm in sample_cols) {
        x[[nm]] <- as.numeric(x[[nm]])
    }

    aggregate_input <- cbind(
        gene_symbol = dat$gene_symbol,
        x,
        stringsAsFactors = FALSE
    )

    collapsed <- aggregate(
        . ~ gene_symbol,
        data = aggregate_input,
        FUN = median,
        na.rm = TRUE
    )

    rownames(collapsed) <- collapsed$gene_symbol
    collapsed$gene_symbol <- NULL

    collapsed <- as.matrix(collapsed)
    storage.mode(collapsed) <- "numeric"

    if (anyDuplicated(rownames(collapsed))) {
        stop(platform, ": duplicate symbols remain after collapse")
    }

    if (anyNA(collapsed)) {
        stop(platform, ": NA values detected after collapse")
    }

    outfile <- file.path(
        outdir,
        paste0(
            "DS001_GSE38900_",
            platform,
            "_cell_composition_input.tsv.gz"
        )
    )

    con <- gzfile(outfile, "wt")

    write.table(
        cbind(gene_symbol = rownames(collapsed), collapsed),
        file = con,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    close(con)

    summary_rows[[platform]] <- data.frame(
        platform = platform,
        original_rows = nrow(dat),
        unique_gene_symbols = length(unique(dat$gene_symbol)),
        duplicated_gene_symbols = length(dup_symbols),
        final_rows = nrow(collapsed),
        sample_count = ncol(collapsed),
        missing_values = sum(is.na(collapsed)),
        matrix_min = min(collapsed),
        matrix_median = median(collapsed),
        matrix_max = max(collapsed)
    )

    message(
        platform,
        ": ",
        nrow(collapsed),
        " genes x ",
        ncol(collapsed),
        " samples"
    )
}

summary_table <- do.call(rbind, summary_rows)

write.table(
    summary_table,
    file.path(
        outdir,
        "DS001_GSE38900_cell_composition_input_preparation_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

if (length(duplicate_rows) > 0) {

    duplicate_table <- do.call(rbind, duplicate_rows)

    write.table(
        duplicate_table,
        file.path(
            outdir,
            "DS001_GSE38900_duplicate_gene_symbols_before_cell_composition.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

message("Cell-composition input preparation complete.")
