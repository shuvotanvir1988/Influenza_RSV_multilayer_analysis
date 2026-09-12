#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(immunedeconv)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

indir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/audit"
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/estimates"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

platforms <- c("GPL10558", "GPL6884")

summary_list <- list()

for (platform in platforms) {

    message("\n=== ", platform, " ===")

    infile <- file.path(
        indir,
        paste0(
            "DS001_GSE38900_",
            platform,
            "_cell_composition_input_HUGO.tsv.gz"
        )
    )

    dat <- read.delim(
        gzfile(infile),
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )

    genes <- dat$gene_symbol

    expr <- as.matrix(
        dat[, setdiff(colnames(dat), "gene_symbol"), drop = FALSE]
    )

    storage.mode(expr) <- "numeric"
    rownames(expr) <- genes

    if (anyDuplicated(rownames(expr))) {
        stop(platform, ": duplicated gene symbols remain")
    }

    if (anyNA(expr)) {
        stop(platform, ": NA values detected")
    }

    message(
        "Input: ",
        nrow(expr),
        " genes x ",
        ncol(expr),
        " samples"
    )

    res <- immunedeconv::deconvolute_mcp_counter(
        expr,
        feature_types = "HUGO_symbols"
    )

    # immunedeconv 2.1.0 stores cell populations as row names.
    if (is.null(rownames(res))) {
        stop(platform, ": MCP-counter returned no row names")
    }

    res_out <- data.frame(
        cell_type = rownames(res),
        res,
        check.names = FALSE,
        row.names = NULL
    )

    outfile <- file.path(
        outdir,
        paste0(
            "DS001_GSE38900_",
            platform,
            "_MCPcounter_scores.tsv"
        )
    )

    write.table(
        res_out,
        outfile,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    numeric_part <- res_out[
        ,
        setdiff(colnames(res_out), "cell_type"),
        drop = FALSE
    ]

    summary_list[[platform]] <- data.frame(
        platform = platform,
        input_genes = nrow(expr),
        samples = ncol(expr),
        estimated_cell_populations = nrow(res_out),
        missing_estimates = sum(is.na(numeric_part)),
        negative_estimates = sum(numeric_part < 0, na.rm = TRUE)
    )

    message(
        "Populations: ",
        paste(res_out$cell_type, collapse = "; ")
    )

    message("Estimated populations: ", nrow(res_out))
    message("Missing estimates: ", sum(is.na(numeric_part)))
}

summary_table <- do.call(rbind, summary_list)

write.table(
    summary_table,
    file.path(
        outdir,
        "DS001_GSE38900_MCPcounter_run_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

message("\nMCP-counter estimation complete.")
