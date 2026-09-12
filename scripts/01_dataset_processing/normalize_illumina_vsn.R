#!/usr/bin/env Rscript

# Script: normalize_illumina_vsn.R
# Project: Influenza A and RSV Systems Immunology
# Dataset: DS001_GSE38900
# Session: 14, revised normalization stage
# Version: 1.0
#
# Purpose:
# Apply variance-stabilizing normalization independently to GPL10558
# and GPL6884 supplementary Illumina signal matrices.
#
# Inputs:
# data/interim/DS001_GSE38900/illumina_raw/
#   DS001_GSE38900_GPL10558_signal_matrix.tsv.gz
#   DS001_GSE38900_GPL6884_signal_matrix.tsv.gz
#
# Outputs:
# data/processed/DS001_GSE38900/preprocessing/
#   DS001_GSE38900_GPL10558_expression_vsn.tsv.gz
#   DS001_GSE38900_GPL6884_expression_vsn.tsv.gz
#
# results/DS001_GSE38900/preprocessing/tables/
#   DS001_GSE38900_vsn_normalization_qc.tsv
#
# logs/session14/
#   DS001_GSE38900_vsn_normalization_summary.txt
#   DS001_GSE38900_GPL10558_vsn_fit.rds
#   DS001_GSE38900_GPL6884_vsn_fit.rds

suppressPackageStartupMessages({
    library(vsn)
})

options(
    stringsAsFactors = FALSE,
    warn = 1
)

root <- file.path(
    Sys.getenv("HOME"),
    "Influenza_RSV_Project"
)

dataset <- "DS001_GSE38900"

input_dir <- file.path(
    root,
    "data",
    "interim",
    dataset,
    "illumina_raw"
)

output_dir <- file.path(
    root,
    "data",
    "processed",
    dataset,
    "preprocessing"
)

results_dir <- file.path(
    root,
    "results",
    dataset,
    "preprocessing",
    "tables"
)

log_dir <- file.path(
    root,
    "logs",
    "session14"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    results_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    log_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

platforms <- c(
    GPL10558 = 36L,
    GPL6884 = 205L
)

expected_probes <- c(
    GPL10558 = 47323L,
    GPL6884 = 48803L
)

safe_quantile <- function(x, probability) {
    as.numeric(
        quantile(
            x,
            probs = probability,
            na.rm = TRUE,
            names = FALSE,
            type = 7
        )
    )
}

qc_rows <- list()
started_at <- Sys.time()

cat(
    "VSN normalization started:",
    format(started_at, tz = "UTC"),
    "\n"
)

for (platform in names(platforms)) {

    cat("\n===", platform, "===\n")

    input_path <- file.path(
        input_dir,
        paste0(
            dataset,
            "_",
            platform,
            "_signal_matrix.tsv.gz"
        )
    )

    if (!file.exists(input_path)) {
        stop(
            paste(
                platform,
                "input file not found:",
                input_path
            )
        )
    }

    input_table <- read.delim(
        gzfile(input_path),
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        quote = "",
        comment.char = "",
        na.strings = c("", "NA")
    )

    input_matrix <- as.matrix(input_table)
    storage.mode(input_matrix) <- "double"

    if (nrow(input_matrix) != expected_probes[[platform]]) {
        stop(
            paste(
                platform,
                "expected",
                expected_probes[[platform]],
                "probes but found",
                nrow(input_matrix)
            )
        )
    }

    if (ncol(input_matrix) != platforms[[platform]]) {
        stop(
            paste(
                platform,
                "expected",
                platforms[[platform]],
                "samples but found",
                ncol(input_matrix)
            )
        )
    }

    if (anyDuplicated(rownames(input_matrix)) > 0) {
        stop(
            paste(
                platform,
                "contains duplicate probe IDs."
            )
        )
    }

    if (anyDuplicated(colnames(input_matrix)) > 0) {
        stop(
            paste(
                platform,
                "contains duplicate sample labels."
            )
        )
    }

    if (anyNA(input_matrix)) {
        stop(
            paste(
                platform,
                "contains missing signal values."
            )
        )
    }

    if (any(!is.finite(input_matrix))) {
        stop(
            paste(
                platform,
                "contains non-finite signal values."
            )
        )
    }

    input_values <- as.numeric(input_matrix)

    cat(
        "Input dimensions:",
        nrow(input_matrix),
        "x",
        ncol(input_matrix),
        "\n"
    )

    cat("Fitting VSN model...\n")

    fit <- vsn2(
        input_matrix,
        verbose = TRUE
    )

    normalized_matrix <- predict(
        fit,
        input_matrix
    )

    if (!all(dim(normalized_matrix) == dim(input_matrix))) {
        stop(
            paste(
                platform,
                "normalized dimensions do not match input."
            )
        )
    }

    if (anyNA(normalized_matrix)) {
        stop(
            paste(
                platform,
                "VSN output contains missing values."
            )
        )
    }

    if (any(!is.finite(normalized_matrix))) {
        stop(
            paste(
                platform,
                "VSN output contains non-finite values."
            )
        )
    }

    rownames(normalized_matrix) <- rownames(input_matrix)
    colnames(normalized_matrix) <- colnames(input_matrix)

    output_path <- file.path(
        output_dir,
        paste0(
            dataset,
            "_",
            platform,
            "_expression_vsn.tsv.gz"
        )
    )

    connection <- gzfile(
        output_path,
        open = "wt"
    )

    write.table(
        normalized_matrix,
        file = connection,
        sep = "\t",
        quote = FALSE,
        col.names = NA,
        row.names = TRUE,
        na = ""
    )

    close(connection)

    fit_path <- file.path(
        log_dir,
        paste0(
            dataset,
            "_",
            platform,
            "_vsn_fit.rds"
        )
    )

    saveRDS(
        fit,
        fit_path
    )

    normalized_values <- as.numeric(normalized_matrix)

    before_sample_medians <- apply(
        input_matrix,
        2,
        median
    )

    after_sample_medians <- apply(
        normalized_matrix,
        2,
        median
    )

    qc_rows[[platform]] <- data.frame(
        dataset = dataset,
        platform = platform,
        probe_count = nrow(input_matrix),
        sample_count = ncol(input_matrix),

        input_min = min(input_values),
        input_q01 = safe_quantile(input_values, 0.01),
        input_median = median(input_values),
        input_q99 = safe_quantile(input_values, 0.99),
        input_max = max(input_values),
        input_negative_values = sum(input_values < 0),

        vsn_min = min(normalized_values),
        vsn_q01 = safe_quantile(normalized_values, 0.01),
        vsn_median = median(normalized_values),
        vsn_q99 = safe_quantile(normalized_values, 0.99),
        vsn_max = max(normalized_values),

        input_sample_median_range =
            max(before_sample_medians) -
            min(before_sample_medians),

        vsn_sample_median_range =
            max(after_sample_medians) -
            min(after_sample_medians),

        missing_output_values =
            sum(is.na(normalized_matrix)),

        output_file = output_path,
        fit_file = fit_path,
        status = "PASS",
        check.names = FALSE
    )

    cat("Output:", output_path, "\n")
    cat(
        "VSN value range:",
        min(normalized_values),
        "to",
        max(normalized_values),
        "\n"
    )
    cat("Status: PASS\n")
}

qc_table <- do.call(
    rbind,
    qc_rows
)

qc_path <- file.path(
    results_dir,
    paste0(
        dataset,
        "_vsn_normalization_qc.tsv"
    )
)

write.table(
    qc_table,
    file = qc_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
)

finished_at <- Sys.time()

summary_path <- file.path(
    log_dir,
    paste0(
        dataset,
        "_vsn_normalization_summary.txt"
    )
)

summary_lines <- c(
    paste(
        "Dataset:",
        dataset
    ),
    paste(
        "Started:",
        format(started_at, tz = "UTC")
    ),
    paste(
        "Finished:",
        format(finished_at, tz = "UTC")
    ),
    paste(
        "Elapsed seconds:",
        round(
            as.numeric(
                difftime(
                    finished_at,
                    started_at,
                    units = "secs"
                )
            ),
            3
        )
    ),
    paste(
        "VSN package version:",
        as.character(
            packageVersion("vsn")
        )
    ),
    paste(
        "Platforms:",
        paste(names(platforms), collapse = ", ")
    ),
    paste(
        "Overall status:",
        if (
            all(qc_table$status == "PASS")
        ) "PASS" else "FAIL"
    ),
    paste(
        "QC table:",
        qc_path
    )
)

writeLines(
    summary_lines,
    con = summary_path
)

cat("\n=== OVERALL ===\n")
cat("Status: PASS\n")
cat("QC table:", qc_path, "\n")
cat("Summary:", summary_path, "\n")
