#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

files <- c(
    GPL10558 = file.path(
        root,
        "results/DS001_GSE38900/cell_composition/estimates",
        "DS001_GSE38900_GPL10558_MCPcounter_scores.tsv"
    ),
    GPL6884 = file.path(
        root,
        "results/DS001_GSE38900/cell_composition/estimates",
        "DS001_GSE38900_GPL6884_MCPcounter_scores.tsv"
    )
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/audit"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

all_summary <- list()

for (platform in names(files)) {

    d <- read.delim(
        files[[platform]],
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )

    if (!"cell_type" %in% colnames(d)) {
        stop(platform, ": cell_type column missing")
    }

    sample_cols <- setdiff(colnames(d), "cell_type")

    for (s in sample_cols) {
        d[[s]] <- as.numeric(d[[s]])
    }

    rows <- lapply(seq_len(nrow(d)), function(i) {

        x <- as.numeric(d[i, sample_cols])

        data.frame(
            platform = platform,
            cell_type = d$cell_type[i],
            samples = length(x),
            missing = sum(is.na(x)),
            minimum = min(x, na.rm = TRUE),
            q25 = unname(quantile(x, 0.25, na.rm = TRUE)),
            median = median(x, na.rm = TRUE),
            mean = mean(x, na.rm = TRUE),
            q75 = unname(quantile(x, 0.75, na.rm = TRUE)),
            maximum = max(x, na.rm = TRUE),
            sd = sd(x, na.rm = TRUE),
            zero_variance = sd(x, na.rm = TRUE) == 0
        )
    })

    all_summary[[platform]] <- do.call(rbind, rows)
}

summary_table <- do.call(rbind, all_summary)

write.table(
    summary_table,
    file.path(
        outdir,
        "DS001_GSE38900_MCPcounter_score_QC.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

print(summary_table)
