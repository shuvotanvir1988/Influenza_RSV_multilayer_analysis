#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

meta_file <- file.path(
    root,
    "data/processed/DS001_GSE38900/DS001_GSE38900_master_sample_metadata_harmonized.tsv"
)

estimate_dir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/estimates"
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/tables"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(
    meta_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

to_bool <- function(x) {
    tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

meta$primary_eligible <-
    to_bool(meta$eligible_primary_infection_vs_control)

meta$direct_eligible <-
    to_bool(meta$eligible_direct_influenza_vs_rsv)

meta$rsv_eligible <-
    to_bool(meta$eligible_rsv_acute_vs_control)

primary_cells <- c(
    "T cells",
    "CD8 T cells",
    "Cytotoxic lymphocytes",
    "NK cells",
    "B lineage",
    "Monocytic lineage",
    "Myeloid dendritic cells",
    "Neutrophils"
)

hedges_g <- function(x1, x0) {

    x1 <- x1[is.finite(x1)]
    x0 <- x0[is.finite(x0)]

    n1 <- length(x1)
    n0 <- length(x0)

    s1 <- sd(x1)
    s0 <- sd(x0)

    pooled_sd <- sqrt(
        ((n1 - 1) * s1^2 + (n0 - 1) * s0^2) /
        (n1 + n0 - 2)
    )

    if (!is.finite(pooled_sd) || pooled_sd == 0) {
        return(c(
            hedges_g = NA,
            ci_low = NA,
            ci_high = NA
        ))
    }

    d <- (mean(x1) - mean(x0)) / pooled_sd

    correction <- 1 - 3 / (4 * (n1 + n0) - 9)

    g <- correction * d

    se_g <- sqrt(
        (n1 + n0) / (n1 * n0) +
        g^2 / (2 * (n1 + n0 - 2))
    )

    c(
        hedges_g = g,
        ci_low = g - 1.96 * se_g,
        ci_high = g + 1.96 * se_g
    )
}

load_scores_long <- function(platform) {

    f <- file.path(
        estimate_dir,
        paste0(
            "DS001_GSE38900_",
            platform,
            "_MCPcounter_scores.tsv"
        )
    )

    if (!file.exists(f)) {
        stop("Missing MCP-counter file: ", f)
    }

    x <- read.delim(
        f,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    if (!"cell_type" %in% colnames(x)) {
        stop(platform, ": cell_type column missing")
    }

    samples <- setdiff(colnames(x), "cell_type")

    out <- do.call(
        rbind,
        lapply(seq_len(nrow(x)), function(i) {
            data.frame(
                cell_type = x$cell_type[i],
                description = samples,
                score = as.numeric(x[i, samples]),
                stringsAsFactors = FALSE
            )
        })
    )

    rownames(out) <- NULL
    out
}

all_results <- list()
group_summaries <- list()

run_contrast <- function(
    platform,
    group1,
    group0,
    contrast_name,
    eligibility_column,
    adjust_age_sex = FALSE
) {

    cat("\n=== ", contrast_name, " ===\n", sep = "")

    scores <- load_scores_long(platform)

    m <- meta[
        meta$platform_id == platform &
        meta[[eligibility_column]],
        ,
        drop = FALSE
    ]

    m <- m[
        m$harmonized_group %in% c(group1, group0),
        ,
        drop = FALSE
    ]

    expected_samples <- unique(m$description)
    score_samples <- unique(scores$description)

    missing_from_scores <- setdiff(
        expected_samples,
        score_samples
    )

    if (length(missing_from_scores) > 0) {
        stop(
            contrast_name,
            ": metadata samples missing from MCP-counter: ",
            paste(missing_from_scores, collapse = "; ")
        )
    }

    d <- merge(
        scores,
        m,
        by = "description",
        all = FALSE
    )

    cat("Metadata samples:", nrow(m), "\n")
    cat("Matched samples:", length(unique(d$description)), "\n")

    print(table(m$harmonized_group))

    d$group_binary <- ifelse(
        d$harmonized_group == group1,
        1,
        0
    )

    results <- list()
    summaries <- list()

    for (ct in unique(d$cell_type)) {

        z <- d[d$cell_type == ct, , drop = FALSE]

        x1 <- z$score[
            z$harmonized_group == group1
        ]

        x0 <- z$score[
            z$harmonized_group == group0
        ]

        effect <- hedges_g(x1, x0)

        summaries[[ct]] <- data.frame(
            platform = platform,
            contrast = contrast_name,
            cell_type = ct,
            primary_immune_population = ct %in% primary_cells,
            group1 = group1,
            group0 = group0,
            n_group1 = length(x1),
            n_group0 = length(x0),
            median_group1 = median(x1),
            median_group0 = median(x0),
            median_difference = median(x1) - median(x0),
            mean_group1 = mean(x1),
            mean_group0 = mean(x0),
            stringsAsFactors = FALSE
        )

        if (adjust_age_sex) {

            z$age_months_harmonized <-
                suppressWarnings(
                    as.numeric(z$age_months_harmonized)
                )

            z$sex_standardized <-
                factor(z$sex_standardized)

            complete <- complete.cases(
                z[, c(
                    "score",
                    "group_binary",
                    "age_months_harmonized",
                    "sex_standardized"
                )]
            )

            zm <- z[complete, , drop = FALSE]

            if (nrow(zm) < 4) {
                stop(
                    contrast_name,
                    " / ",
                    ct,
                    ": insufficient complete cases"
                )
            }

            fit <- lm(
                score ~ group_binary +
                    age_months_harmonized +
                    sex_standardized,
                data = zm
            )

            coef_table <- summary(fit)$coefficients

            if (!"group_binary" %in% rownames(coef_table)) {
                stop(
                    contrast_name,
                    " / ",
                    ct,
                    ": group coefficient unavailable"
                )
            }

            beta <- coef_table[
                "group_binary",
                "Estimate"
            ]

            se <- coef_table[
                "group_binary",
                "Std. Error"
            ]

            p <- coef_table[
                "group_binary",
                "Pr(>|t|)"
            ]

            model_n <- nrow(zm)

            model_type <-
                "age_sex_adjusted_linear_model"

        } else {

            fit <- lm(
                score ~ group_binary,
                data = z
            )

            coef_table <- summary(fit)$coefficients

            beta <- coef_table[
                "group_binary",
                "Estimate"
            ]

            se <- coef_table[
                "group_binary",
                "Std. Error"
            ]

            p <- coef_table[
                "group_binary",
                "Pr(>|t|)"
            ]

            model_n <- nrow(z)

            model_type <- "unadjusted_linear_model"
        }

        wilcox_p <- wilcox.test(
            x1,
            x0,
            exact = FALSE
        )$p.value

        results[[ct]] <- data.frame(
            platform = platform,
            contrast = contrast_name,
            cell_type = ct,
            primary_immune_population = ct %in% primary_cells,
            group1 = group1,
            group0 = group0,
            n_group1 = length(x1),
            n_group0 = length(x0),
            hedges_g = unname(effect["hedges_g"]),
            hedges_g_ci_low = unname(effect["ci_low"]),
            hedges_g_ci_high = unname(effect["ci_high"]),
            model_type = model_type,
            model_n = model_n,
            beta_group1_vs_group0 = beta,
            beta_se = se,
            beta_ci_low = beta - 1.96 * se,
            beta_ci_high = beta + 1.96 * se,
            model_p = p,
            wilcoxon_p = wilcox_p,
            stringsAsFactors = FALSE
        )
    }

    results <- do.call(rbind, results)
    summaries <- do.call(rbind, summaries)

    results$model_FDR <- p.adjust(
        results$model_p,
        method = "BH"
    )

    results$wilcoxon_FDR <- p.adjust(
        results$wilcoxon_p,
        method = "BH"
    )

    all_results[[contrast_name]] <<- results
    group_summaries[[contrast_name]] <<- summaries
}

run_contrast(
    platform = "GPL10558",
    group1 = "RSV acute",
    group0 = "Healthy control",
    contrast_name = "GPL10558_RSV_vs_control",
    eligibility_column = "rsv_eligible",
    adjust_age_sex = FALSE
)

run_contrast(
    platform = "GPL6884",
    group1 = "RSV acute",
    group0 = "Healthy control",
    contrast_name = "GPL6884_RSV_vs_control",
    eligibility_column = "rsv_eligible",
    adjust_age_sex = TRUE
)

run_contrast(
    platform = "GPL6884",
    group1 = "Influenza A acute",
    group0 = "Healthy control",
    contrast_name = "GPL6884_Influenza_vs_control",
    eligibility_column = "primary_eligible",
    adjust_age_sex = TRUE
)

run_contrast(
    platform = "GPL6884",
    group1 = "Influenza A acute",
    group0 = "RSV acute",
    contrast_name = "GPL6884_Influenza_vs_RSV",
    eligibility_column = "direct_eligible",
    adjust_age_sex = TRUE
)

result_table <- do.call(rbind, all_results)
summary_table <- do.call(rbind, group_summaries)

rownames(result_table) <- NULL
rownames(summary_table) <- NULL

write.table(
    result_table,
    file.path(
        outdir,
        "DS001_GSE38900_MCPcounter_disease_comparisons.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    summary_table,
    file.path(
        outdir,
        "DS001_GSE38900_MCPcounter_group_summaries.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n=== PRIMARY IMMUNE RESULTS ===\n")

primary <- result_table[
    result_table$primary_immune_population,
]

print(
    primary[, c(
        "contrast",
        "cell_type",
        "n_group1",
        "n_group0",
        "hedges_g",
        "hedges_g_ci_low",
        "hedges_g_ci_high",
        "beta_group1_vs_group0",
        "model_p",
        "model_FDR"
    )],
    row.names = FALSE
)

cat("\nDisease-group MCP-counter analysis complete.\n")
