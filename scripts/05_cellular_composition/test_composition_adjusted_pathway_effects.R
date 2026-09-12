#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

domain_file <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/tables",
    "DS001_GSE38900_domain_sample_scores.tsv"
)

meta_file <- file.path(
    root,
    "data/processed/DS001_GSE38900",
    "DS001_GSE38900_master_sample_metadata_harmonized.tsv"
)

pc_file <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/sensitivity",
    "DS001_GSE38900_MCPcounter_composition_PCs.tsv"
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/sensitivity"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

domains <- read.delim(
    domain_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

meta <- read.delim(
    meta_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

pcs <- read.delim(
    pc_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

to_bool <- function(x) {
    tolower(trimws(as.character(x))) %in%
        c("true", "t", "1", "yes")
}

meta$primary_eligible <- to_bool(
    meta$eligible_primary_infection_vs_control
)

meta$direct_eligible <- to_bool(
    meta$eligible_direct_influenza_vs_rsv
)

meta$rsv_eligible <- to_bool(
    meta$eligible_rsv_acute_vs_control
)

contrast_specs <- list(

    GPL10558_RSV_vs_control = list(
        platform = "GPL10558",
        group1 = "RSV acute",
        group0 = "Healthy control",
        eligibility = "rsv_eligible",
        pc_names = paste0("PC", 1:4),
        adjust_age_sex = FALSE
    ),

    GPL6884_RSV_vs_control = list(
        platform = "GPL6884",
        group1 = "RSV acute",
        group0 = "Healthy control",
        eligibility = "rsv_eligible",
        pc_names = paste0("PC", 1:5),
        adjust_age_sex = TRUE
    ),

    GPL6884_Influenza_vs_control = list(
        platform = "GPL6884",
        group1 = "Influenza A acute",
        group0 = "Healthy control",
        eligibility = "primary_eligible",
        pc_names = paste0("PC", 1:5),
        adjust_age_sex = TRUE
    ),

    GPL6884_Influenza_vs_RSV = list(
        platform = "GPL6884",
        group1 = "Influenza A acute",
        group0 = "RSV acute",
        eligibility = "direct_eligible",
        pc_names = paste0("PC", 1:5),
        adjust_age_sex = TRUE
    )
)

results <- list()

for (contrast_name in names(contrast_specs)) {

    spec <- contrast_specs[[contrast_name]]

    cat("\n=== ", contrast_name, " ===\n", sep = "")

    m <- meta[
        meta$platform_id == spec$platform &
        meta[[spec$eligibility]] &
        meta$harmonized_group %in%
            c(spec$group1, spec$group0),
        ,
        drop = FALSE
    ]

    d <- domains[
        domains$platform == spec$platform,
        ,
        drop = FALSE
    ]

    p <- pcs[
        pcs$platform == spec$platform,
        ,
        drop = FALSE
    ]

    z <- merge(
        d,
        m[, c(
            "description",
            "harmonized_group",
            "age_months_harmonized",
            "sex_standardized"
        )],
        by.x = "sample",
        by.y = "description",
        all = FALSE
    )

    z <- merge(
        z,
        p,
        by = c("platform", "sample"),
        all = FALSE
    )

    z$group_binary <- ifelse(
        z$harmonized_group == spec$group1,
        1,
        0
    )

    cat(
        "Samples:",
        length(unique(z$sample)),
        "\n"
    )

    rows <- list()

    for (
        domain_name in
        unique(z$preregistered_domain)
    ) {

        q <- z[
            z$preregistered_domain == domain_name,
            ,
            drop = FALSE
        ]

        if (spec$adjust_age_sex) {

            q$age_months_harmonized <-
                suppressWarnings(
                    as.numeric(
                        q$age_months_harmonized
                    )
                )

            q$sex_standardized <-
                factor(q$sex_standardized)

            needed <- c(
                "domain_score",
                "group_binary",
                "age_months_harmonized",
                "sex_standardized",
                spec$pc_names
            )

            q <- q[
                complete.cases(q[, needed]),
                ,
                drop = FALSE
            ]

            base_formula <- as.formula(
                paste(
                    "domain_score ~ group_binary +",
                    "age_months_harmonized +",
                    "sex_standardized"
                )
            )

            adjusted_formula <- as.formula(
                paste(
                    "domain_score ~ group_binary +",
                    "age_months_harmonized +",
                    "sex_standardized +",
                    paste(
                        spec$pc_names,
                        collapse = " + "
                    )
                )
            )

        } else {

            needed <- c(
                "domain_score",
                "group_binary",
                spec$pc_names
            )

            q <- q[
                complete.cases(q[, needed]),
                ,
                drop = FALSE
            ]

            base_formula <-
                domain_score ~ group_binary

            adjusted_formula <- as.formula(
                paste(
                    "domain_score ~ group_binary +",
                    paste(
                        spec$pc_names,
                        collapse = " + "
                    )
                )
            )
        }

        base_fit <- lm(
            base_formula,
            data = q
        )

        adjusted_fit <- lm(
            adjusted_formula,
            data = q
        )

        base_coef <-
            summary(base_fit)$coefficients[
                "group_binary",
            ]

        adj_coef <-
            summary(adjusted_fit)$coefficients[
                "group_binary",
            ]

        beta_base <-
            base_coef["Estimate"]

        beta_adjusted <-
            adj_coef["Estimate"]

        attenuation <- if (
            is.finite(beta_base) &&
            abs(beta_base) > 1e-12
        ) {
            100 *
            (
                1 -
                abs(beta_adjusted) /
                abs(beta_base)
            )
        } else {
            NA_real_
        }

        sign_preserved <-
            sign(beta_base) ==
            sign(beta_adjusted)

        rows[[domain_name]] <- data.frame(
            contrast = contrast_name,
            platform = spec$platform,
            group1 = spec$group1,
            group0 = spec$group0,
            preregistered_domain =
                domain_name,
            n = nrow(q),
            composition_PCs =
                paste(
                    spec$pc_names,
                    collapse = ";"
                ),

            beta_base =
                unname(beta_base),

            beta_base_se =
                unname(
                    base_coef["Std. Error"]
                ),

            p_base =
                unname(
                    base_coef["Pr(>|t|)"]
                ),

            beta_adjusted =
                unname(beta_adjusted),

            beta_adjusted_se =
                unname(
                    adj_coef["Std. Error"]
                ),

            p_adjusted =
                unname(
                    adj_coef["Pr(>|t|)"]
                ),

            attenuation_percent =
                attenuation,

            sign_preserved =
                sign_preserved,

            stringsAsFactors = FALSE
        )
    }

    r <- do.call(rbind, rows)
    rownames(r) <- NULL

    r$FDR_base <- p.adjust(
        r$p_base,
        method = "BH"
    )

    r$FDR_adjusted <- p.adjust(
        r$p_adjusted,
        method = "BH"
    )

    results[[contrast_name]] <- r
}

final <- do.call(rbind, results)
rownames(final) <- NULL

write.table(
    final,
    file.path(
        outdir,
        "DS001_GSE38900_composition_adjusted_pathway_effects.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n=== COMPOSITION-ADJUSTED PATHWAY EFFECTS ===\n\n")

show <- final[, c(
    "contrast",
    "preregistered_domain",
    "n",
    "beta_base",
    "FDR_base",
    "beta_adjusted",
    "FDR_adjusted",
    "attenuation_percent",
    "sign_preserved"
)]

print(show, row.names = FALSE)

cat("\nAnalysis complete.\n")
