#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

expr_files <- c(
    GPL10558 = file.path(
        root,
        "results/DS001_GSE38900/cell_composition/audit",
        "DS001_GSE38900_GPL10558_cell_composition_input_HUGO.tsv.gz"
    ),
    GPL6884 = file.path(
        root,
        "results/DS001_GSE38900/cell_composition/audit",
        "DS001_GSE38900_GPL6884_cell_composition_input_HUGO.tsv.gz"
    )
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

focus_genes <- c(
    "ISG15",
    "IRF7",
    "IFIH1",
    "STAT1",
    "EIF2AK2",
    "OAS3",
    "OASL",
    "TRIM25"
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

meta$primary_eligible <-
    to_bool(meta$eligible_primary_infection_vs_control)

meta$direct_eligible <-
    to_bool(meta$eligible_direct_influenza_vs_rsv)

meta$rsv_eligible <-
    to_bool(meta$eligible_rsv_acute_vs_control)

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

expr_long <- list()

for (platform in names(expr_files)) {

    d <- read.delim(
        gzfile(expr_files[[platform]]),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    d$gene_symbol <- toupper(trimws(d$gene_symbol))

    d <- d[
        d$gene_symbol %in% focus_genes,
        ,
        drop = FALSE
    ]

    sample_cols <- setdiff(
        colnames(d),
        "gene_symbol"
    )

    tmp <- do.call(
        rbind,
        lapply(seq_len(nrow(d)), function(i) {
            data.frame(
                platform = platform,
                gene_symbol = d$gene_symbol[i],
                sample = sample_cols,
                expression = as.numeric(d[i, sample_cols]),
                stringsAsFactors = FALSE
            )
        })
    )

    expr_long[[platform]] <- tmp
}

expr_long <- do.call(rbind, expr_long)
rownames(expr_long) <- NULL

results <- list()

for (contrast_name in names(contrast_specs)) {

    spec <- contrast_specs[[contrast_name]]

    m <- meta[
        meta$platform_id == spec$platform &
        meta[[spec$eligibility]] &
        meta$harmonized_group %in%
            c(spec$group1, spec$group0),
        ,
        drop = FALSE
    ]

    e <- expr_long[
        expr_long$platform == spec$platform,
        ,
        drop = FALSE
    ]

    p <- pcs[
        pcs$platform == spec$platform,
        ,
        drop = FALSE
    ]

    z <- merge(
        e,
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

    rows <- list()

    for (gene in focus_genes) {

        q <- z[
            z$gene_symbol == gene,
            ,
            drop = FALSE
        ]

        if (nrow(q) == 0) {
            next
        }

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
                "expression",
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

            base_formula <-
                expression ~
                group_binary +
                age_months_harmonized +
                sex_standardized

            adjusted_formula <- as.formula(
                paste(
                    "expression ~ group_binary +",
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
                "expression",
                "group_binary",
                spec$pc_names
            )

            q <- q[
                complete.cases(q[, needed]),
                ,
                drop = FALSE
            ]

            base_formula <-
                expression ~ group_binary

            adjusted_formula <- as.formula(
                paste(
                    "expression ~ group_binary +",
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

        adj_fit <- lm(
            adjusted_formula,
            data = q
        )

        base_coef <-
            summary(base_fit)$coefficients[
                "group_binary",
            ]

        adj_coef <-
            summary(adj_fit)$coefficients[
                "group_binary",
            ]

        beta_base <-
            unname(base_coef["Estimate"])

        beta_adjusted <-
            unname(adj_coef["Estimate"])

        attenuation <- if (
            abs(beta_base) > 1e-12
        ) {
            100 * (
                1 -
                abs(beta_adjusted) /
                abs(beta_base)
            )
        } else {
            NA_real_
        }

        rows[[gene]] <- data.frame(
            contrast = contrast_name,
            platform = spec$platform,
            gene_symbol = gene,
            n = nrow(q),
            beta_base = beta_base,
            p_base =
                unname(
                    base_coef["Pr(>|t|)"]
                ),
            beta_adjusted = beta_adjusted,
            p_adjusted =
                unname(
                    adj_coef["Pr(>|t|)"]
                ),
            attenuation_percent = attenuation,
            sign_preserved =
                sign(beta_base) ==
                sign(beta_adjusted),
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
        "DS001_GSE38900_interferon_focus_gene_composition_effects.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n=== INTERFERON GENE ROBUSTNESS ===\n")

print(
    final[, c(
        "contrast",
        "gene_symbol",
        "beta_base",
        "FDR_base",
        "beta_adjusted",
        "FDR_adjusted",
        "attenuation_percent",
        "sign_preserved"
    )],
    row.names = FALSE
)

cat("\nGene-level analysis complete.\n")
