#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(limma)
})

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

mapping_file <- file.path(
    root,
    "results/DS001_GSE38900/differential_expression/tables",
    "DS001_GSE38900_expression_metadata_sample_mapping.tsv"
)

expr10558_file <- file.path(
    root,
    "data/processed/DS001_GSE38900/analysis_ready",
    "GPL10558_gene_expression.tsv.gz"
)

expr6884_file <- file.path(
    root,
    "data/processed/DS001_GSE38900/analysis_ready",
    "GPL6884_gene_expression.tsv.gz"
)

primary_dir <- file.path(
    root,
    "results/DS001_GSE38900/differential_expression/tables"
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/differential_expression/sensitivity"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

mapping <- read.delim(
    mapping_file,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

annotation_cols <- c(
    "gene_symbol",
    "entrez_gene_id",
    "probe_id"
)

save_results <- function(fit, coef_name, annotation, outfile) {

    tt <- topTable(
        fit,
        coef=coef_name,
        number=Inf,
        adjust.method="BH",
        sort.by="P"
    )

    tt$gene_symbol <- rownames(tt)

    tt <- merge(
        tt,
        annotation,
        by="gene_symbol",
        all.x=TRUE,
        sort=FALSE
    )

    write.table(
        tt,
        outfile,
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    tt
}

cat("============================================================\n")
cat("SESSION 16 — DE SENSITIVITY ANALYSES\n")
cat("============================================================\n")

# ============================================================
# 1. GPL10558 SEX-ADJUSTED RSV VS CONTROL
# ============================================================

cat("\n============================================================\n")
cat("1. GPL10558 SEX-ADJUSTED RSV VS CONTROL\n")
cat("============================================================\n")

x10558 <- read.delim(
    gzfile(expr10558_file),
    stringsAsFactors=FALSE,
    check.names=FALSE
)

ann10558 <- x10558[, annotation_cols]

rownames(x10558) <- x10558$gene_symbol

x10558 <- x10558[, !colnames(x10558) %in% annotation_cols]

map10558 <- mapping[
    mapping$platform_id == "GPL10558",
]

map10558 <- map10558[
    match(
        colnames(x10558),
        map10558$expression_sample_id
    ),
]

stopifnot(
    all(
        colnames(x10558) ==
        map10558$expression_sample_id
    )
)

keep10558 <- map10558$harmonized_group %in% c(
    "Healthy control",
    "RSV acute"
)

expr10558 <- as.matrix(
    x10558[, keep10558]
)

mode(expr10558) <- "numeric"

meta10558 <- map10558[keep10558,]

group10558 <- factor(
    meta10558$harmonized_group,
    levels=c(
        "Healthy control",
        "RSV acute"
    )
)

sex10558 <- factor(
    meta10558$sex_standardized,
    levels=c(
        "Female",
        "Male"
    )
)

design10558 <- model.matrix(
    ~ 0 + group10558 + sex10558
)

colnames(design10558) <- c(
    "Healthy_control",
    "RSV_acute",
    "sex_Male"
)

cat("Design dimensions:", dim(design10558), "\n")
cat("Design rank:", qr(design10558)$rank, "\n")

if (
    qr(design10558)$rank != ncol(design10558)
) {
    stop("GPL10558 sensitivity design not full rank.")
}

contrast10558 <- makeContrasts(
    RSV_vs_control =
        RSV_acute - Healthy_control,
    levels=design10558
)

fit10558 <- lmFit(
    expr10558,
    design10558
)

fit10558 <- contrasts.fit(
    fit10558,
    contrast10558
)

fit10558 <- eBayes(fit10558)

sens10558 <- save_results(
    fit10558,
    "RSV_vs_control",
    ann10558,
    file.path(
        outdir,
        "GPL10558_sex_adjusted_RSV_vs_control.tsv"
    )
)

primary10558 <- read.delim(
    file.path(
        primary_dir,
        "GPL10558_RSV_acute_vs_control_full_DE.tsv"
    ),
    stringsAsFactors=FALSE
)

cmp10558 <- merge(
    primary10558[
        , c(
            "gene_symbol",
            "logFC",
            "adj.P.Val"
        )
    ],
    sens10558[
        , c(
            "gene_symbol",
            "logFC",
            "adj.P.Val"
        )
    ],
    by="gene_symbol",
    suffixes=c(
        "_primary",
        "_sex_adjusted"
    )
)

r10558 <- cor(
    cmp10558$logFC_primary,
    cmp10558$logFC_sex_adjusted,
    method="pearson"
)

rho10558 <- cor(
    cmp10558$logFC_primary,
    cmp10558$logFC_sex_adjusted,
    method="spearman"
)

cat("Pearson logFC correlation:", r10558, "\n")
cat("Spearman logFC correlation:", rho10558, "\n")

# ============================================================
# 2. GPL6884 UNADJUSTED VS AGE+SEX-ADJUSTED
# ============================================================

cat("\n============================================================\n")
cat("2. GPL6884 UNADJUSTED VS ADJUSTED MODELS\n")
cat("============================================================\n")

x6884 <- read.delim(
    gzfile(expr6884_file),
    stringsAsFactors=FALSE,
    check.names=FALSE
)

ann6884 <- x6884[, annotation_cols]

rownames(x6884) <- x6884$gene_symbol

x6884 <- x6884[, !colnames(x6884) %in% annotation_cols]

map6884 <- mapping[
    mapping$platform_id == "GPL6884",
]

map6884 <- map6884[
    match(
        colnames(x6884),
        map6884$expression_sample_id
    ),
]

stopifnot(
    all(
        colnames(x6884) ==
        map6884$expression_sample_id
    )
)

expr6884 <- as.matrix(x6884)
mode(expr6884) <- "numeric"

group6884 <- factor(
    map6884$harmonized_group,
    levels=c(
        "Healthy control",
        "RSV acute",
        "Influenza A acute",
        "HRV acute",
        "RSV recovery"
    )
)

design6884_unadj <- model.matrix(
    ~ 0 + group6884
)

colnames(design6884_unadj) <- c(
    "Healthy_control",
    "RSV_acute",
    "Influenza_A_acute",
    "HRV_acute",
    "RSV_recovery"
)

contrast6884 <- makeContrasts(

    RSV_vs_control =
        RSV_acute - Healthy_control,

    Influenza_vs_control =
        Influenza_A_acute - Healthy_control,

    Influenza_vs_RSV =
        Influenza_A_acute - RSV_acute,

    levels=design6884_unadj
)

fit6884_unadj <- lmFit(
    expr6884,
    design6884_unadj
)

fit6884_unadj <- contrasts.fit(
    fit6884_unadj,
    contrast6884
)

fit6884_unadj <- eBayes(
    fit6884_unadj
)

unadj_rsv <- save_results(
    fit6884_unadj,
    "RSV_vs_control",
    ann6884,
    file.path(
        outdir,
        "GPL6884_unadjusted_RSV_vs_control.tsv"
    )
)

unadj_flu <- save_results(
    fit6884_unadj,
    "Influenza_vs_control",
    ann6884,
    file.path(
        outdir,
        "GPL6884_unadjusted_Influenza_vs_control.tsv"
    )
)

unadj_flu_rsv <- save_results(
    fit6884_unadj,
    "Influenza_vs_RSV",
    ann6884,
    file.path(
        outdir,
        "GPL6884_unadjusted_Influenza_vs_RSV.tsv"
    )
)

primary_files <- list(

    RSV_vs_control =
        "GPL6884_RSV_acute_vs_control_full_DE.tsv",

    Influenza_vs_control =
        "GPL6884_InfluenzaA_vs_control_full_DE.tsv",

    Influenza_vs_RSV =
        "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv"
)

unadjusted_tables <- list(
    RSV_vs_control=unadj_rsv,
    Influenza_vs_control=unadj_flu,
    Influenza_vs_RSV=unadj_flu_rsv
)

comparison_rows <- list()

for (name in names(primary_files)) {

    primary <- read.delim(
        file.path(
            primary_dir,
            primary_files[[name]]
        ),
        stringsAsFactors=FALSE
    )

    unadjusted <- unadjusted_tables[[name]]

    cmp <- merge(
        primary[
            , c(
                "gene_symbol",
                "logFC",
                "adj.P.Val"
            )
        ],
        unadjusted[
            , c(
                "gene_symbol",
                "logFC",
                "adj.P.Val"
            )
        ],
        by="gene_symbol",
        suffixes=c(
            "_adjusted",
            "_unadjusted"
        )
    )

    pearson <- cor(
        cmp$logFC_adjusted,
        cmp$logFC_unadjusted,
        method="pearson"
    )

    spearman <- cor(
        cmp$logFC_adjusted,
        cmp$logFC_unadjusted,
        method="spearman"
    )

    same_direction <- mean(
        sign(cmp$logFC_adjusted) ==
        sign(cmp$logFC_unadjusted)
    )

    comparison_rows[[length(comparison_rows)+1]] <- data.frame(
        comparison=name,
        pearson_logFC=pearson,
        spearman_logFC=spearman,
        direction_concordance=same_direction,
        adjusted_FDR_0_05=sum(
            cmp$adj.P.Val_adjusted < 0.05
        ),
        unadjusted_FDR_0_05=sum(
            cmp$adj.P.Val_unadjusted < 0.05
        )
    )

    cat(
        name,
        "\nPearson:", pearson,
        "\nSpearman:", spearman,
        "\nDirection concordance:", same_direction,
        "\n\n"
    )
}

adjustment_summary <- do.call(
    rbind,
    comparison_rows
)

write.table(
    adjustment_summary,
    file.path(
        outdir,
        "GPL6884_unadjusted_vs_adjusted_summary.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# 3. AGE-RESTRICTED INFLUENZA VS RSV
# ============================================================

cat("\n============================================================\n")
cat("3. AGE-RESTRICTED INFLUENZA VS RSV\n")
cat("============================================================\n")

age <- as.numeric(
    map6884$age_months_harmonized
)

# Restrict to age range represented in both disease groups.
#
# Influenza observed range:
# 0.97–18.83 months
#
# RSV observed range:
# 0.37–22.33 months
#
# Therefore use intersection:
# 0.97–18.83 months
#
keep_age <- (
    map6884$harmonized_group %in% c(
        "Influenza A acute",
        "RSV acute"
    )
    &
    age >= 0.97
    &
    age <= 18.83
)

expr_age <- expr6884[, keep_age]

meta_age <- map6884[keep_age,]

age_age <- as.numeric(
    meta_age$age_months_harmonized
)

group_age <- factor(
    meta_age$harmonized_group,
    levels=c(
        "RSV acute",
        "Influenza A acute"
    )
)

sex_age <- factor(
    meta_age$sex_standardized,
    levels=c(
        "Female",
        "Male"
    )
)

cat("Restricted sample counts:\n")
print(table(group_age))

cat("Restricted age summary by group:\n")
print(
    tapply(
        age_age,
        group_age,
        summary
    )
)

design_age <- model.matrix(
    ~ 0 + group_age + age_age + sex_age
)

colnames(design_age) <- c(
    "RSV_acute",
    "Influenza_A_acute",
    "age_months",
    "sex_Male"
)

cat("Design dimensions:", dim(design_age), "\n")
cat("Design rank:", qr(design_age)$rank, "\n")

if (
    qr(design_age)$rank != ncol(design_age)
) {
    stop("Age-restricted model is not full rank.")
}

contrast_age <- makeContrasts(
    Influenza_vs_RSV =
        Influenza_A_acute - RSV_acute,
    levels=design_age
)

fit_age <- lmFit(
    expr_age,
    design_age
)

fit_age <- contrasts.fit(
    fit_age,
    contrast_age
)

fit_age <- eBayes(
    fit_age
)

age_results <- save_results(
    fit_age,
    "Influenza_vs_RSV",
    ann6884,
    file.path(
        outdir,
        "GPL6884_age_restricted_Influenza_vs_RSV.tsv"
    )
)

primary_flu_rsv <- read.delim(
    file.path(
        primary_dir,
        "GPL6884_InfluenzaA_vs_RSVacute_full_DE.tsv"
    ),
    stringsAsFactors=FALSE
)

cmp_age <- merge(
    primary_flu_rsv[
        , c(
            "gene_symbol",
            "logFC",
            "adj.P.Val"
        )
    ],
    age_results[
        , c(
            "gene_symbol",
            "logFC",
            "adj.P.Val"
        )
    ],
    by="gene_symbol",
    suffixes=c(
        "_primary",
        "_age_restricted"
    )
)

pearson_age <- cor(
    cmp_age$logFC_primary,
    cmp_age$logFC_age_restricted,
    method="pearson"
)

spearman_age <- cor(
    cmp_age$logFC_primary,
    cmp_age$logFC_age_restricted,
    method="spearman"
)

direction_age <- mean(
    sign(cmp_age$logFC_primary) ==
    sign(cmp_age$logFC_age_restricted)
)

cat("Primary vs age-restricted Pearson:", pearson_age, "\n")
cat("Primary vs age-restricted Spearman:", spearman_age, "\n")
cat("Direction concordance:", direction_age, "\n")

# ============================================================
# MASTER SENSITIVITY SUMMARY
# ============================================================

master <- data.frame(

    analysis=c(
        "GPL10558 primary vs sex-adjusted",
        "GPL6884 RSV adjusted vs unadjusted",
        "GPL6884 Influenza adjusted vs unadjusted",
        "GPL6884 Influenza-vs-RSV adjusted vs unadjusted",
        "GPL6884 Influenza-vs-RSV primary vs age-restricted"
    ),

    pearson_logFC=c(
        r10558,
        adjustment_summary$pearson_logFC[
            adjustment_summary$comparison ==
            "RSV_vs_control"
        ],
        adjustment_summary$pearson_logFC[
            adjustment_summary$comparison ==
            "Influenza_vs_control"
        ],
        adjustment_summary$pearson_logFC[
            adjustment_summary$comparison ==
            "Influenza_vs_RSV"
        ],
        pearson_age
    ),

    spearman_logFC=c(
        rho10558,
        adjustment_summary$spearman_logFC[
            adjustment_summary$comparison ==
            "RSV_vs_control"
        ],
        adjustment_summary$spearman_logFC[
            adjustment_summary$comparison ==
            "Influenza_vs_control"
        ],
        adjustment_summary$spearman_logFC[
            adjustment_summary$comparison ==
            "Influenza_vs_RSV"
        ],
        spearman_age
    )
)

write.table(
    master,
    file.path(
        outdir,
        "DS001_DE_sensitivity_summary.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n============================================================\n")
cat("MASTER SENSITIVITY SUMMARY\n")
cat("============================================================\n")

print(master)

cat("\nSensitivity analyses completed successfully.\n")
