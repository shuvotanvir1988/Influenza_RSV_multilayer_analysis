suppressPackageStartupMessages({
    library(limma)
})

ROOT <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

activity_file <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/bulk",
    "GPL6884_COLLECTRI_ULM_ACTIVITY_MIN5_v1.0.tsv.gz"
)

mapping_file <- file.path(
    ROOT,
    "results/DS001_GSE38900/differential_expression/tables",
    "DS001_GSE38900_expression_metadata_sample_mapping.tsv"
)

coverage_file <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/qc",
    "GPL6884_COLLECTRI_TARGET_COVERAGE_v1.0.tsv"
)

out_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/bulk/differential_activity"
)

design_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/design"
)

qc_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/qc"
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(design_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)

cat("=== LOAD TF ACTIVITY ===\n")

act <- read.delim(
    activity_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

if (!"expression_sample_id" %in% colnames(act)) {
    stop("Missing expression_sample_id column in activity matrix.")
}

sample_ids <- act$expression_sample_id
tf_names <- setdiff(colnames(act), "expression_sample_id")

activity_mat <- t(
    as.matrix(
        act[, tf_names, drop=FALSE]
    )
)

mode(activity_mat) <- "numeric"

colnames(activity_mat) <- sample_ids
rownames(activity_mat) <- tf_names

cat("Activity matrix dimensions (TF x sample):\n")
print(dim(activity_mat))

if (anyNA(activity_mat)) {
    stop("Missing values detected in TF activity matrix.")
}

cat("\n=== LOAD FROZEN SAMPLE MAPPING ===\n")

map <- read.delim(
    mapping_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

map6884 <- map[
    map$platform_id == "GPL6884",
    ,
    drop=FALSE
]

if (nrow(map6884) != 205) {
    stop(
        paste(
            "Expected 205 GPL6884 mapping rows; found",
            nrow(map6884)
        )
    )
}

if (!all(colnames(activity_mat) == map6884$expression_sample_id)) {
    stop("TF activity sample order does not match frozen GPL6884 mapping.")
}

map6884$age_months_harmonized <- as.numeric(
    map6884$age_months_harmonized
)

if (any(is.na(map6884$age_months_harmonized))) {
    stop("GPL6884 contains missing harmonized ages.")
}

if (any(is.na(map6884$sex_standardized))) {
    stop("GPL6884 contains missing standardized sex.")
}

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

sex6884 <- factor(
    map6884$sex_standardized,
    levels=c(
        "Female",
        "Male"
    )
)

age6884 <- map6884$age_months_harmonized

design6884 <- model.matrix(
    ~ 0 + group6884 + age6884 + sex6884
)

colnames(design6884) <- c(
    "Healthy_control",
    "RSV_acute",
    "Influenza_A_acute",
    "HRV_acute",
    "RSV_recovery",
    "age_months",
    "sex_Male"
)

cat("\nGPL6884 sample counts:\n")
print(table(group6884))

cat("\nGPL6884 sex counts:\n")
print(table(sex6884))

cat("\nGPL6884 age summary:\n")
print(summary(age6884))

cat("\nGPL6884 design dimensions:\n")
print(dim(design6884))

cat("\nGPL6884 design rank:\n")
print(qr(design6884)$rank)

if (qr(design6884)$rank != ncol(design6884)) {
    stop("GPL6884 regulatory design matrix is not full rank.")
}

write.table(
    design6884,
    file=file.path(
        design_dir,
        "GPL6884_regulatory_primary_design_matrix_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    col.names=NA
)

contrasts6884 <- makeContrasts(
    RSV_acute_vs_control =
        RSV_acute - Healthy_control,

    Influenza_A_vs_control =
        Influenza_A_acute - Healthy_control,

    Influenza_A_vs_RSV_acute =
        Influenza_A_acute - RSV_acute,

    levels=design6884
)

cat("\n=== FIT LIMMA MODEL ===\n")

fit <- lmFit(
    activity_mat,
    design6884
)

fit <- contrasts.fit(
    fit,
    contrasts6884
)

fit <- eBayes(fit)

coverage <- read.delim(
    coverage_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

coverage <- coverage[
    coverage$measurable_targets >= 5,
    ,
    drop=FALSE
]

if (!all(rownames(activity_mat) %in% coverage$source)) {
    stop("Not all inferred TFs are present in frozen GPL6884 coverage table.")
}

save_contrast <- function(
    fit,
    coef_name,
    filename,
    contrast_label
) {
    tt <- topTable(
        fit,
        coef=coef_name,
        number=Inf,
        adjust.method="BH",
        sort.by="none"
    )

    tt$TF <- rownames(tt)

    result <- data.frame(
        TF = tt$TF,
        contrast = contrast_label,
        delta_activity = tt$logFC,
        average_activity = tt$AveExpr,
        moderated_t = tt$t,
        P_value = tt$P.Value,
        FDR = tt$adj.P.Val,
        B = tt$B,
        stringsAsFactors=FALSE
    )

    result <- merge(
        result,
        coverage[
            ,
            c(
                "source",
                "unique_targets_total",
                "measurable_targets",
                "coverage_fraction",
                "eligible_min5",
                "eligible_min10"
            )
        ],
        by.x="TF",
        by.y="source",
        all.x=TRUE,
        sort=FALSE
    )

    result <- result[
        match(rownames(activity_mat), result$TF),
        ,
        drop=FALSE
    ]

    result$direction <- ifelse(
        result$delta_activity > 0,
        "higher_activity_in_numerator",
        ifelse(
            result$delta_activity < 0,
            "lower_activity_in_numerator",
            "no_difference"
        )
    )

    write.table(
        result,
        file=filename,
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    result
}

rsv <- save_contrast(
    fit,
    "RSV_acute_vs_control",
    file.path(
        out_dir,
        "GPL6884_RSV_acute_vs_control_TF_activity_limma_v1.0.tsv"
    ),
    "RSV acute vs healthy control"
)

flu <- save_contrast(
    fit,
    "Influenza_A_vs_control",
    file.path(
        out_dir,
        "GPL6884_InfluenzaA_vs_control_TF_activity_limma_v1.0.tsv"
    ),
    "Influenza A acute vs healthy control"
)

flu_rsv <- save_contrast(
    fit,
    "Influenza_A_vs_RSV_acute",
    file.path(
        out_dir,
        "GPL6884_InfluenzaA_vs_RSVacute_TF_activity_limma_v1.0.tsv"
    ),
    "Influenza A acute vs RSV acute"
)

summarize_contrast <- function(x, label) {
    data.frame(
        contrast=label,
        TFs_tested=nrow(x),
        FDR_lt_0.05=sum(x$FDR < 0.05, na.rm=TRUE),
        FDR_lt_0.05_positive=sum(
            x$FDR < 0.05 & x$delta_activity > 0,
            na.rm=TRUE
        ),
        FDR_lt_0.05_negative=sum(
            x$FDR < 0.05 & x$delta_activity < 0,
            na.rm=TRUE
        ),
        median_abs_delta_activity=median(
            abs(x$delta_activity),
            na.rm=TRUE
        ),
        stringsAsFactors=FALSE
    )
}

summary_df <- rbind(
    summarize_contrast(
        rsv,
        "RSV acute vs healthy control"
    ),
    summarize_contrast(
        flu,
        "Influenza A acute vs healthy control"
    ),
    summarize_contrast(
        flu_rsv,
        "Influenza A acute vs RSV acute"
    )
)

write.table(
    summary_df,
    file=file.path(
        out_dir,
        "GPL6884_TF_activity_limma_summary_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

qc <- data.frame(
    platform="GPL6884",
    activity_TFs=nrow(activity_mat),
    activity_samples=ncol(activity_mat),
    design_rows=nrow(design6884),
    design_columns=ncol(design6884),
    design_rank=qr(design6884)$rank,
    contrasts=3,
    missing_activity_values=sum(is.na(activity_mat)),
    stringsAsFactors=FALSE
)

write.table(
    qc,
    file=file.path(
        qc_dir,
        "GPL6884_TF_ACTIVITY_LIMMA_QC_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n=== DIFFERENTIAL TF-ACTIVITY SUMMARY ===\n")
print(summary_df)

cat("\n=== QC ===\n")
print(qc)

cat("\nOutputs written to:\n")
cat(out_dir, "\n")
cat("\nGPL6884 DIFFERENTIAL TF-ACTIVITY ANALYSIS COMPLETE.\n")
