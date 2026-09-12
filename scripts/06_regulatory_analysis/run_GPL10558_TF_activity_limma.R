suppressPackageStartupMessages({
    library(limma)
})

ROOT <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

activity_file <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/replication",
    "GPL10558_COLLECTRI_ULM_ACTIVITY_MIN5_v1.0.tsv.gz"
)

mapping_file <- file.path(
    ROOT,
    "results/DS001_GSE38900/differential_expression/tables",
    "DS001_GSE38900_expression_metadata_sample_mapping.tsv"
)

coverage_file <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/qc",
    "GPL10558_COLLECTRI_TARGET_COVERAGE_v1.0.tsv"
)

out_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/replication/differential_activity"
)

qc_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/qc"
)

design_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/design"
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(design_dir, recursive=TRUE, showWarnings=FALSE)

cat("=== LOAD GPL10558 TF ACTIVITY ===\n")

act <- read.delim(
    activity_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

if (!"expression_sample_id" %in% colnames(act)) {
    stop("Missing expression_sample_id column.")
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

if (anyNA(activity_mat)) {
    stop("Missing values detected in GPL10558 TF activity matrix.")
}

cat("Activity matrix dimensions (TF x sample):\n")
print(dim(activity_mat))

cat("\n=== LOAD FROZEN SAMPLE MAPPING ===\n")

map <- read.delim(
    mapping_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

map10558 <- map[
    map$platform_id == "GPL10558",
    ,
    drop=FALSE
]

if (nrow(map10558) != 36) {
    stop(
        paste(
            "Expected 36 GPL10558 samples; found",
            nrow(map10558)
        )
    )
}

if (!all(colnames(activity_mat) == map10558$expression_sample_id)) {
    stop("GPL10558 activity sample order does not match frozen mapping.")
}

group10558 <- factor(
    map10558$harmonized_group,
    levels=c(
        "Healthy control",
        "RSV acute"
    )
)

design10558 <- model.matrix(
    ~ 0 + group10558
)

colnames(design10558) <- c(
    "Healthy_control",
    "RSV_acute"
)

cat("\nGPL10558 sample counts:\n")
print(table(group10558))

cat("\nGPL10558 design dimensions:\n")
print(dim(design10558))

cat("\nGPL10558 design rank:\n")
print(qr(design10558)$rank)

if (qr(design10558)$rank != ncol(design10558)) {
    stop("GPL10558 regulatory design matrix is not full rank.")
}

write.table(
    design10558,
    file=file.path(
        design_dir,
        "GPL10558_regulatory_primary_design_matrix_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    col.names=NA
)

contrast10558 <- makeContrasts(
    RSV_acute_vs_control =
        RSV_acute - Healthy_control,
    levels=design10558
)

cat("\n=== FIT GPL10558 LIMMA MODEL ===\n")

fit <- lmFit(
    activity_mat,
    design10558
)

fit <- contrasts.fit(
    fit,
    contrast10558
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
    stop("Not all GPL10558 inferred TFs are present in >=5-target coverage table.")
}

tt <- topTable(
    fit,
    coef="RSV_acute_vs_control",
    number=Inf,
    adjust.method="BH",
    sort.by="none"
)

tt$TF <- rownames(tt)

result <- data.frame(
    TF=tt$TF,
    contrast="RSV acute vs healthy control",
    delta_activity=tt$logFC,
    average_activity=tt$AveExpr,
    moderated_t=tt$t,
    P_value=tt$P.Value,
    FDR=tt$adj.P.Val,
    B=tt$B,
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
    "higher_activity_in_RSV",
    ifelse(
        result$delta_activity < 0,
        "lower_activity_in_RSV",
        "no_difference"
    )
)

out_file <- file.path(
    out_dir,
    "GPL10558_RSV_acute_vs_control_TF_activity_limma_v1.0.tsv"
)

write.table(
    result,
    file=out_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

summary_df <- data.frame(
    contrast="RSV acute vs healthy control",
    TFs_tested=nrow(result),
    FDR_lt_0.05=sum(result$FDR < 0.05, na.rm=TRUE),
    FDR_lt_0.05_positive=sum(
        result$FDR < 0.05 & result$delta_activity > 0,
        na.rm=TRUE
    ),
    FDR_lt_0.05_negative=sum(
        result$FDR < 0.05 & result$delta_activity < 0,
        na.rm=TRUE
    ),
    median_abs_delta_activity=median(
        abs(result$delta_activity),
        na.rm=TRUE
    ),
    stringsAsFactors=FALSE
)

write.table(
    summary_df,
    file=file.path(
        out_dir,
        "GPL10558_TF_activity_limma_summary_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

qc <- data.frame(
    platform="GPL10558",
    activity_TFs=nrow(activity_mat),
    activity_samples=ncol(activity_mat),
    design_rows=nrow(design10558),
    design_columns=ncol(design10558),
    design_rank=qr(design10558)$rank,
    contrasts=1,
    missing_activity_values=sum(is.na(activity_mat)),
    stringsAsFactors=FALSE
)

write.table(
    qc,
    file=file.path(
        qc_dir,
        "GPL10558_TF_ACTIVITY_LIMMA_QC_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n=== GPL10558 DIFFERENTIAL TF-ACTIVITY SUMMARY ===\n")
print(summary_df)

cat("\n=== QC ===\n")
print(qc)

cat("\nWritten:\n")
cat(out_file, "\n")
cat("\nGPL10558 DIFFERENTIAL TF-ACTIVITY ANALYSIS COMPLETE.\n")
