suppressPackageStartupMessages({
    library(limma)
})

ROOT <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

activity_file <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/rnaseq_validation",
    "GSE155925_COLLECTRI_ULM_ACTIVITY_MIN5_v1.0.tsv.gz"
)

meta_file <- file.path(
    ROOT,
    "results/rnaseq_validation/tables",
    "GSE155925_sample_eligibility_FROZEN_v1.0.tsv"
)

driver_file <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/tables",
    "REGULATORY_HIGH_CONFIDENCE_DRIVERS_v1.0.tsv"
)

out_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/rnaseq_validation"
)

qc_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/qc"
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)

act <- read.delim(
    activity_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

m <- read.delim(
    meta_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

drivers <- read.delim(
    driver_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

m <- m[
    tolower(as.character(m$primary_eligible)) == "true",
    ,
    drop=FALSE
]

sample_ids <- act$count_matrix_label
tf_names <- setdiff(colnames(act), "count_matrix_label")

A <- t(as.matrix(act[, tf_names, drop=FALSE]))
mode(A) <- "numeric"
colnames(A) <- sample_ids
rownames(A) <- tf_names

# Reorder metadata exactly to activity columns
m <- m[match(colnames(A), m$count_matrix_label), , drop=FALSE]

if (any(is.na(m$count_matrix_label))) {
    stop("Metadata alignment failure.")
}

if (!all(colnames(A) == m$count_matrix_label)) {
    stop("Activity/metadata sample order mismatch.")
}

m$primary_group <- relevel(
    factor(m$primary_group),
    ref="VIRUS_NEGATIVE"
)

m$sex <- factor(m$sex)
m$hospital_batch <- factor(m$hospital_batch)
m$enrollment_year_batch <- factor(m$enrollment_year_batch)
m$age_z <- as.numeric(scale(as.numeric(m$age_months)))

design <- model.matrix(
    ~ age_z + sex + hospital_batch +
      enrollment_year_batch + primary_group,
    data=m
)

cat("=== GSE155925 REGULATORY DESIGN ===\n")
print(dim(design))
cat("Rank:", qr(design)$rank, "\n")
print(table(m$primary_group))

if (qr(design)$rank != ncol(design)) {
    stop("GSE155925 regulatory design is not full rank.")
}

fit <- lmFit(A, design)
fit <- eBayes(fit)

coef_name <- grep(
    "^primary_groupRSV_ONLY$",
    colnames(design),
    value=TRUE
)

if (length(coef_name) != 1) {
    stop("Could not uniquely identify RSV_ONLY coefficient.")
}

tt <- topTable(
    fit,
    coef=coef_name,
    number=Inf,
    adjust.method="BH",
    sort.by="none"
)

tt$TF <- rownames(tt)

res <- data.frame(
    TF=tt$TF,
    rnaseq_delta_activity=tt$logFC,
    rnaseq_average_activity=tt$AveExpr,
    rnaseq_moderated_t=tt$t,
    rnaseq_P=tt$P.Value,
    rnaseq_FDR=tt$adj.P.Val,
    stringsAsFactors=FALSE
)

val <- merge(
    drivers,
    res,
    by="TF",
    all.x=TRUE,
    sort=FALSE
)

val <- val[match(drivers$TF, val$TF), , drop=FALSE]

val$evaluable_rnaseq <- !is.na(val$rnaseq_delta_activity)
val$direction_concordant <- (
    sign(val$rsv_delta) == sign(val$rnaseq_delta_activity)
)
val$nominal_support <- val$rnaseq_P < 0.05
val$FDR_support <- val$rnaseq_FDR < 0.05

write.table(
    res,
    file=file.path(
        out_dir,
        "GSE155925_all_TF_activity_limma_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

write.table(
    val,
    file=file.path(
        out_dir,
        "GSE155925_frozen37_driver_validation_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

tier1 <- grepl(
    "^Tier1_",
    val$RSV_replication_tier
)

summary <- data.frame(
    frozen_drivers=37,
    evaluable=sum(val$evaluable_rnaseq, na.rm=TRUE),
    direction_concordant=sum(val$evaluable_rnaseq & val$direction_concordant, na.rm=TRUE),
    nominal_support=sum(val$evaluable_rnaseq & val$nominal_support, na.rm=TRUE),
    FDR_support=sum(val$evaluable_rnaseq & val$FDR_support, na.rm=TRUE),
    tier1_drivers=sum(tier1, na.rm=TRUE),
    tier1_evaluable=sum(tier1 & val$evaluable_rnaseq, na.rm=TRUE),
    tier1_direction_concordant=sum(tier1 & val$evaluable_rnaseq & val$direction_concordant, na.rm=TRUE),
    tier1_nominal_support=sum(tier1 & val$evaluable_rnaseq & val$nominal_support, na.rm=TRUE),
    tier1_FDR_support=sum(tier1 & val$evaluable_rnaseq & val$FDR_support, na.rm=TRUE),
    stringsAsFactors=FALSE
)

write.table(
    summary,
    file=file.path(
        out_dir,
        "GSE155925_frozen37_validation_summary_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

qc <- data.frame(
    samples=ncol(A),
    TFs=nrow(A),
    design_rows=nrow(design),
    design_columns=ncol(design),
    design_rank=qr(design)$rank,
    missing_activity_values=sum(is.na(A)),
    stringsAsFactors=FALSE
)

write.table(
    qc,
    file=file.path(
        qc_dir,
        "GSE155925_TF_ACTIVITY_LIMMA_QC_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n=== FROZEN 37-DRIVER VALIDATION SUMMARY ===\n")
print(summary)

cat("\n=== TOP VALIDATED DRIVERS ===\n")
show <- val[
    order(val$rnaseq_FDR, val$rnaseq_P),
    c(
        "TF",
        "regulatory_class",
        "rsv_delta",
        "rnaseq_delta_activity",
        "rnaseq_P",
        "rnaseq_FDR",
        "direction_concordant",
        "RSV_replication_tier",
        "frozen170_targets"
    )
]
print(head(show, 37))

cat("\nGSE155925 REGULATORY VALIDATION COMPLETE.\n")
