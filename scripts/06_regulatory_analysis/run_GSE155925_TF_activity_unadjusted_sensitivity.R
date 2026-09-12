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
    "results/regulatory_driver_analysis/rnaseq_validation/diagnostics"
)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

act <- read.delim(activity_file, check.names=FALSE, stringsAsFactors=FALSE)
m <- read.delim(meta_file, check.names=FALSE, stringsAsFactors=FALSE)
drivers <- read.delim(driver_file, check.names=FALSE, stringsAsFactors=FALSE)

m <- m[tolower(as.character(m$primary_eligible)) == "true", , drop=FALSE]

sample_ids <- act$count_matrix_label
tf_names <- setdiff(colnames(act), "count_matrix_label")

A <- t(as.matrix(act[, tf_names, drop=FALSE]))
mode(A) <- "numeric"
colnames(A) <- sample_ids
rownames(A) <- tf_names

m <- m[match(colnames(A), m$count_matrix_label), , drop=FALSE]
if (!all(colnames(A) == m$count_matrix_label)) stop("Alignment failure.")

m$primary_group <- relevel(factor(m$primary_group), ref="VIRUS_NEGATIVE")

# Unadjusted sensitivity model
design0 <- model.matrix(~ primary_group, data=m)
fit0 <- eBayes(lmFit(A, design0))
coef0 <- grep("^primary_groupRSV_ONLY$", colnames(design0), value=TRUE)

tt0 <- topTable(fit0, coef=coef0, number=Inf, adjust.method="BH", sort.by="none")
tt0$TF <- rownames(tt0)

res0 <- data.frame(
    TF=tt0$TF,
    unadjusted_delta_activity=tt0$logFC,
    unadjusted_P=tt0$P.Value,
    unadjusted_FDR=tt0$adj.P.Val,
    stringsAsFactors=FALSE
)

val0 <- merge(
    drivers[,c("TF","rsv_delta","RSV_replication_tier","frozen170_targets")],
    res0, by="TF", all.x=TRUE, sort=FALSE
)

val0 <- val0[match(drivers$TF,val0$TF),,drop=FALSE]
val0$direction_concordant <- sign(val0$rsv_delta) == sign(val0$unadjusted_delta_activity)

summary0 <- data.frame(
    frozen_drivers=nrow(val0),
    direction_concordant=sum(val0$direction_concordant,na.rm=TRUE),
    nominal_support=sum(val0$unadjusted_P < 0.05,na.rm=TRUE),
    FDR_support=sum(val0$unadjusted_FDR < 0.05,na.rm=TRUE),
    stringsAsFactors=FALSE
)

write.table(
    res0,
    file=file.path(out_dir,"GSE155925_all_TF_activity_unadjusted_v1.0.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

write.table(
    val0,
    file=file.path(out_dir,"GSE155925_frozen37_unadjusted_sensitivity_v1.0.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

write.table(
    summary0,
    file=file.path(out_dir,"GSE155925_frozen37_unadjusted_summary_v1.0.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
)

cat("=== UNADJUSTED SENSITIVITY SUMMARY ===\n")
print(summary0)
