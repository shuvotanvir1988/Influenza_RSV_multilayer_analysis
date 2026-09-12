suppressPackageStartupMessages({
    library(limma)
    library(readr)
    library(dplyr)
    library(stringr)
})

ROOT <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

meta_file <- file.path(
    ROOT,
    "results/scrnaseq_validation/GSE283746/pseudobulk",
    "GSE283746_primary_pseudobulk_sample_metadata_FROZEN_v1.0.tsv"
)

driver_file <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/tables",
    "REGULATORY_HIGH_CONFIDENCE_DRIVERS_v1.0.tsv"
)

ulm_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/scrnaseq_validation/ULM"
)

out_dir <- file.path(
    ROOT,
    "results/regulatory_driver_analysis/scrnaseq_validation/differential_activity"
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

meta_all <- read_tsv(meta_file, show_col_types=FALSE)
drivers <- read.delim(driver_file, check.names=FALSE, stringsAsFactors=FALSE)

celltypes <- c(
    "Naive CD4+ T",
    "Naive CD8+ T",
    "TrB",
    "B",
    "NK",
    "Memory CD4+ T",
    "CD14+ Monocyte",
    "Cytotoxic T",
    "Treg",
    "Memory B",
    "CD16+ Monocyte"
)

safe_name <- function(x) {
    x |>
        str_replace_all("\\+", "plus") |>
        str_replace_all("/", "_") |>
        str_replace_all(" ", "_")
}

all_results <- list()
summary_rows <- list()

for (ct in celltypes) {

    safe <- safe_name(ct)

    activity_file <- file.path(
        ulm_dir,
        paste0(
            "GSE283746_",
            safe,
            "_COLLECTRI_ULM_ACTIVITY_MIN5_v1.0.tsv.gz"
        )
    )

    act <- read.delim(
        activity_file,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )

    sample_ids <- act$Sample
    tf_names <- setdiff(colnames(act), "Sample")

    A <- t(as.matrix(act[, tf_names, drop=FALSE]))
    mode(A) <- "numeric"
    colnames(A) <- sample_ids
    rownames(A) <- tf_names

    meta <- meta_all |>
        filter(GroupedAnnotation == ct) |>
        distinct(Sample, .keep_all=TRUE)

    meta <- as.data.frame(meta)
    meta <- meta[match(colnames(A), meta$Sample), , drop=FALSE]

    if (!all(colnames(A) == meta$Sample)) {
        stop(paste(ct, "activity/metadata alignment failure"))
    }

    meta$RSV_status <- factor(
        meta$`Combined Condition`,
        levels=c("Healthy", "RSV")
    )

    meta$sex <- factor(
        meta$Sex,
        levels=c("F", "M")
    )

    meta$age_z <- as.numeric(scale(meta$Age))

    design <- model.matrix(
        ~ age_z + sex + RSV_status,
        data=meta
    )

    if (qr(design)$rank != ncol(design)) {
        stop(paste(ct, "regulatory limma design not full rank"))
    }

    fit <- eBayes(lmFit(A, design))

    coef_name <- grep(
        "^RSV_statusRSV$",
        colnames(design),
        value=TRUE
    )

    if (length(coef_name) != 1) {
        stop(paste(ct, "RSV coefficient not uniquely identified"))
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
        cell_type=ct,
        delta_activity=tt$logFC,
        moderated_t=tt$t,
        P_value=tt$P.Value,
        FDR=tt$adj.P.Val,
        stringsAsFactors=FALSE
    )

    val <- merge(
        drivers[
            ,
            c(
                "TF",
                "rsv_delta",
                "RSV_replication_tier",
                "frozen170_targets"
            )
        ],
        res,
        by="TF",
        all.x=TRUE,
        sort=FALSE
    )

    val <- val[match(drivers$TF, val$TF), , drop=FALSE]

    val$evaluable <- !is.na(val$delta_activity)
    val$direction_concordant <- (
        sign(val$rsv_delta) == sign(val$delta_activity)
    )
    val$nominal_support <- val$P_value < 0.05
    val$FDR_support <- val$FDR < 0.05

    all_results[[length(all_results)+1]] <- val

    summary_rows[[length(summary_rows)+1]] <- data.frame(
        cell_type=ct,
        frozen_drivers=nrow(val),
        evaluable=sum(val$evaluable, na.rm=TRUE),
        direction_concordant=sum(
            val$evaluable & val$direction_concordant,
            na.rm=TRUE
        ),
        nominal_support=sum(
            val$evaluable & val$nominal_support,
            na.rm=TRUE
        ),
        FDR_support=sum(
            val$evaluable & val$FDR_support,
            na.rm=TRUE
        ),
        stringsAsFactors=FALSE
    )
}

master <- do.call(rbind, all_results)
summary_df <- do.call(rbind, summary_rows)

write.table(
    master,
    file=file.path(
        out_dir,
        "GSE283746_frozen37_celltype_validation_master_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

write.table(
    summary_df,
    file=file.path(
        out_dir,
        "GSE283746_frozen37_celltype_validation_summary_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("=== CELL-TYPE VALIDATION SUMMARY ===\n")
print(summary_df)

cat("\n=== TOP CELL TYPES BY DIRECTION CONCORDANCE ===\n")
print(
    summary_df[
        order(
            -summary_df$direction_concordant,
            -summary_df$nominal_support
        ),
    ]
)

cat("\nGSE283746 CELL-TYPE REGULATORY VALIDATION COMPLETE.\n")
