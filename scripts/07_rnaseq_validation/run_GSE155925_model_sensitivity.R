suppressPackageStartupMessages({
    library(DESeq2)
    library(readr)
    library(dplyr)
})

count_file <- paste0(
    "data/rnaseq_validation/GSE155925/raw_counts/",
    "GSE155925_Raw_counts_matrix.txt.gz"
)

meta_file <- paste0(
    "results/rnaseq_validation/tables/",
    "GSE155925_sample_eligibility_FROZEN_v1.0.tsv"
)

outdir <- "results/rnaseq_validation/GSE155925_DE/model_sensitivity"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

x <- read.delim(
    count_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

rownames(x) <- x[[1]]
x[[1]] <- NULL

m <- read_tsv(meta_file, show_col_types=FALSE)
m <- m %>% filter(primary_eligible == TRUE)

rownames(m) <- m$count_matrix_label
x <- x[, rownames(m), drop=FALSE]

keep <- rowSums(x >= 10) >= 5
x <- x[keep, , drop=FALSE]

m$primary_group <- relevel(
    factor(m$primary_group),
    ref="VIRUS_NEGATIVE"
)

m$sex <- factor(m$sex)
m$hospital_batch <- factor(m$hospital_batch)
m$enrollment_year_batch <- factor(m$enrollment_year_batch)
m$age_z <- as.numeric(scale(m$age_months))

models <- list(
    M0_unadjusted =
        ~ primary_group,

    M1_age =
        ~ age_z + primary_group,

    M2_batch =
        ~ hospital_batch +
          enrollment_year_batch +
          primary_group,

    M3_full =
        ~ age_z + sex +
          hospital_batch +
          enrollment_year_batch +
          primary_group
)

summary_list <- list()

for (nm in names(models)) {

    cat("\n====================================\n")
    cat(nm, "\n")
    cat("====================================\n")

    design_formula <- models[[nm]]

    X <- model.matrix(design_formula, data=m)

    cat("Parameters:", ncol(X), "\n")
    cat("Rank:", qr(X)$rank, "\n")
    cat("Full rank:", qr(X)$rank == ncol(X), "\n")
    cat("Condition number:", kappa(X), "\n")

    if (qr(X)$rank != ncol(X)) {
        cat("SKIPPING rank-deficient model\n")
        next
    }

    dds <- DESeqDataSetFromMatrix(
        countData=round(as.matrix(x)),
        colData=m,
        design=design_formula
    )

    dds <- DESeq(dds, quiet=TRUE)

    res <- results(
        dds,
        contrast=c(
            "primary_group",
            "RSV_ONLY",
            "VIRUS_NEGATIVE"
        ),
        alpha=0.05
    )

    z <- as.data.frame(res)
    z$gene_id <- rownames(z)

    z <- z %>%
        select(
            gene_id,
            baseMean,
            log2FoldChange,
            lfcSE,
            stat,
            pvalue,
            padj
        )

    write_tsv(
        z,
        file.path(
            outdir,
            paste0("GSE155925_", nm, "_RSV_vs_negative.tsv")
        )
    )

    summary_list[[nm]] <- data.frame(
        model=nm,
        FDR05=sum(z$padj < 0.05, na.rm=TRUE),
        FDR05_absLFC05=sum(
            z$padj < 0.05 &
            abs(z$log2FoldChange) >= 0.5,
            na.rm=TRUE
        ),
        median_abs_LFC=median(
            abs(z$log2FoldChange),
            na.rm=TRUE
        )
    )
}

summary_df <- bind_rows(summary_list)

write_tsv(
    summary_df,
    file.path(outdir, "GSE155925_model_sensitivity_summary.tsv")
)

cat("\n=== MODEL SUMMARY ===\n")
print(summary_df)
