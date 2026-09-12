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

outdir <- "results/rnaseq_validation/GSE155925_DE/age_sensitivity"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

x <- read.delim(
    count_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

rownames(x) <- x[[1]]
x[[1]] <- NULL

m <- read_tsv(meta_file, show_col_types=FALSE) %>%
    filter(primary_eligible == TRUE)

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

run_model <- function(meta, design_formula, label) {

    meta <- as.data.frame(meta)

    rownames(meta) <- meta$count_matrix_label

    xx <- x[, rownames(meta), drop=FALSE]

    # Same frozen gene filter rule
    keep <- rowSums(xx >= 10) >= 5
    xx <- xx[keep, , drop=FALSE]

    meta$primary_group <- relevel(
        factor(meta$primary_group),
        ref="VIRUS_NEGATIVE"
    )

    if ("age_z" %in% all.vars(design_formula)) {
        meta$age_z <- as.numeric(scale(meta$age_months))
    }

    X <- model.matrix(
        design_formula,
        data=meta
    )

    cat("\n====================================\n")
    cat(label, "\n")
    cat("====================================\n")

    cat("Samples:", nrow(meta), "\n")
    print(table(meta$primary_group))

    cat("\nAge summary:\n")
    print(
        tapply(
            meta$age_months,
            meta$primary_group,
            summary
        )
    )

    cat(
        "\nDesign parameters:", ncol(X),
        " rank:", qr(X)$rank,
        " full rank:",
        qr(X)$rank == ncol(X),
        "\n"
    )

    if (qr(X)$rank != ncol(X)) {
        stop(paste(label, "design is rank deficient"))
    }

    dds <- DESeqDataSetFromMatrix(
        countData=round(as.matrix(xx)),
        colData=meta,
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
            paste0(label, "_RSV_vs_negative.tsv")
        )
    )

    cat(
        "FDR<0.05:",
        sum(z$padj < 0.05, na.rm=TRUE),
        "\n"
    )

    cat(
        "FDR<0.05 & |LFC|>=0.5:",
        sum(
            z$padj < 0.05 &
            abs(z$log2FoldChange) >= 0.5,
            na.rm=TRUE
        ),
        "\n"
    )
}

# ============================================================
# A. Age <=30, all sexes
# ============================================================

mA <- m %>%
    filter(age_months <= 30)

run_model(
    mA,
    ~ primary_group,
    "A_age30_allsex"
)

# ============================================================
# B. Age <=30, male only
# ============================================================

mB <- m %>%
    filter(
        age_months <= 30,
        sex == "Male"
    )

run_model(
    mB,
    ~ age_z + primary_group,
    "B_age30_maleonly"
)

cat("\nAge-sensitivity analyses complete.\n")
