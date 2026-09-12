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

outdir <- "results/rnaseq_validation/GSE155925_DE"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------
# Counts
# ------------------------------------------------------------

x <- read.delim(
    count_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

rownames(x) <- x[[1]]
x[[1]] <- NULL

# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------

m <- read_tsv(meta_file, show_col_types=FALSE)

m <- m %>%
    filter(primary_eligible == TRUE)

rownames(m) <- m$count_matrix_label

# Exact ordering
x <- x[, rownames(m), drop=FALSE]

stopifnot(identical(colnames(x), rownames(m)))

# ------------------------------------------------------------
# Frozen filter
# ------------------------------------------------------------

keep <- rowSums(x >= 10) >= 5
x_f <- x[keep, , drop=FALSE]

cat("Genes before filter:", nrow(x), "\n")
cat("Genes after filter :", nrow(x_f), "\n")

# ------------------------------------------------------------
# Covariates
# ------------------------------------------------------------

m$primary_group <- relevel(
    factor(m$primary_group),
    ref="VIRUS_NEGATIVE"
)

m$sex <- factor(m$sex)
m$hospital_batch <- factor(m$hospital_batch)
m$enrollment_year_batch <- factor(m$enrollment_year_batch)

m$age_z <- as.numeric(scale(m$age_months))

# ------------------------------------------------------------
# Primary adjusted DESeq2
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData=round(as.matrix(x_f)),
    colData=m,
    design=~ age_z + sex + hospital_batch +
            enrollment_year_batch + primary_group
)

dds <- DESeq(dds)

cat("\n=== RESULTS NAMES ===\n")
print(resultsNames(dds))

res <- results(
    dds,
    contrast=c(
        "primary_group",
        "RSV_ONLY",
        "VIRUS_NEGATIVE"
    ),
    alpha=0.05
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_df <- res_df %>%
    select(
        gene_id,
        baseMean,
        log2FoldChange,
        lfcSE,
        stat,
        pvalue,
        padj
    ) %>%
    arrange(padj)

write_tsv(
    res_df,
    file.path(
        outdir,
        "GSE155925_RSV_vs_negative_DESeq2_adjusted.tsv"
    )
)

# ------------------------------------------------------------
# Normalized counts
# ------------------------------------------------------------

norm <- counts(dds, normalized=TRUE)

write.table(
    cbind(gene_id=rownames(norm), norm),
    file.path(
        outdir,
        "GSE155925_normalized_counts.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ------------------------------------------------------------
# VST
# ------------------------------------------------------------

vsd <- vst(dds, blind=FALSE)

write.table(
    cbind(
        gene_id=rownames(assay(vsd)),
        assay(vsd)
    ),
    file.path(
        outdir,
        "GSE155925_VST_expression.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n=== DE SUMMARY ===\n")
print(summary(res))

cat("\nFDR <0.05:",
    sum(res$padj < 0.05, na.rm=TRUE),
    "\n")

cat("FDR <0.05 & |log2FC| >=0.5:",
    sum(
        res$padj < 0.05 &
        abs(res$log2FoldChange) >= 0.5,
        na.rm=TRUE
    ),
    "\n")

cat("FDR <0.05 & |log2FC| >=1:",
    sum(
        res$padj < 0.05 &
        abs(res$log2FoldChange) >= 1,
        na.rm=TRUE
    ),
    "\n")

saveRDS(
    dds,
    file.path(outdir, "GSE155925_primary_dds.rds")
)

cat("\nAnalysis complete.\n")
