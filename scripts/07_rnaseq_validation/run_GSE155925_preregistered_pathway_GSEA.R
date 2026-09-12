suppressPackageStartupMessages({
    library(fgsea)
    library(readr)
    library(dplyr)
    library(tibble)
})

de_file <- paste0(
    "results/rnaseq_validation/GSE155925_DE/",
    "GSE155925_RSV_vs_negative_DESeq2_adjusted.tsv"
)

membership_file <- paste0(
    "results/rnaseq_validation/validation/pathway_validation/",
    "GSE155925_preregistered_pathway_membership_FROZEN_v1.0.tsv.gz"
)

outdir <- "results/rnaseq_validation/validation/pathway_validation/"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. DESeq2 ranked statistics
# ============================================================

de <- read_tsv(de_file, show_col_types=FALSE)

cat("=== DE TABLE ===\n")
cat("Rows:", nrow(de), "\n")

if (!all(c("gene_id", "stat") %in% colnames(de))) {
    stop("Required DESeq2 columns gene_id/stat are missing.")
}

de <- de %>%
    mutate(
        gene_symbol = sub(":.*$", "", gene_id)
    )

rank_df <- de %>%
    filter(
        !is.na(stat),
        !is.na(gene_symbol),
        gene_symbol != ""
    ) %>%
    group_by(gene_symbol) %>%
    slice_max(
        order_by=abs(stat),
        n=1,
        with_ties=FALSE
    ) %>%
    ungroup()

ranks <- rank_df$stat
names(ranks) <- rank_df$gene_symbol
ranks <- sort(ranks, decreasing=TRUE)

cat("\n=== RANKED LIST ===\n")
cat("Genes:", length(ranks), "\n")
cat("Positive:", sum(ranks > 0), "\n")
cat("Negative:", sum(ranks < 0), "\n")
cat("Zero:", sum(ranks == 0), "\n")

# ============================================================
# 2. Frozen preregistered pathways
# ============================================================

m <- read_tsv(
    membership_file,
    show_col_types=FALSE
)

pathways <- split(
    m$gene_symbol,
    m$pathway
)

pathways <- lapply(
    pathways,
    function(x) unique(x[!is.na(x) & x != ""])
)

cat("\n=== FROZEN PATHWAYS ===\n")
cat("Unique pathways:", length(pathways), "\n")

# ============================================================
# 3. RNA-seq gene-set coverage
# ============================================================

coverage <- tibble(
    pathway=names(pathways),
    frozen_size=sapply(pathways, length),
    RNAseq_size=sapply(
        pathways,
        function(x) sum(x %in% names(ranks))
    )
) %>%
    mutate(
        coverage_fraction=RNAseq_size/frozen_size
    )

cat("\n=== RNA-SEQ PATHWAY COVERAGE ===\n")
print(summary(coverage$coverage_fraction))

cat(
    "Pathways with >=15 evaluable genes:",
    sum(coverage$RNAseq_size >= 15),
    "\n"
)

# ============================================================
# 4. Preranked GSEA
# ============================================================

set.seed(240810)

fg <- fgseaMultilevel(
    pathways=pathways,
    stats=ranks,
    minSize=15,
    maxSize=500,
    eps=0
)

fg <- as_tibble(fg) %>%
    arrange(padj)

fg_out <- fg %>%
    mutate(
        leadingEdge=sapply(
            leadingEdge,
            function(x) paste(x, collapse=";")
        )
    )

# ============================================================
# 5. Restore frozen microarray annotations
# ============================================================

annotation <- m %>%
    select(
        pathway,
        preregistered_domain,
        biological_domain,
        NES_RSV_GPL10558,
        padj_RSV_GPL10558,
        NES_RSV_GPL6884,
        padj_RSV_GPL6884,
        RSV_replicated
    ) %>%
    distinct()

final <- annotation %>%
    left_join(coverage, by="pathway") %>%
    left_join(fg_out, by="pathway") %>%
    mutate(
        RNAseq_significant =
            !is.na(padj) & padj < 0.05,

        concordant_GPL10558 =
            !is.na(NES) &
            !is.na(NES_RSV_GPL10558) &
            sign(NES) == sign(NES_RSV_GPL10558),

        concordant_GPL6884 =
            !is.na(NES) &
            !is.na(NES_RSV_GPL6884) &
            sign(NES) == sign(NES_RSV_GPL6884),

        concordant_both_microarrays =
            concordant_GPL10558 &
            concordant_GPL6884
    )

# ============================================================
# 6. Save complete results
# ============================================================

write_tsv(
    final,
    paste0(
        outdir,
        "GSE155925_preregistered_pathway_GSEA_M3.tsv"
    )
)

write_tsv(
    coverage,
    paste0(
        outdir,
        "GSE155925_preregistered_pathway_RNAseq_coverage.tsv"
    )
)

# ============================================================
# 7. Summary
# ============================================================

cat("\n=== GSEA RESULTS ===\n")
cat("Annotation rows:", nrow(final), "\n")
cat("Unique pathways:", n_distinct(final$pathway), "\n")
cat("GSEA-evaluable pathways:", sum(!is.na(final$NES)), "\n")

cat(
    "RNA-seq FDR <0.05:",
    sum(final$RNAseq_significant, na.rm=TRUE),
    "\n"
)

cat(
    "Direction concordant with GPL10558:",
    sum(final$concordant_GPL10558, na.rm=TRUE),
    "/",
    sum(!is.na(final$NES)),
    "\n"
)

cat(
    "Direction concordant with GPL6884:",
    sum(final$concordant_GPL6884, na.rm=TRUE),
    "/",
    sum(!is.na(final$NES)),
    "\n"
)

cat(
    "Concordant with BOTH microarrays:",
    sum(final$concordant_both_microarrays, na.rm=TRUE),
    "/",
    sum(!is.na(final$NES)),
    "\n"
)

cat("\n=== SIGNIFICANT PATHWAYS BY DOMAIN ===\n")

print(
    final %>%
        filter(RNAseq_significant) %>%
        count(preregistered_domain, sort=TRUE)
)

cat("\n=== TOP 30 RNA-SEQ PATHWAYS ===\n")

print(
    final %>%
        arrange(padj) %>%
        select(
            pathway,
            preregistered_domain,
            NES,
            padj,
            NES_RSV_GPL10558,
            NES_RSV_GPL6884,
            RSV_replicated,
            concordant_both_microarrays
        ) %>%
        head(30),
    n=30
)

cat("\nAnalysis complete.\n")
