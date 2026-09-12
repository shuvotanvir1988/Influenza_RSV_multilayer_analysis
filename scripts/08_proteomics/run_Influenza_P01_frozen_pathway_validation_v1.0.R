suppressPackageStartupMessages({
    library(data.table)
    library(fgsea)
})

ROOT <- "results/proteomics_validation/Influenza-P01"

PRIMARY <- file.path(
    ROOT,
    "differential_proteomics",
    "Influenza-P01_primary_adjusted_infection_limma_v1.0.tsv"
)

ANNOT <- file.path(
    ROOT,
    "analysis_inputs",
    "Influenza-P01_human_assay_annotation_v1.0.tsv"
)

HYP <- paste0(
    "results/proteomics_validation/frozen_targets/",
    "PROTEOMICS_PATHWAY_HYPOTHESES_FROZEN_v1.0.tsv"
)

REACTOME <- paste0(
    "results/DS001_GSE38900/pathway_analysis/gene_sets/",
    "MSigDB_Reactome_Homo_sapiens.tsv.gz"
)

OUTDIR <- file.path(
    ROOT,
    "pathway_validation"
)

dir.create(
    OUTDIR,
    recursive=TRUE,
    showWarnings=FALSE
)

# ============================================================
# 1. READ INPUTS
# ============================================================

res <- fread(PRIMARY)

ann <- fread(ANNOT)

hyp58 <- fread(HYP)

gs <- as.data.table(
    read.delim(
        gzfile(REACTOME),
        sep="\t",
        header=TRUE,
        stringsAsFactors=FALSE,
        check.names=FALSE
    )
)

cat("=== INPUT DIMENSIONS ===\n")
cat("Differential assays:", nrow(res), "\n")
cat("Annotation rows:", nrow(ann), "\n")
cat("Frozen pathway rows:", nrow(hyp58), "\n")
cat(
    "Unique frozen pathways:",
    uniqueN(hyp58$pathway),
    "\n"
)
cat(
    "Reactome membership rows:",
    nrow(gs),
    "\n\n"
)

if (nrow(hyp58) != 58) {
    stop(
        paste(
            "Expected 58 frozen pathway rows;",
            "found",
            nrow(hyp58)
        )
    )
}

if (uniqueN(hyp58$pathway) != 56) {
    stop(
        paste(
            "Expected 56 unique frozen pathways;",
            "found",
            uniqueN(hyp58$pathway)
        )
    )
}

# ============================================================
# 2. AUDIT DUPLICATE DOMAIN ASSIGNMENTS
# ============================================================

duplicate_audit <- hyp58[
    duplicated(pathway) |
    duplicated(pathway, fromLast=TRUE)
][
    order(pathway, preregistered_domain)
]

cat("=== DUPLICATE PATHWAY DOMAIN ASSIGNMENTS ===\n")

print(
    duplicate_audit[
        ,
        .(
            pathway,
            preregistered_domain,
            biological_domain,
            NES_Influenza_GPL6884
        )
    ]
)

dup_audit_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_pathway_duplicate_domain_audit_v1.0.tsv"
)

fwrite(
    duplicate_audit,
    dup_audit_file,
    sep="\t"
)

# ============================================================
# 3. CREATE UNIQUE PATHWAY HYPOTHESIS TABLE
# ============================================================
#
# Each biological pathway is tested exactly once.
# Domain assignments are handled separately downstream.
# ============================================================

check_cols <- c(
    "NES_RSV_GPL10558",
    "padj_RSV_GPL10558",
    "NES_RSV_GPL6884",
    "padj_RSV_GPL6884",
    "NES_Influenza_GPL6884",
    "padj_Influenza_GPL6884",
    "NES_Influenza_vs_RSV_GPL6884",
    "padj_Influenza_vs_RSV_GPL6884",
    "RSV_replicated",
    "shared_infection_response",
    "pathogen_differential",
    "pathogen_shift",
    "biological_domain"
)

for (p in unique(hyp58$pathway)) {

    z <- hyp58[
        pathway == p
    ]

    if (nrow(z) > 1) {

        for (cc in check_cols) {

            vals <- unique(
                as.character(
                    z[[cc]]
                )
            )

            if (length(vals) != 1) {

                stop(
                    paste(
                        "Duplicate pathway differs in",
                        cc,
                        ":",
                        p
                    )
                )
            }
        }
    }
}

hyp56 <- hyp58[
    !duplicated(pathway)
][
    ,
    preregistered_domain := NULL
]

if (
    nrow(hyp56) != 56 ||
    uniqueN(hyp56$pathway) != 56
) {
    stop(
        "Failed to construct 56-pathway unique hypothesis table."
    )
}

hyp56_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_56_unique_pathway_hypotheses_v1.0.tsv"
)

fwrite(
    hyp56,
    hyp56_file,
    sep="\t"
)

# Domain map retains all 58 original assignments.
domain_map <- unique(
    hyp58[
        ,
        .(
            pathway,
            preregistered_domain
        )
    ]
)

if (nrow(domain_map) != 58) {
    stop(
        "Expected 58 unique pathway-domain assignments."
    )
}

domain_map_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_58_pathway_domain_assignments_v1.0.tsv"
)

fwrite(
    domain_map,
    domain_map_file,
    sep="\t"
)

# ============================================================
# 4. IDENTIFY REQUIRED COLUMNS
# ============================================================

find_col <- function(
    x,
    candidates,
    label
) {

    hit <- candidates[
        candidates %in% names(x)
    ]

    if (length(hit) == 0) {

        stop(
            paste0(
                "Could not identify ",
                label,
                ". Available columns:\n",
                paste(
                    names(x),
                    collapse=", "
                )
            )
        )
    }

    hit[1]
}

res_id <- find_col(
    res,
    c(
        "SeqId",
        "seqid",
        "SeqID",
        "assay_id"
    ),
    "assay identifier"
)

ann_id <- find_col(
    ann,
    c(
        "SeqId",
        "seqid",
        "SeqID",
        "assay_id"
    ),
    "annotation assay identifier"
)

t_col <- find_col(
    res,
    c(
        "t",
        "moderated_t",
        "t_statistic",
        "tstat"
    ),
    "moderated t-statistic"
)

symbol_col <- find_col(
    ann,
    c(
        "EntrezGeneSymbol",
        "gene_symbol",
        "GeneSymbol"
    ),
    "gene-symbol annotation"
)

cat("\n=== COLUMN RESOLUTION ===\n")
cat("Result ID:", res_id, "\n")
cat("Annotation ID:", ann_id, "\n")
cat("t-statistic:", t_col, "\n")
cat("Gene symbol:", symbol_col, "\n\n")

# ============================================================
# 5. JOIN STATISTICS TO ANNOTATION
# ============================================================

ann_small <- ann[
    ,
    .(
        assay_join_id=get(ann_id),
        raw_gene_symbol=get(symbol_col)
    )
]

tmp <- merge(
    res,
    ann_small,
    by.x=res_id,
    by.y="assay_join_id",
    all.x=TRUE
)

if (nrow(tmp) != nrow(res)) {
    stop(
        "Annotation join changed assay-row count."
    )
}

# ============================================================
# 6. EXPAND COMPOUND SYMBOLS
# ============================================================

tmp[
    ,
    raw_gene_symbol :=
        as.character(
            raw_gene_symbol
        )
]

tmp <- tmp[
    !is.na(raw_gene_symbol) &
    trimws(raw_gene_symbol) != ""
]

id_cols <- setdiff(
    names(tmp),
    "raw_gene_symbol"
)

expanded <- tmp[
    ,
    .(
        gene_symbol=trimws(
            unlist(
                strsplit(
                    raw_gene_symbol,
                    "\\|"
                )
            )
        )
    ),
    by=id_cols
]

expanded <- expanded[
    !is.na(gene_symbol) &
    gene_symbol != ""
]

# ============================================================
# 7. COLLAPSE TO GENE-LEVEL PROTEOMIC RANK
# ============================================================
#
# Predeclared rule:
# median moderated limma t across all assays per gene.
# ============================================================

gene_rank <- expanded[
    ,
    .(
        protein_median_t=
            median(
                get(t_col),
                na.rm=TRUE
            ),

        n_assays=
            uniqueN(
                get(res_id)
            )
    ),
    by=gene_symbol
]

gene_rank <- gene_rank[
    is.finite(
        protein_median_t
    )
]

setorder(
    gene_rank,
    -protein_median_t,
    gene_symbol
)

if (
    anyDuplicated(
        gene_rank$gene_symbol
    )
) {
    stop(
        "Duplicated gene symbols after gene collapse."
    )
}

ranks <- gene_rank$protein_median_t
names(ranks) <- gene_rank$gene_symbol

ranks <- sort(
    ranks,
    decreasing=TRUE
)

gene_rank_file <- file.path(
    OUTDIR,
    "Influenza-P01_proteome_gene_level_median_t_rank_v1.0.tsv"
)

fwrite(
    gene_rank,
    gene_rank_file,
    sep="\t"
)

cat("=== GENE-LEVEL RANK ===\n")
cat(
    "Unique ranked genes:",
    length(ranks),
    "\n"
)

cat(
    "Genes represented by >1 assay:",
    sum(
        gene_rank$n_assays > 1
    ),
    "\n\n"
)

# ============================================================
# 8. RECOVER MEMBERSHIP FOR 56 UNIQUE PATHWAYS
# ============================================================

membership <- gs[
    gs_name %in%
        hyp56$pathway,
    .(
        pathway=gs_name,
        gene_symbol
    )
]

membership <- unique(
    membership
)

found <- unique(
    membership$pathway
)

missing_pathways <- setdiff(
    hyp56$pathway,
    found
)

cat("=== PATHWAY MEMBERSHIP ===\n")
cat(
    "Unique frozen pathways:",
    nrow(hyp56),
    "\n"
)
cat(
    "Recovered:",
    length(found),
    "\n"
)
cat(
    "Missing:",
    length(missing_pathways),
    "\n"
)

if (
    length(missing_pathways) > 0
) {

    print(
        missing_pathways
    )

    stop(
        "Frozen pathway membership incomplete."
    )
}

membership_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_56_unique_pathway_membership_v1.0.tsv"
)

fwrite(
    membership,
    membership_file,
    sep="\t"
)

# ============================================================
# 9. PROTEOMIC COVERAGE
# ============================================================

coverage <- membership[
    ,
    .(
        frozen_pathway_size=
            uniqueN(
                gene_symbol
            ),

        measured_genes=
            uniqueN(
                gene_symbol[
                    gene_symbol %in%
                    names(ranks)
                ]
            )
    ),
    by=pathway
]

coverage[
    ,
    coverage_fraction :=
        measured_genes /
        frozen_pathway_size
]

coverage[
    ,
    coverage_percent :=
        100 *
        coverage_fraction
]

coverage <- merge(
    coverage,
    hyp56[
        ,
        .(
            pathway,
            transcript_NES=
                NES_Influenza_GPL6884
        )
    ],
    by="pathway",
    all.x=TRUE
)

coverage_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_56_pathway_proteomic_coverage_v1.0.tsv"
)

fwrite(
    coverage,
    coverage_file,
    sep="\t"
)

# ============================================================
# 10. BUILD PATHWAY LIST
# ============================================================

pathways <- split(
    membership$gene_symbol,
    membership$pathway
)

pathways <- lapply(
    pathways,
    function(x) {

        unique(
            x[
                x %in%
                names(ranks)
            ]
        )
    }
)

measured_sizes <- lengths(
    pathways
)

testable <- (
    measured_sizes >= 10
)

cat("\n=== PATHWAY TESTABILITY ===\n")

cat(
    "Unique frozen pathways:",
    length(pathways),
    "\n"
)

cat(
    "Testable at >=10 measured genes:",
    sum(testable),
    "\n"
)

cat(
    "Not testable:",
    sum(!testable),
    "\n"
)

print(
    summary(
        measured_sizes
    )
)

testability <- data.table(
    pathway=names(
        measured_sizes
    ),
    measured_size=as.integer(
        measured_sizes
    ),
    primary_testable=as.logical(
        testable
    )
)

testability_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_56_pathway_testability_v1.0.tsv"
)

fwrite(
    testability,
    testability_file,
    sep="\t"
)

# ============================================================
# 11. FGSEA
# ============================================================

set.seed(
    20260812
)

fg <- fgseaMultilevel(
    pathways=
        pathways[testable],

    stats=
        ranks,

    minSize=
        10,

    maxSize=
        1000,

    eps=
        0
)

fg <- as.data.table(
    fg
)

fg[
    ,
    leadingEdge :=
        vapply(
            leadingEdge,
            function(x) {
                paste(
                    x,
                    collapse=";"
                )
            },
            character(1)
        )
]

# ============================================================
# 12. PRIMARY BH FAMILY = UNIQUE TESTED PATHWAYS
# ============================================================

fg[
    ,
    proteomic_padj_primary :=
        p.adjust(
            pval,
            method="BH"
        )
]

fg <- merge(
    fg,
    hyp56[
        ,
        .(
            pathway,
            transcript_NES=
                NES_Influenza_GPL6884,

            transcript_padj=
                padj_Influenza_GPL6884,

            biological_domain,

            shared_infection_response,

            pathogen_differential,

            pathogen_shift
        )
    ],
    by="pathway",
    all.x=TRUE
)

fg[
    ,
    direction_concordant :=
        sign(NES) ==
        sign(
            transcript_NES
        )
]

fg[
    ,
    proteomic_FDR05 :=
        proteomic_padj_primary <
        0.05
]

fg[
    ,
    concordant_FDR05 :=
        direction_concordant &
        proteomic_FDR05
]

fg[
    ,
    discordant_FDR05 :=
        (!direction_concordant) &
        proteomic_FDR05
]

fg[
    ,
    abs_NES_for_sort := abs(NES)
]

setorder(
    fg,
    proteomic_padj_primary,
    -abs_NES_for_sort
)

fg[
    ,
    abs_NES_for_sort := NULL
]

fg_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_56_unique_pathway_fgsea_v1.0.tsv"
)

fwrite(
    fg,
    fg_file,
    sep="\t"
)

# ============================================================
# 13. PRIMARY PATHWAY VALIDATION STATISTICS
# ============================================================

direction_test <- binom.test(
    sum(
        fg$direction_concordant
    ),
    nrow(fg),
    p=0.5,
    alternative="greater"
)

rho_test <- suppressWarnings(
    cor.test(
        fg$transcript_NES,
        fg$NES,
        method="spearman",
        exact=FALSE
    )
)

summary_table <- data.table(
    metric=c(
        "frozen_rows",
        "unique_frozen_pathways",
        "testable_unique_pathways",
        "direction_concordant",
        "direction_concordance_fraction",
        "direction_binomial_p",
        "NES_spearman_rho",
        "NES_spearman_p",
        "proteomic_FDR05",
        "concordant_FDR05",
        "discordant_FDR05"
    ),

    value=c(
        58,
        56,
        nrow(fg),

        sum(
            fg$direction_concordant
        ),

        mean(
            fg$direction_concordant
        ),

        direction_test$p.value,

        unname(
            rho_test$estimate
        ),

        rho_test$p.value,

        sum(
            fg$proteomic_FDR05
        ),

        sum(
            fg$concordant_FDR05
        ),

        sum(
            fg$discordant_FDR05
        )
    )
)

summary_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_pathway_validation_summary_v1.0.tsv"
)

fwrite(
    summary_table,
    summary_file,
    sep="\t"
)

# ============================================================
# 14. DOMAIN SUMMARY
# ============================================================
#
# Re-expand unique tested pathway results to the original
# 58 pathway-domain assignments.
#
# Thus the two cross-domain pathways legitimately contribute
# to both Metabolism and Translation summaries, while still
# being tested only once statistically.
# ============================================================

domain_results <- merge(
    domain_map,
    fg,
    by="pathway",
    all.x=FALSE,
    all.y=FALSE,
    allow.cartesian=TRUE
)

domain_summary <- domain_results[
    ,
    .(
        n_pathway_assignments=
            .N,

        n_unique_pathways=
            uniqueN(
                pathway
            ),

        n_direction_concordant=
            sum(
                direction_concordant
            ),

        direction_concordance_fraction=
            mean(
                direction_concordant
            ),

        n_proteomic_FDR05=
            sum(
                proteomic_FDR05
            ),

        n_concordant_FDR05=
            sum(
                concordant_FDR05
            ),

        n_discordant_FDR05=
            sum(
                discordant_FDR05
            ),

        median_transcript_NES=
            median(
                transcript_NES,
                na.rm=TRUE
            ),

        median_protein_NES=
            median(
                NES,
                na.rm=TRUE
            )
    ),
    by=preregistered_domain
]

setorder(
    domain_summary,
    -n_concordant_FDR05,
    -direction_concordance_fraction
)

domain_file <- file.path(
    OUTDIR,
    "Influenza-P01_FROZEN_7_domain_validation_summary_v1.0.tsv"
)

fwrite(
    domain_summary,
    domain_file,
    sep="\t"
)

# ============================================================
# 15. REPORT
# ============================================================

cat("\n=== PRIMARY PATHWAY VALIDATION ===\n")

print(
    summary_table
)

cat("\n=== DOMAIN SUMMARY ===\n")

print(
    domain_summary
)

cat("\n=== CONCORDANT FDR-SUPPORTED PATHWAYS ===\n")

concordant <- fg[
    concordant_FDR05 == TRUE,
    .(
        pathway,
        transcript_NES,
        NES,
        proteomic_padj_primary
    )
]

if (
    nrow(concordant) > 0
) {

    print(
        concordant
    )

} else {

    cat("None\n")
}

cat("\n=== DISCORDANT FDR-SUPPORTED PATHWAYS ===\n")

discordant <- fg[
    discordant_FDR05 == TRUE,
    .(
        pathway,
        transcript_NES,
        NES,
        proteomic_padj_primary
    )
]

if (
    nrow(discordant) > 0
) {

    print(
        discordant
    )

} else {

    cat("None\n")
}

cat("\n=== OUTPUT FILES ===\n")

for (
    f in c(
        dup_audit_file,
        hyp56_file,
        domain_map_file,
        gene_rank_file,
        membership_file,
        coverage_file,
        testability_file,
        fg_file,
        summary_file,
        domain_file
    )
) {

    cat(
        f,
        "\n"
    )
}

cat(
    "\nSTATUS: FROZEN PATHWAY VALIDATION COMPLETE\n"
)
