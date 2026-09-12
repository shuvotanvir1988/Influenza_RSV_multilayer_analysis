#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

rep_file <- file.path(
    root,
    "results/DS001_GSE38900/pathway_analysis/tables/preregistered_domains",
    "DS001_GSE38900_preregistered_domain_representative_pathways.tsv"
)

reactome_file <- file.path(
    root,
    "results/DS001_GSE38900/pathway_analysis/gene_sets",
    "MSigDB_Reactome_Homo_sapiens.tsv.gz"
)

coverage_file <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/audit",
    "DS001_GSE38900_representative_pathway_gene_coverage.tsv"
)

expr_files <- c(
    GPL10558 = file.path(
        root,
        "results/DS001_GSE38900/cell_composition/audit",
        "DS001_GSE38900_GPL10558_cell_composition_input_HUGO.tsv.gz"
    ),
    GPL6884 = file.path(
        root,
        "results/DS001_GSE38900/cell_composition/audit",
        "DS001_GSE38900_GPL6884_cell_composition_input_HUGO.tsv.gz"
    )
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/tables"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

target_domains <- c(
    "Interferon",
    "Antigen_Presentation",
    "Neutrophil_Monocyte",
    "Adaptive_Immunity",
    "Metabolism",
    "Translation"
)

rep <- read.delim(rep_file, check.names = FALSE)
rep <- rep[
    rep$preregistered_domain %in% target_domains,
    ,
    drop = FALSE
]

rep <- unique(
    rep[, c(
        "pathway",
        "biological_domain",
        "preregistered_domain"
    )]
)

rep$pathway <- toupper(trimws(rep$pathway))

reactome <- read.delim(
    gzfile(reactome_file),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

reactome$gs_name <- toupper(trimws(reactome$gs_name))
reactome$gene_symbol <- toupper(trimws(reactome$gene_symbol))

coverage <- read.delim(
    coverage_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

all_scores <- list()
all_domains <- list()

for (platform in names(expr_files)) {

    message("Scoring ", platform)

    d <- read.delim(
        gzfile(expr_files[[platform]]),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    genes <- toupper(trimws(d$gene_symbol))
    expr <- as.matrix(
        d[, setdiff(colnames(d), "gene_symbol"), drop = FALSE]
    )

    storage.mode(expr) <- "numeric"
    rownames(expr) <- genes

    # Gene-wise z-score across samples
    z <- t(scale(t(expr)))

    z[!is.finite(z)] <- NA_real_

    platform_scores <- list()

    for (i in seq_len(nrow(rep))) {

        pw <- rep$pathway[i]

        geneset <- unique(
            reactome$gene_symbol[
                reactome$gs_name == pw
            ]
        )

        genes_present <- intersect(
            geneset,
            rownames(z)
        )

        if (length(genes_present) < 3) {
            next
        }

        score <- colMeans(
            z[genes_present, , drop = FALSE],
            na.rm = TRUE
        )

        covrow <- coverage[
            coverage$platform == platform &
            coverage$pathway == pw,
            ,
            drop = FALSE
        ]

        coverage_pct <- if (nrow(covrow) == 1) {
            covrow$coverage_percent
        } else {
            NA_real_
        }

        platform_scores[[pw]] <- data.frame(
            platform = platform,
            pathway = pw,
            preregistered_domain =
                rep$preregistered_domain[i],
            biological_domain =
                rep$biological_domain[i],
            coverage_percent = coverage_pct,
            sample = names(score),
            pathway_score = as.numeric(score),
            high_coverage = coverage_pct >= 80,
            stringsAsFactors = FALSE
        )
    }

    ps <- do.call(rbind, platform_scores)
    rownames(ps) <- NULL

    all_scores[[platform]] <- ps

    # Domain-level score:
    # mean of high-coverage pathway scores within each sample/domain.
    high <- ps[ps$high_coverage, , drop = FALSE]

    domain_score <- aggregate(
        pathway_score ~
            platform +
            preregistered_domain +
            sample,
        data = high,
        FUN = mean,
        na.rm = TRUE
    )

    colnames(domain_score)[
        colnames(domain_score) == "pathway_score"
    ] <- "domain_score"

    all_domains[[platform]] <- domain_score
}

pathway_scores <- do.call(rbind, all_scores)
domain_scores <- do.call(rbind, all_domains)

write.table(
    pathway_scores,
    file.path(
        outdir,
        "DS001_GSE38900_representative_pathway_sample_scores.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    domain_scores,
    file.path(
        outdir,
        "DS001_GSE38900_domain_sample_scores.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n=== PATHWAY SCORE COUNTS ===\n")
print(
    with(
        pathway_scores,
        table(platform, preregistered_domain)
    )
)

cat("\n=== DOMAIN SCORE COUNTS ===\n")
print(
    with(
        domain_scores,
        table(platform, preregistered_domain)
    )
)

cat("\nScoring complete.\n")
