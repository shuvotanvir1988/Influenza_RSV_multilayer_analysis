#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

rank_dir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/ranked_lists"
)

geneset_dir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/gene_sets"
)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("SESSION 17 — GENE SET COVERAGE AUDIT\n")
cat("============================================================\n")

rank_files <- c(

  GPL10558_RSV_acute_vs_control =
    "GPL10558_RSV_acute_vs_control_ranked_by_limma_t.tsv",

  GPL6884_RSV_acute_vs_control =
    "GPL6884_RSV_acute_vs_control_ranked_by_limma_t.tsv",

  GPL6884_InfluenzaA_vs_control =
    "GPL6884_InfluenzaA_vs_control_ranked_by_limma_t.tsv",

  GPL6884_InfluenzaA_vs_RSVacute =
    "GPL6884_InfluenzaA_vs_RSVacute_ranked_by_limma_t.tsv"
)

geneset_files <- c(

  Hallmark =
    "MSigDB_Hallmark_Homo_sapiens.tsv.gz",

  Reactome =
    "MSigDB_Reactome_Homo_sapiens.tsv.gz",

  GO_BP =
    "MSigDB_GO_BP_Homo_sapiens.tsv.gz"
)

# ------------------------------------------------------------
# Load ranked universes
# ------------------------------------------------------------

ranked_universes <- list()

for (nm in names(rank_files)) {

  path <- file.path(
    rank_dir,
    rank_files[[nm]]
  )

  x <- fread(path)

  genes <- unique(
    trimws(
      as.character(x$gene_symbol)
    )
  )

  genes <- genes[
    !is.na(genes) &
    genes != ""
  ]

  ranked_universes[[nm]] <- genes
}

# ------------------------------------------------------------
# Audit each gene-set collection
# ------------------------------------------------------------

summary_rows <- list()
pathway_rows <- list()

for (collection_name in names(geneset_files)) {

  gs_path <- file.path(
    geneset_dir,
    geneset_files[[collection_name]]
  )

  gs <- fread(gs_path)

  gs <- gs[
    !is.na(gene_symbol) &
    gene_symbol != ""
  ]

  collection_genes <- unique(gs$gene_symbol)

  pathways <- split(
    gs$gene_symbol,
    gs$gs_name
  )

  pathways <- lapply(
    pathways,
    unique
  )

  for (contrast_name in names(ranked_universes)) {

    universe <- ranked_universes[[contrast_name]]

    overlap <- intersect(
      collection_genes,
      universe
    )

    overall_coverage <- (
      length(overlap) /
      length(collection_genes)
    ) * 100

    pathway_stats <- lapply(
      names(pathways),
      function(pathway_name) {

        original_genes <- pathways[[pathway_name]]

        measured <- intersect(
          original_genes,
          universe
        )

        data.frame(
          collection = collection_name,
          contrast = contrast_name,
          pathway = pathway_name,
          genes_in_msigdb = length(original_genes),
          genes_measured = length(measured),
          pathway_coverage_percent =
            100 * length(measured) /
            length(original_genes),
          stringsAsFactors = FALSE
        )
      }
    )

    pathway_stats <- bind_rows(pathway_stats)

    pathway_rows[[length(pathway_rows) + 1]] <-
      pathway_stats

    summary_rows[[length(summary_rows) + 1]] <-
      data.frame(
        collection = collection_name,
        contrast = contrast_name,
        ranked_gene_universe =
          length(universe),
        collection_unique_genes =
          length(collection_genes),
        collection_genes_measured =
          length(overlap),
        collection_coverage_percent =
          overall_coverage,
        total_pathways =
          length(pathways),
        pathways_with_at_least_10_measured_genes =
          sum(pathway_stats$genes_measured >= 10),
        pathways_with_at_least_15_measured_genes =
          sum(pathway_stats$genes_measured >= 15),
        pathways_with_at_least_25_measured_genes =
          sum(pathway_stats$genes_measured >= 25),
        stringsAsFactors = FALSE
      )
  }
}

summary_df <- bind_rows(summary_rows)

pathway_df <- bind_rows(pathway_rows)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

summary_file <- file.path(
  outdir,
  "DS001_GSE38900_gene_set_coverage_summary.tsv"
)

pathway_file <- file.path(
  outdir,
  "DS001_GSE38900_pathway_level_coverage.tsv.gz"
)

fwrite(
  as.data.table(summary_df),
  summary_file,
  sep = "\t"
)

fwrite(
  as.data.table(pathway_df),
  pathway_file,
  sep = "\t"
)

cat("\n============================================================\n")
cat("GENE SET COVERAGE SUMMARY\n")
cat("============================================================\n")

print(summary_df)

cat("\nSaved:\n")
cat(summary_file, "\n")
cat(pathway_file, "\n")

cat("\nSESSION 17 GENE SET COVERAGE AUDIT COMPLETE\n")
