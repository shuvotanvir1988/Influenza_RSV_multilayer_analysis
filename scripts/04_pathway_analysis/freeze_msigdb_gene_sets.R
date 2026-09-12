#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(msigdbr)
  library(dplyr)
  library(data.table)
})

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/gene_sets"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("SESSION 17 — FREEZE MSIGDB GENE SET RESOURCES\n")
cat("============================================================\n\n")

cat("R version:\n")
cat(R.version.string, "\n\n")

cat("Package versions:\n")
cat("msigdbr:", as.character(packageVersion("msigdbr")), "\n")
cat("fgsea:",
    if (requireNamespace("fgsea", quietly = TRUE))
      as.character(packageVersion("fgsea"))
    else
      "NOT INSTALLED",
    "\n")
cat("data.table:", as.character(packageVersion("data.table")), "\n")
cat("dplyr:", as.character(packageVersion("dplyr")), "\n\n")

# ------------------------------------------------------------
# Inspect available MSigDB collections
# ------------------------------------------------------------

collections <- msigdbr_collections()

collections_file <- file.path(
  outdir,
  "MSigDB_available_collections.tsv"
)

fwrite(
  as.data.table(collections),
  collections_file,
  sep = "\t"
)

cat("Available collections saved:\n")
cat(collections_file, "\n\n")

# ------------------------------------------------------------
# Retrieve human gene sets
# ------------------------------------------------------------

hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

reactome <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
)

gobp <- msigdbr(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
)

# ------------------------------------------------------------
# Keep reproducibility-relevant fields
# ------------------------------------------------------------

keep_columns <- c(
  "gs_name",
  "gene_symbol",
  "ncbi_gene",
  "ensembl_gene",
  "gs_collection",
  "gs_subcollection",
  "gs_id",
  "gs_description",
  "gs_url",
  "db_version"
)

trim_columns <- function(x) {

  available <- intersect(
    keep_columns,
    colnames(x)
  )

  x %>%
    select(all_of(available)) %>%
    distinct()
}

hallmark_out <- trim_columns(hallmark)
reactome_out <- trim_columns(reactome)
gobp_out <- trim_columns(gobp)

# ------------------------------------------------------------
# Write frozen local snapshots
# ------------------------------------------------------------

hallmark_file <- file.path(
  outdir,
  "MSigDB_Hallmark_Homo_sapiens.tsv.gz"
)

reactome_file <- file.path(
  outdir,
  "MSigDB_Reactome_Homo_sapiens.tsv.gz"
)

gobp_file <- file.path(
  outdir,
  "MSigDB_GO_BP_Homo_sapiens.tsv.gz"
)

fwrite(
  as.data.table(hallmark_out),
  hallmark_file,
  sep = "\t"
)

fwrite(
  as.data.table(reactome_out),
  reactome_file,
  sep = "\t"
)

fwrite(
  as.data.table(gobp_out),
  gobp_file,
  sep = "\t"
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summarize_collection <- function(name, x) {

  version_value <- if ("db_version" %in% colnames(x)) {
    paste(unique(x$db_version), collapse = ";")
  } else {
    NA_character_
  }

  data.frame(
    collection = name,
    rows = nrow(x),
    unique_pathways = n_distinct(x$gs_name),
    unique_genes = n_distinct(x$gene_symbol),
    db_version = version_value,
    stringsAsFactors = FALSE
  )
}

summary <- bind_rows(
  summarize_collection("Hallmark", hallmark_out),
  summarize_collection("Reactome", reactome_out),
  summarize_collection("GO_BP", gobp_out)
)

summary_file <- file.path(
  outdir,
  "MSigDB_gene_set_freeze_summary.tsv"
)

fwrite(
  as.data.table(summary),
  summary_file,
  sep = "\t"
)

cat("============================================================\n")
cat("GENE SET FREEZE SUMMARY\n")
cat("============================================================\n")

print(summary)

cat("\nSaved:\n")
cat(hallmark_file, "\n")
cat(reactome_file, "\n")
cat(gobp_file, "\n")
cat(summary_file, "\n")

cat("\nSESSION 17 MSIGDB FREEZE COMPLETE\n")
