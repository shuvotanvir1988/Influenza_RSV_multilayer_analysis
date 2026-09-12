#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

# ============================================================
# SESSION 17 — HALLMARK PRERANKED ENRICHMENT
# ============================================================

set.seed(1701)

project_root <- function() {
  env_root <- Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = "")

  if (nzchar(env_root)) {
    env_root <- normalizePath(
      env_root,
      winslash = "/",
      mustWork = TRUE
    )

    if (
      dir.exists(file.path(env_root, "scripts")) &&
      dir.exists(file.path(env_root, "results"))
    ) {
      return(env_root)
    }

    stop(
      "INFLUENZA_RSV_PROJECT_ROOT does not point to a valid project root."
    )
  }

  current <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )

  repeat {
    if (
      dir.exists(file.path(current, "scripts")) &&
      dir.exists(file.path(current, "results"))
    ) {
      return(current)
    }

    parent <- dirname(current)

    if (identical(parent, current)) {
      break
    }

    current <- parent
  }

  stop(
    paste(
      "Could not locate Influenza_RSV_Project root.",
      "Run from inside the project directory or set",
      "INFLUENZA_RSV_PROJECT_ROOT."
    )
  )
}

root <- project_root()

rank_dir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/ranked_lists"
)

geneset_file <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/gene_sets",
  "MSigDB_Hallmark_Homo_sapiens.tsv.gz"
)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/fgsea_hallmark"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Frozen analysis decisions
# ------------------------------------------------------------

MIN_SIZE <- 15
MAX_SIZE <- 500
FDR_THRESHOLD <- 0.05
RANDOM_SEED <- 1701
RANKING_METRIC <- "moderated_limma_t"
METHOD <- "fgseaMultilevel"
COLLECTION <- "MSigDB_Hallmark"

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

cat("============================================================\n")
cat("SESSION 17 — HALLMARK PRERANKED ENRICHMENT\n")
cat("============================================================\n")

cat("\nFrozen parameters:\n")
cat("Ranking metric:", RANKING_METRIC, "\n")
cat("Method:", METHOD, "\n")
cat("Minimum pathway size:", MIN_SIZE, "\n")
cat("Maximum pathway size:", MAX_SIZE, "\n")
cat("Primary FDR threshold:", FDR_THRESHOLD, "\n")
cat("MSigDB collection:", COLLECTION, "\n")
cat("Random seed:", RANDOM_SEED, "\n")

# ============================================================
# Validate required files
# ============================================================

if (!file.exists(geneset_file)) {
  stop(
    "Frozen Hallmark gene-set file not found:\n",
    geneset_file
  )
}

for (nm in names(rank_files)) {

  path <- file.path(
    rank_dir,
    rank_files[[nm]]
  )

  if (!file.exists(path)) {
    stop(
      "Ranked gene-list file not found for ",
      nm,
      ":\n",
      path
    )
  }
}

# ============================================================
# Load frozen Hallmark gene sets
# ============================================================

gs <- fread(
  geneset_file
)

required_gs_cols <- c(
  "gs_name",
  "gene_symbol"
)

missing_gs_cols <- setdiff(
  required_gs_cols,
  names(gs)
)

if (length(missing_gs_cols) > 0) {
  stop(
    "Gene-set file missing columns: ",
    paste(missing_gs_cols, collapse = ", ")
  )
}

gs <- gs[
  !is.na(gs_name) &
  gs_name != "" &
  !is.na(gene_symbol) &
  gene_symbol != ""
]

pathways <- split(
  gs$gene_symbol,
  gs$gs_name
)

pathways <- lapply(
  pathways,
  unique
)

cat(
  "\nFrozen Hallmark pathways loaded:",
  length(pathways),
  "\n"
)

if (length(pathways) != 50) {
  warning(
    "Expected 50 Hallmark pathways, observed ",
    length(pathways)
  )
}

# ============================================================
# Output containers
# ============================================================

summary_rows <- list()

# ============================================================
# Run enrichment independently for each contrast
# ============================================================

for (contrast_name in names(rank_files)) {

  cat("\n============================================================\n")
  cat(contrast_name, "\n")
  cat("============================================================\n")

  rank_file <- file.path(
    rank_dir,
    rank_files[[contrast_name]]
  )

  x <- fread(
    rank_file
  )

  required_rank_cols <- c(
    "gene_symbol",
    "t"
  )

  missing_rank_cols <- setdiff(
    required_rank_cols,
    names(x)
  )

  if (length(missing_rank_cols) > 0) {
    stop(
      contrast_name,
      ": ranked file missing columns: ",
      paste(missing_rank_cols, collapse = ", ")
    )
  }

  # ----------------------------------------------------------
  # Clean ranking
  # ----------------------------------------------------------

  x[, gene_symbol := trimws(
    as.character(gene_symbol)
  )]

  x[, t := as.numeric(t)]

  x <- x[
    !is.na(gene_symbol) &
    gene_symbol != "" &
    !is.na(t) &
    is.finite(t)
  ]

  if (anyDuplicated(x$gene_symbol)) {
    duplicated_genes <- unique(
      x$gene_symbol[
        duplicated(x$gene_symbol)
      ]
    )

    stop(
      contrast_name,
      ": duplicated gene symbols detected: ",
      paste(
        head(duplicated_genes, 20),
        collapse = ", "
      )
    )
  }

  stats <- x$t
  names(stats) <- x$gene_symbol

  stats <- sort(
    stats,
    decreasing = TRUE
  )

  cat(
    "Ranked genes:",
    length(stats),
    "\n"
  )

  cat(
    "Maximum t:",
    max(stats),
    "\n"
  )

  cat(
    "Minimum t:",
    min(stats),
    "\n"
  )

  # ----------------------------------------------------------
  # Pathway coverage within this ranked universe
  # ----------------------------------------------------------

  measured_sizes <- vapply(
    pathways,
    function(g) {
      length(
        intersect(
          g,
          names(stats)
        )
      )
    },
    numeric(1)
  )

  testable_before_fgsea <- sum(
    measured_sizes >= MIN_SIZE &
    measured_sizes <= MAX_SIZE
  )

  cat(
    "Hallmark pathways passing size rule:",
    testable_before_fgsea,
    "\n"
  )

  # ----------------------------------------------------------
  # Run fgseaMultilevel
  # ----------------------------------------------------------

  set.seed(RANDOM_SEED)

  result <- fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    eps = 0
  )

  result <- as.data.table(
    result
  )

  if (nrow(result) == 0) {
    stop(
      contrast_name,
      ": fgsea returned zero pathways."
    )
  }

  # ----------------------------------------------------------
  # Convert leading-edge list to text
  # ----------------------------------------------------------

  result[, leadingEdge := vapply(
    leadingEdge,
    function(z) {
      paste(
        z,
        collapse = ";"
      )
    },
    character(1)
  )]

  # ----------------------------------------------------------
  # Derived fields
  # ----------------------------------------------------------

  result[, abs_NES := abs(NES)]

  result[, enrichment_direction :=
    fifelse(
      NES > 0,
      "positive",
      fifelse(
        NES < 0,
        "negative",
        "zero"
      )
    )
  ]

  result[, significant_FDR05 :=
    !is.na(padj) &
    padj < FDR_THRESHOLD
  ]

  # ----------------------------------------------------------
  # Stable sorting
  #
  # Primary: smallest FDR
  # Secondary: largest absolute NES
  # ----------------------------------------------------------

  setorder(
    result,
    padj,
    -abs_NES,
    pathway
  )

  # ----------------------------------------------------------
  # Save complete result
  # ----------------------------------------------------------

  result_file <- file.path(
    outdir,
    paste0(
      contrast_name,
      "_Hallmark_fgsea.tsv"
    )
  )

  fwrite(
    result,
    result_file,
    sep = "\t"
  )

  # ----------------------------------------------------------
  # Save significant-only result
  # ----------------------------------------------------------

  significant_result <- result[
    significant_FDR05 == TRUE
  ]

  significant_file <- file.path(
    outdir,
    paste0(
      contrast_name,
      "_Hallmark_fgsea_FDR05.tsv"
    )
  )

  fwrite(
    significant_result,
    significant_file,
    sep = "\t"
  )

  # ----------------------------------------------------------
  # Summary statistics
  # ----------------------------------------------------------

  n_tested <- nrow(result)

  n_sig <- sum(
    result$significant_FDR05,
    na.rm = TRUE
  )

  n_positive_sig <- sum(
    result$significant_FDR05 &
    result$NES > 0,
    na.rm = TRUE
  )

  n_negative_sig <- sum(
    result$significant_FDR05 &
    result$NES < 0,
    na.rm = TRUE
  )

  n_na_padj <- sum(
    is.na(result$padj)
  )

  cat(
    "Hallmark pathways tested:",
    n_tested,
    "\n"
  )

  cat(
    "Significant FDR < 0.05:",
    n_sig,
    "\n"
  )

  cat(
    "  Positive NES:",
    n_positive_sig,
    "\n"
  )

  cat(
    "  Negative NES:",
    n_negative_sig,
    "\n"
  )

  cat(
    "Pathways with NA adjusted P value:",
    n_na_padj,
    "\n"
  )

  # ----------------------------------------------------------
  # Show top pathways
  # ----------------------------------------------------------

  cat("\nTop pathways by FDR:\n")

  top_n <- min(
    15,
    nrow(result)
  )

  print(
    result[
      1:top_n,
      .(
        pathway,
        size,
        NES,
        abs_NES,
        pval,
        padj,
        significant_FDR05
      )
    ]
  )

  # ----------------------------------------------------------
  # Show strongest positive NES
  # ----------------------------------------------------------

  positive_view <- copy(
    result[
      !is.na(NES) &
      NES > 0
    ]
  )

  setorder(
    positive_view,
    -NES
  )

  cat("\nStrongest positive NES pathways:\n")

  if (nrow(positive_view) > 0) {

    print(
      positive_view[
        1:min(10, .N),
        .(
          pathway,
          size,
          NES,
          padj
        )
      ]
    )

  } else {

    cat("None\n")
  }

  # ----------------------------------------------------------
  # Show strongest negative NES
  # ----------------------------------------------------------

  negative_view <- copy(
    result[
      !is.na(NES) &
      NES < 0
    ]
  )

  setorder(
    negative_view,
    NES
  )

  cat("\nStrongest negative NES pathways:\n")

  if (nrow(negative_view) > 0) {

    print(
      negative_view[
        1:min(10, .N),
        .(
          pathway,
          size,
          NES,
          padj
        )
      ]
    )

  } else {

    cat("None\n")
  }

  cat("\nSaved:\n")
  cat(result_file, "\n")
  cat(significant_file, "\n")

  # ----------------------------------------------------------
  # Add summary row
  # ----------------------------------------------------------

  summary_rows[[length(summary_rows) + 1]] <-
    data.frame(

      contrast = contrast_name,

      ranked_genes =
        length(stats),

      pathways_loaded =
        length(pathways),

      pathways_passing_size_rule =
        testable_before_fgsea,

      pathways_tested =
        n_tested,

      significant_FDR05 =
        n_sig,

      significant_positive_NES =
        n_positive_sig,

      significant_negative_NES =
        n_negative_sig,

      NA_adjusted_p_values =
        n_na_padj,

      minSize =
        MIN_SIZE,

      maxSize =
        MAX_SIZE,

      FDR_threshold =
        FDR_THRESHOLD,

      ranking_metric =
        RANKING_METRIC,

      method =
        METHOD,

      collection =
        COLLECTION,

      random_seed =
        RANDOM_SEED,

      stringsAsFactors = FALSE
    )
}

# ============================================================
# Combined summary
# ============================================================

summary <- rbindlist(
  summary_rows,
  fill = TRUE
)

summary_file <- file.path(
  outdir,
  "DS001_GSE38900_Hallmark_fgsea_summary.tsv"
)

fwrite(
  summary,
  summary_file,
  sep = "\t"
)

cat("\n============================================================\n")
cat("HALLMARK FGSEA SUMMARY\n")
cat("============================================================\n")

print(
  summary
)

cat("\nSaved summary:\n")
cat(summary_file, "\n")

# ============================================================
# Session metadata
# ============================================================

metadata_file <- file.path(
  outdir,
  "DS001_GSE38900_Hallmark_fgsea_parameters.tsv"
)

metadata <- data.table(
  parameter = c(
    "analysis_session",
    "ranking_metric",
    "fgsea_method",
    "minSize",
    "maxSize",
    "FDR_threshold",
    "MSigDB_collection",
    "random_seed",
    "R_version",
    "fgsea_version"
  ),

  value = c(
    "Session_17",
    RANKING_METRIC,
    METHOD,
    as.character(MIN_SIZE),
    as.character(MAX_SIZE),
    as.character(FDR_THRESHOLD),
    COLLECTION,
    as.character(RANDOM_SEED),
    R.version.string,
    as.character(
      packageVersion("fgsea")
    )
  )
)

fwrite(
  metadata,
  metadata_file,
  sep = "\t"
)

cat("\nSaved parameter record:\n")
cat(metadata_file, "\n")

cat("\n============================================================\n")
cat("SESSION 17 HALLMARK ENRICHMENT COMPLETE\n")
cat("============================================================\n")
