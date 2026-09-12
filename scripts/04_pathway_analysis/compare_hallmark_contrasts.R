#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

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

indir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/fgsea_hallmark"
)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/hallmark_comparison"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

files <- c(
  RSV_GPL10558 =
    "GPL10558_RSV_acute_vs_control_Hallmark_fgsea.tsv",

  RSV_GPL6884 =
    "GPL6884_RSV_acute_vs_control_Hallmark_fgsea.tsv",

  Influenza_GPL6884 =
    "GPL6884_InfluenzaA_vs_control_Hallmark_fgsea.tsv",

  Influenza_vs_RSV_GPL6884 =
    "GPL6884_InfluenzaA_vs_RSVacute_Hallmark_fgsea.tsv"
)

cat("============================================================\n")
cat("SESSION 17 — HALLMARK CONTRAST COMPARISON\n")
cat("============================================================\n")

tables <- list()

for (nm in names(files)) {

  path <- file.path(indir, files[[nm]])

  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }

  x <- fread(path)

  required <- c(
    "pathway",
    "NES",
    "pval",
    "padj"
  )

  missing <- setdiff(required, names(x))

  if (length(missing) > 0) {
    stop(
      nm,
      ": missing columns: ",
      paste(missing, collapse=", ")
    )
  }

  x <- x[
    ,
    .(
      pathway,
      NES,
      pval,
      padj
    )
  ]

  setnames(
    x,
    c("NES", "pval", "padj"),
    paste0(
      c("NES_", "pval_", "padj_"),
      nm
    )
  )

  tables[[nm]] <- x
}

# ------------------------------------------------------------
# Merge all contrasts
# ------------------------------------------------------------

merged <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by="pathway",
      all=TRUE
    )
  },
  tables
)

# ------------------------------------------------------------
# Significance flags
# ------------------------------------------------------------

merged[
  ,
  sig_RSV_GPL10558 :=
    padj_RSV_GPL10558 < 0.05
]

merged[
  ,
  sig_RSV_GPL6884 :=
    padj_RSV_GPL6884 < 0.05
]

merged[
  ,
  sig_Influenza_GPL6884 :=
    padj_Influenza_GPL6884 < 0.05
]

merged[
  ,
  sig_Influenza_vs_RSV_GPL6884 :=
    padj_Influenza_vs_RSV_GPL6884 < 0.05
]

# ------------------------------------------------------------
# RSV cross-platform concordance
# ------------------------------------------------------------

merged[
  ,
  RSV_same_direction :=
    sign(NES_RSV_GPL10558) ==
    sign(NES_RSV_GPL6884)
]

merged[
  ,
  RSV_significant_both :=
    sig_RSV_GPL10558 &
    sig_RSV_GPL6884
]

merged[
  ,
  RSV_replicated :=
    RSV_significant_both &
    RSV_same_direction
]

# ------------------------------------------------------------
# Shared influenza/RSV response
# GPL6884 only: directly comparable platform
# ------------------------------------------------------------

merged[
  ,
  shared_same_direction :=
    sign(NES_RSV_GPL6884) ==
    sign(NES_Influenza_GPL6884)
]

merged[
  ,
  shared_significant_both :=
    sig_RSV_GPL6884 &
    sig_Influenza_GPL6884
]

merged[
  ,
  shared_infection_response :=
    shared_significant_both &
    shared_same_direction
]

# ------------------------------------------------------------
# Differential pathogen response
# ------------------------------------------------------------

merged[
  ,
  pathogen_differential :=
    sig_Influenza_vs_RSV_GPL6884
]

merged[
  ,
  pathogen_preference :=
    fifelse(
      pathogen_differential &
      NES_Influenza_vs_RSV_GPL6884 > 0,
      "Influenza_enriched",
      fifelse(
        pathogen_differential &
        NES_Influenza_vs_RSV_GPL6884 < 0,
        "RSV_enriched",
        "Not_significant"
      )
    )
]

# ------------------------------------------------------------
# Correlations
# ------------------------------------------------------------

rsv_cor <- cor(
  merged$NES_RSV_GPL10558,
  merged$NES_RSV_GPL6884,
  method="spearman",
  use="complete.obs"
)

infection_cor <- cor(
  merged$NES_RSV_GPL6884,
  merged$NES_Influenza_GPL6884,
  method="spearman",
  use="complete.obs"
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary <- data.table(
  metric = c(
    "Hallmark pathways",
    "RSV significant GPL10558",
    "RSV significant GPL6884",
    "RSV significant on both platforms",
    "RSV replicated same direction",
    "RSV NES Spearman correlation",
    "Shared RSV/influenza significant same direction",
    "RSV vs influenza NES Spearman correlation",
    "Significant influenza-vs-RSV pathways",
    "Influenza-enriched differential pathways",
    "RSV-enriched differential pathways"
  ),

  value = c(
    nrow(merged),

    sum(
      merged$sig_RSV_GPL10558,
      na.rm=TRUE
    ),

    sum(
      merged$sig_RSV_GPL6884,
      na.rm=TRUE
    ),

    sum(
      merged$RSV_significant_both,
      na.rm=TRUE
    ),

    sum(
      merged$RSV_replicated,
      na.rm=TRUE
    ),

    rsv_cor,

    sum(
      merged$shared_infection_response,
      na.rm=TRUE
    ),

    infection_cor,

    sum(
      merged$pathogen_differential,
      na.rm=TRUE
    ),

    sum(
      merged$pathogen_preference ==
        "Influenza_enriched",
      na.rm=TRUE
    ),

    sum(
      merged$pathogen_preference ==
        "RSV_enriched",
      na.rm=TRUE
    )
  )
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

merged_file <- file.path(
  outdir,
  "DS001_GSE38900_Hallmark_contrast_comparison.tsv"
)

summary_file <- file.path(
  outdir,
  "DS001_GSE38900_Hallmark_comparison_summary.tsv"
)

fwrite(
  merged,
  merged_file,
  sep="\t"
)

fwrite(
  summary,
  summary_file,
  sep="\t"
)

cat("\n============================================================\n")
cat("HALLMARK COMPARISON SUMMARY\n")
cat("============================================================\n")

print(summary)

cat("\n============================================================\n")
cat("REPLICATED RSV PATHWAYS\n")
cat("============================================================\n")

print(
  merged[
    RSV_replicated == TRUE,
    .(
      pathway,
      NES_RSV_GPL10558,
      NES_RSV_GPL6884,
      padj_RSV_GPL10558,
      padj_RSV_GPL6884
    )
  ][
    order(-abs(NES_RSV_GPL6884))
  ]
)

cat("\n============================================================\n")
cat("PATHOGEN-DIFFERENTIAL PATHWAYS\n")
cat("============================================================\n")

print(
  merged[
    pathogen_differential == TRUE,
    .(
      pathway,
      NES_RSV_GPL6884,
      NES_Influenza_GPL6884,
      NES_Influenza_vs_RSV_GPL6884,
      padj_Influenza_vs_RSV_GPL6884,
      pathogen_preference
    )
  ][
    order(padj_Influenza_vs_RSV_GPL6884)
  ]
)

cat("\nSaved:\n")
cat(merged_file, "\n")
cat(summary_file, "\n")

cat("\n============================================================\n")
cat("SESSION 17 HALLMARK COMPARISON COMPLETE\n")
cat("============================================================\n")
