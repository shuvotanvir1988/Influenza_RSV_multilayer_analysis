#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

infile <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/reactome_comparison",
  "DS001_GSE38900_Reactome_contrast_comparison.tsv.gz"
)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/preregistered_domains"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("SESSION 17 — PREREGISTERED BIOLOGICAL DOMAIN SUMMARY\n")
cat("============================================================\n")

x <- fread(infile)

# ------------------------------------------------------------
# Domain definitions
# ------------------------------------------------------------

domain_rules <- list(

  Interferon = c(
    "INTERFERON",
    "IFN_STIMULATED",
    "ANTIVIRAL"
  ),

  Innate_Cytokine = c(
    "INTERLEUKIN",
    "CYTOKINE",
    "TNF",
    "NFKB",
    "NF_KB",
    "TLR",
    "INFLAMMASOME",
    "NLRP3",
    "INNATE_IMMUNE"
  ),

  Antigen_Presentation = c(
    "ANTIGEN_PROCESSING",
    "CROSS_PRESENTATION",
    "MHC_CLASS",
    "PRESENTATION"
  ),

  Neutrophil_Monocyte = c(
    "NEUTROPHIL",
    "MONOCYTE",
    "MYELOID",
    "DEGRANULATION"
  ),

  Adaptive_Immunity = c(
    "T_CELL",
    "B_CELL",
    "TCR",
    "BCR",
    "IMMUNOGLOBULIN",
    "LYMPHOCYTE"
  ),

  Metabolism = c(
    "METABOLISM",
    "RESPIRATORY_ELECTRON",
    "OXIDATIVE",
    "MITOCHONDRIAL",
    "FATTY_ACID",
    "LIPID",
    "CHOLESTEROL",
    "PEROXISOME",
    "AMINO_ACID"
  ),

  Translation = c(
    "TRANSLATION",
    "RIBOSOM",
    "RRNA",
    "NONSENSE_MEDIATED",
    "NMD",
    "SRP_DEPENDENT",
    "EIF2AK4",
    "PROTEIN_TARGETING"
  )
)

matches_domain <- function(pathway, patterns) {

  p <- toupper(pathway)

  any(
    vapply(
      patterns,
      function(pattern) {
        grepl(pattern, p)
      },
      logical(1)
    )
  )
}

# ------------------------------------------------------------
# Build long domain membership table
# ------------------------------------------------------------

domain_rows <- list()

for (domain_name in names(domain_rules)) {

  patterns <- domain_rules[[domain_name]]

  idx <- vapply(
    x$pathway,
    matches_domain,
    logical(1),
    patterns = patterns
  )

  temp <- copy(
    x[idx]
  )

  if (nrow(temp) == 0) {
    next
  }

  temp[, preregistered_domain := domain_name]

  domain_rows[[length(domain_rows) + 1]] <- temp
}

domain_table <- rbindlist(
  domain_rows,
  fill = TRUE
)

# ------------------------------------------------------------
# Domain-level quantitative summary
# ------------------------------------------------------------

summary_rows <- list()

for (domain_name in names(domain_rules)) {

  d <- domain_table[
    preregistered_domain == domain_name
  ]

  if (nrow(d) == 0) {
    next
  }

  summary_rows[[length(summary_rows) + 1]] <-
    data.table(

      preregistered_domain = domain_name,

      pathways_in_domain = nrow(d),

      RSV_sig_GPL10558 =
        sum(
          d$sig_RSV_GPL10558,
          na.rm = TRUE
        ),

      RSV_sig_GPL6884 =
        sum(
          d$sig_RSV_GPL6884,
          na.rm = TRUE
        ),

      RSV_replicated =
        sum(
          d$RSV_replicated,
          na.rm = TRUE
        ),

      shared_infection_response =
        sum(
          d$shared_infection_response,
          na.rm = TRUE
        ),

      pathogen_differential =
        sum(
          d$pathogen_differential,
          na.rm = TRUE
        ),

      influenza_shifted =
        sum(
          d$pathogen_shift == "Influenza_shifted",
          na.rm = TRUE
        ),

      RSV_shifted =
        sum(
          d$pathogen_shift == "RSV_shifted",
          na.rm = TRUE
        ),

      median_NES_RSV_GPL10558 =
        median(
          d$NES_RSV_GPL10558,
          na.rm = TRUE
        ),

      median_NES_RSV_GPL6884 =
        median(
          d$NES_RSV_GPL6884,
          na.rm = TRUE
        ),

      median_NES_Influenza_GPL6884 =
        median(
          d$NES_Influenza_GPL6884,
          na.rm = TRUE
        ),

      median_NES_Influenza_vs_RSV =
        median(
          d$NES_Influenza_vs_RSV_GPL6884,
          na.rm = TRUE
        )
    )
}

summary <- rbindlist(
  summary_rows,
  fill = TRUE
)

# ------------------------------------------------------------
# Save full domain membership
# ------------------------------------------------------------

membership_file <- file.path(
  outdir,
  "DS001_GSE38900_preregistered_domain_pathways.tsv.gz"
)

summary_file <- file.path(
  outdir,
  "DS001_GSE38900_preregistered_domain_summary.tsv"
)

fwrite(
  domain_table,
  membership_file,
  sep = "\t"
)

fwrite(
  summary,
  summary_file,
  sep = "\t"
)

# ------------------------------------------------------------
# Representative pathways per domain
# ------------------------------------------------------------

representative_rows <- list()

for (domain_name in unique(domain_table$preregistered_domain)) {

  d <- copy(
    domain_table[
      preregistered_domain == domain_name
    ]
  )

  d[, direct_abs_NES :=
      abs(NES_Influenza_vs_RSV_GPL6884)
  ]

  setorder(
    d,
    padj_Influenza_vs_RSV_GPL6884,
    -direct_abs_NES
  )

  representative_rows[[length(representative_rows) + 1]] <-
    d[
      1:min(10, .N)
    ]
}

representative <- rbindlist(
  representative_rows,
  fill = TRUE
)

representative_file <- file.path(
  outdir,
  "DS001_GSE38900_preregistered_domain_representative_pathways.tsv"
)

fwrite(
  representative,
  representative_file,
  sep = "\t"
)

# ------------------------------------------------------------
# Terminal output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PREREGISTERED DOMAIN SUMMARY\n")
cat("============================================================\n")

print(summary)

cat("\n============================================================\n")
cat("REPRESENTATIVE PATHWAYS BY DOMAIN\n")
cat("============================================================\n")

for (domain_name in unique(representative$preregistered_domain)) {

  cat("\n------------------------------\n")
  cat(domain_name, "\n")
  cat("------------------------------\n")

  print(
    representative[
      preregistered_domain == domain_name,
      .(
        pathway,
        NES_RSV_GPL10558,
        NES_RSV_GPL6884,
        NES_Influenza_GPL6884,
        NES_Influenza_vs_RSV_GPL6884,
        padj_Influenza_vs_RSV_GPL6884,
        pathogen_shift
      )
    ]
  )
}

cat("\nSaved:\n")
cat(membership_file, "\n")
cat(summary_file, "\n")
cat(representative_file, "\n")

cat("\n============================================================\n")
cat("SESSION 17 PREREGISTERED DOMAIN SUMMARY COMPLETE\n")
cat("============================================================\n")
