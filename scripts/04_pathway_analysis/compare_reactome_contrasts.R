#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# SESSION 17 — REACTOME CROSS-CONTRAST COMPARISON
# AND BIOLOGICAL DOMAIN CLASSIFICATION
# ============================================================

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
  "results/DS001_GSE38900/pathway_analysis/tables/fgsea_reactome"
)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/reactome_comparison"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

FDR_THRESHOLD <- 0.05

cat("============================================================\n")
cat("SESSION 17 — REACTOME CONTRAST COMPARISON\n")
cat("============================================================\n")

cat("\nFrozen interpretation settings:\n")
cat("FDR threshold:", FDR_THRESHOLD, "\n")
cat("Direction based on NES sign\n")
cat("Negative NES in Influenza-vs-RSV = RSV-relative shift\n")
cat("Positive NES in Influenza-vs-RSV = Influenza-relative shift\n")

# ============================================================
# Input files
# ============================================================

files <- c(

  RSV_GPL10558 =
    "GPL10558_RSV_acute_vs_control_Reactome_fgsea.tsv",

  RSV_GPL6884 =
    "GPL6884_RSV_acute_vs_control_Reactome_fgsea.tsv",

  Influenza_GPL6884 =
    "GPL6884_InfluenzaA_vs_control_Reactome_fgsea.tsv",

  Influenza_vs_RSV_GPL6884 =
    "GPL6884_InfluenzaA_vs_RSVacute_Reactome_fgsea.tsv"
)

tables <- list()

# ============================================================
# Load and standardize each Reactome result
# ============================================================

for (nm in names(files)) {

  path <- file.path(
    indir,
    files[[nm]]
  )

  if (!file.exists(path)) {
    stop(
      "Missing Reactome result file:\n",
      path
    )
  }

  x <- fread(path)

  required <- c(
    "pathway",
    "size",
    "NES",
    "pval",
    "padj"
  )

  missing <- setdiff(
    required,
    names(x)
  )

  if (length(missing) > 0) {
    stop(
      nm,
      ": missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  x <- x[
    ,
    .(
      pathway,
      size,
      NES,
      pval,
      padj
    )
  ]

  setnames(
    x,
    c(
      "size",
      "NES",
      "pval",
      "padj"
    ),
    paste0(
      c(
        "size_",
        "NES_",
        "pval_",
        "padj_"
      ),
      nm
    )
  )

  tables[[nm]] <- x
}

# ============================================================
# Merge all contrasts
# ============================================================

merged <- Reduce(
  function(x, y) {

    merge(
      x,
      y,
      by = "pathway",
      all = TRUE
    )

  },
  tables
)

cat(
  "\nMerged Reactome pathways:",
  nrow(merged),
  "\n"
)

# ============================================================
# Significance flags
# ============================================================

merged[
  ,
  sig_RSV_GPL10558 :=
    !is.na(padj_RSV_GPL10558) &
    padj_RSV_GPL10558 < FDR_THRESHOLD
]

merged[
  ,
  sig_RSV_GPL6884 :=
    !is.na(padj_RSV_GPL6884) &
    padj_RSV_GPL6884 < FDR_THRESHOLD
]

merged[
  ,
  sig_Influenza_GPL6884 :=
    !is.na(padj_Influenza_GPL6884) &
    padj_Influenza_GPL6884 < FDR_THRESHOLD
]

merged[
  ,
  sig_Influenza_vs_RSV_GPL6884 :=
    !is.na(padj_Influenza_vs_RSV_GPL6884) &
    padj_Influenza_vs_RSV_GPL6884 < FDR_THRESHOLD
]

# ============================================================
# RSV cross-platform replication
# ============================================================

merged[
  ,
  RSV_same_direction :=
    !is.na(NES_RSV_GPL10558) &
    !is.na(NES_RSV_GPL6884) &
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

# ============================================================
# Shared infection response
#
# Compare influenza vs control and RSV vs control on GPL6884.
# ============================================================

merged[
  ,
  shared_same_direction :=
    !is.na(NES_RSV_GPL6884) &
    !is.na(NES_Influenza_GPL6884) &
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

# ============================================================
# Pathogen-differential response
# ============================================================

merged[
  ,
  pathogen_differential :=
    sig_Influenza_vs_RSV_GPL6884
]

merged[
  ,
  pathogen_shift :=
    fifelse(
      pathogen_differential &
      NES_Influenza_vs_RSV_GPL6884 > 0,
      "Influenza_shifted",

      fifelse(
        pathogen_differential &
        NES_Influenza_vs_RSV_GPL6884 < 0,
        "RSV_shifted",
        "Not_significant"
      )
    )
]

# ============================================================
# Biological-domain classification
#
# These categories correspond to the preregistered biological
# interpretation framework plus mechanistic domains emerging
# from the Reactome analysis.
#
# Classification is keyword-based and does NOT create new
# enrichment tests.
# ============================================================

classify_domain <- function(pathway) {

  p <- toupper(pathway)

  # ----------------------------------------------------------
  # Interferon / antiviral
  # ----------------------------------------------------------

  if (
    grepl(
      "INTERFERON|ANTIVIRAL|ISG15|OAS|RIG_I|MDA5|DDX58|IFIH1",
      p
    )
  ) {
    return("Interferon_antiviral")
  }

  # ----------------------------------------------------------
  # Cytokine / inflammatory signaling
  # ----------------------------------------------------------

  if (
    grepl(
      "INTERLEUKIN|CYTOKINE|TNF|NFKB|NF_KB|JAK|STAT|CHEMOKINE",
      p
    )
  ) {
    return("Cytokine_inflammation")
  }

  # ----------------------------------------------------------
  # Neutrophil / monocyte / myeloid
  # ----------------------------------------------------------

  if (
    grepl(
      "NEUTROPHIL|MONOCYTE|MYELOID|DEGRANULATION|GRANULOCYTE",
      p
    )
  ) {
    return("Myeloid_neutrophil")
  }

  # ----------------------------------------------------------
  # Antigen processing / adaptive immunity
  # ----------------------------------------------------------

  if (
    grepl(
      "ANTIGEN|MHC|TCR|T_CELL|B_CELL|IMMUNOGLOBULIN|LYMPHOCYTE",
      p
    )
  ) {
    return("Antigen_adaptive_immunity")
  }

  # ----------------------------------------------------------
  # Translation / ribosome / RNA quality control
  # ----------------------------------------------------------

  if (
    grepl(
      "TRANSLATION|RIBOSOM|RRNA|NONSENSE_MEDIATED|NMD|SRP_|EIF|PROTEIN_TARGETING",
      p
    )
  ) {
    return("Translation_ribosome")
  }

  # ----------------------------------------------------------
  # Cell cycle / proliferation
  # ----------------------------------------------------------

  if (
    grepl(
      "CELL_CYCLE|MITOTIC|MITOSIS|DNA_REPLICATION|CHROMATID|G1_|G2_|S_PHASE|CHECKPOINT",
      p
    )
  ) {
    return("Cell_cycle_proliferation")
  }

  # ----------------------------------------------------------
  # Mitochondrial / oxidative metabolism
  # ----------------------------------------------------------

  if (
    grepl(
      "MITOCHONDR|RESPIRATORY_ELECTRON|OXIDATIVE|TCA|TRICARBOXYLIC",
      p
    )
  ) {
    return("Mitochondrial_energy")
  }

  # ----------------------------------------------------------
  # Lipid / amino acid / metabolic programs
  # ----------------------------------------------------------

  if (
    grepl(
      "METABOLISM|FATTY_ACID|LIPID|CHOLESTEROL|AMINO_ACID|PEROXISOM|SELENOAMINO",
      p
    )
  ) {
    return("Metabolism")
  }

  # ----------------------------------------------------------
  # Stress / proteostasis
  # ----------------------------------------------------------

  if (
    grepl(
      "STARVATION|STRESS|UNFOLDED|AUTOPHAG|QUALITY_CONTROL|PROTEOSTASIS",
      p
    )
  ) {
    return("Stress_proteostasis")
  }

  # ----------------------------------------------------------
  # Apoptosis / cell death
  # ----------------------------------------------------------

  if (
    grepl(
      "APOPTOSIS|CELL_DEATH|CASPASE|NECRO",
      p
    )
  ) {
    return("Cell_death")
  }

  return("Other")
}

merged[
  ,
  biological_domain :=
    vapply(
      pathway,
      classify_domain,
      character(1)
    )
]

# ============================================================
# Correlations
# ============================================================

rsv_spearman <- cor(
  merged$NES_RSV_GPL10558,
  merged$NES_RSV_GPL6884,
  method = "spearman",
  use = "complete.obs"
)

infection_spearman <- cor(
  merged$NES_RSV_GPL6884,
  merged$NES_Influenza_GPL6884,
  method = "spearman",
  use = "complete.obs"
)

# ============================================================
# Global summary
# ============================================================

summary <- data.table(

  metric = c(

    "Reactome pathways represented",

    "RSV significant GPL10558",

    "RSV significant GPL6884",

    "RSV significant on both platforms",

    "RSV replicated same direction",

    "RSV NES Spearman correlation",

    "Shared RSV/influenza significant same direction",

    "RSV vs influenza NES Spearman correlation",

    "Significant influenza-vs-RSV pathways",

    "Influenza-shifted differential pathways",

    "RSV-shifted differential pathways"
  ),

  value = c(

    nrow(merged),

    sum(
      merged$sig_RSV_GPL10558,
      na.rm = TRUE
    ),

    sum(
      merged$sig_RSV_GPL6884,
      na.rm = TRUE
    ),

    sum(
      merged$RSV_significant_both,
      na.rm = TRUE
    ),

    sum(
      merged$RSV_replicated,
      na.rm = TRUE
    ),

    rsv_spearman,

    sum(
      merged$shared_infection_response,
      na.rm = TRUE
    ),

    infection_spearman,

    sum(
      merged$pathogen_differential,
      na.rm = TRUE
    ),

    sum(
      merged$pathogen_shift ==
        "Influenza_shifted",
      na.rm = TRUE
    ),

    sum(
      merged$pathogen_shift ==
        "RSV_shifted",
      na.rm = TRUE
    )
  )
)

# ============================================================
# Domain-level summary
# ============================================================

domains <- sort(
  unique(
    merged$biological_domain
  )
)

domain_rows <- list()

for (domain in domains) {

  x <- merged[
    biological_domain == domain
  ]

  domain_rows[[length(domain_rows) + 1]] <-
    data.table(

      biological_domain = domain,

      total_pathways =
        nrow(x),

      RSV_sig_GPL10558 =
        sum(
          x$sig_RSV_GPL10558,
          na.rm = TRUE
        ),

      RSV_sig_GPL6884 =
        sum(
          x$sig_RSV_GPL6884,
          na.rm = TRUE
        ),

      RSV_replicated =
        sum(
          x$RSV_replicated,
          na.rm = TRUE
        ),

      shared_infection_response =
        sum(
          x$shared_infection_response,
          na.rm = TRUE
        ),

      pathogen_differential =
        sum(
          x$pathogen_differential,
          na.rm = TRUE
        ),

      influenza_shifted =
        sum(
          x$pathogen_shift ==
            "Influenza_shifted",
          na.rm = TRUE
        ),

      RSV_shifted =
        sum(
          x$pathogen_shift ==
            "RSV_shifted",
          na.rm = TRUE
        )
    )
}

domain_summary <- rbindlist(
  domain_rows,
  fill = TRUE
)

setorder(
  domain_summary,
  -pathogen_differential,
  -RSV_replicated
)

# ============================================================
# Representative replicated RSV pathways
# ============================================================

replicated_rsv <- merged[
  RSV_replicated == TRUE
]

replicated_rsv[
  ,
  mean_abs_NES_RSV :=
    (
      abs(NES_RSV_GPL10558) +
      abs(NES_RSV_GPL6884)
    ) / 2
]

setorder(
  replicated_rsv,
  -mean_abs_NES_RSV
)

# ============================================================
# Representative shared infection pathways
# ============================================================

shared <- merged[
  shared_infection_response == TRUE
]

shared[
  ,
  mean_abs_NES_infection :=
    (
      abs(NES_RSV_GPL6884) +
      abs(NES_Influenza_GPL6884)
    ) / 2
]

setorder(
  shared,
  -mean_abs_NES_infection
)

# ============================================================
# Differential pathways
# ============================================================

differential <- merged[
  pathogen_differential == TRUE
]

differential[
  ,
  abs_direct_NES :=
    abs(
      NES_Influenza_vs_RSV_GPL6884
    )
]

setorder(
  differential,
  padj_Influenza_vs_RSV_GPL6884,
  -abs_direct_NES
)

# ============================================================
# Save complete comparison
# ============================================================

comparison_file <- file.path(
  outdir,
  "DS001_GSE38900_Reactome_contrast_comparison.tsv.gz"
)

summary_file <- file.path(
  outdir,
  "DS001_GSE38900_Reactome_comparison_summary.tsv"
)

domain_file <- file.path(
  outdir,
  "DS001_GSE38900_Reactome_domain_summary.tsv"
)

replication_file <- file.path(
  outdir,
  "DS001_GSE38900_Reactome_RSV_replicated.tsv"
)

shared_file <- file.path(
  outdir,
  "DS001_GSE38900_Reactome_shared_infection_response.tsv"
)

differential_file <- file.path(
  outdir,
  "DS001_GSE38900_Reactome_pathogen_differential.tsv"
)

fwrite(
  merged,
  comparison_file,
  sep = "\t"
)

fwrite(
  summary,
  summary_file,
  sep = "\t"
)

fwrite(
  domain_summary,
  domain_file,
  sep = "\t"
)

fwrite(
  replicated_rsv,
  replication_file,
  sep = "\t"
)

fwrite(
  shared,
  shared_file,
  sep = "\t"
)

fwrite(
  differential,
  differential_file,
  sep = "\t"
)

# ============================================================
# Terminal output
# ============================================================

cat("\n============================================================\n")
cat("REACTOME COMPARISON SUMMARY\n")
cat("============================================================\n")

print(summary)

cat("\n============================================================\n")
cat("BIOLOGICAL DOMAIN SUMMARY\n")
cat("============================================================\n")

print(domain_summary)

cat("\n============================================================\n")
cat("TOP 25 REPLICATED RSV PATHWAYS\n")
cat("============================================================\n")

print(
  replicated_rsv[
    1:min(25, .N),
    .(
      pathway,
      biological_domain,
      NES_RSV_GPL10558,
      NES_RSV_GPL6884,
      padj_RSV_GPL10558,
      padj_RSV_GPL6884
    )
  ]
)

cat("\n============================================================\n")
cat("TOP 25 SHARED INFLUENZA-RSV PATHWAYS\n")
cat("============================================================\n")

print(
  shared[
    1:min(25, .N),
    .(
      pathway,
      biological_domain,
      NES_RSV_GPL6884,
      NES_Influenza_GPL6884,
      padj_RSV_GPL6884,
      padj_Influenza_GPL6884
    )
  ]
)

cat("\n============================================================\n")
cat("TOP 30 PATHOGEN-DIFFERENTIAL PATHWAYS\n")
cat("============================================================\n")

print(
  differential[
    1:min(30, .N),
    .(
      pathway,
      biological_domain,
      NES_RSV_GPL6884,
      NES_Influenza_GPL6884,
      NES_Influenza_vs_RSV_GPL6884,
      padj_Influenza_vs_RSV_GPL6884,
      pathogen_shift
    )
  ]
)

cat("\nSaved:\n")
cat(comparison_file, "\n")
cat(summary_file, "\n")
cat(domain_file, "\n")
cat(replication_file, "\n")
cat(shared_file, "\n")
cat(differential_file, "\n")

cat("\n============================================================\n")
cat("SESSION 17 REACTOME COMPARISON COMPLETE\n")
cat("============================================================\n")
