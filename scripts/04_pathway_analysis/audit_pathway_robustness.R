#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

indir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/fgsea_reactome"
)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/robustness_audit"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("SESSION 17 — PATHWAY ROBUSTNESS AUDIT\n")
cat("============================================================\n")

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

res <- lapply(
  files,
  function(f) {
    fread(
      file.path(
        indir,
        f
      )
    )
  }
)

# ------------------------------------------------------------
# Theme definitions
# ------------------------------------------------------------

themes <- list(

  Interferon = c(
    "INTERFERON",
    "IFN_"
  ),

  Neutrophil_Myeloid = c(
    "NEUTROPHIL",
    "MYELOID",
    "DEGRANULATION"
  ),

  Antigen_Presentation = c(
    "ANTIGEN_PROCESSING",
    "CROSS_PRESENTATION",
    "MHC_CLASS",
    "PRESENTATION"
  ),

  Translation = c(
    "TRANSLATION",
    "RIBOSOM",
    "RRNA",
    "NONSENSE_MEDIATED",
    "SRP_DEPENDENT",
    "EIF2AK4"
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
  )
)

matches_theme <- function(pathway, patterns) {

  p <- toupper(pathway)

  any(
    vapply(
      patterns,
      function(pattern) {
        grepl(
          pattern,
          p,
          fixed = TRUE
        )
      },
      logical(1)
    )
  )
}

# ------------------------------------------------------------
# Extract theme-level enrichment
# ------------------------------------------------------------

theme_rows <- list()

for (contrast_name in names(res)) {

  x <- copy(
    res[[contrast_name]]
  )

  for (theme_name in names(themes)) {

    idx <- vapply(
      x$pathway,
      matches_theme,
      logical(1),
      patterns = themes[[theme_name]]
    )

    z <- copy(
      x[idx]
    )

    if (nrow(z) == 0) {
      next
    }

    z[, contrast := contrast_name]
    z[, theme := theme_name]

    theme_rows[[length(theme_rows) + 1]] <- z
  }
}

theme_table <- rbindlist(
  theme_rows,
  fill = TRUE
)

fwrite(
  theme_table,
  file.path(
    outdir,
    "DS001_GSE38900_major_theme_Reactome_pathways.tsv.gz"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# Leading-edge/core enrichment gene audit
# ------------------------------------------------------------

if (!"leadingEdge" %in% names(theme_table)) {

  stop(
    "leadingEdge column is absent from fgsea output. ",
    "Cannot perform leading-edge robustness audit."
  )
}

leading_rows <- list()

for (i in seq_len(nrow(theme_table))) {

  le <- theme_table$leadingEdge[[i]]

  if (length(le) == 0) {
    next
  }

  leading_rows[[length(leading_rows) + 1]] <-
    data.table(
      contrast = theme_table$contrast[i],
      theme = theme_table$theme[i],
      pathway = theme_table$pathway[i],
      NES = theme_table$NES[i],
      padj = theme_table$padj[i],
      leading_edge_gene = le
    )
}

leading <- rbindlist(
  leading_rows,
  fill = TRUE
)

fwrite(
  leading,
  file.path(
    outdir,
    "DS001_GSE38900_major_theme_leading_edge_genes.tsv.gz"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# Gene recurrence within each theme/contrast
# ------------------------------------------------------------

gene_recurrence <- leading[
  ,
  .(
    pathways_containing_gene =
      uniqueN(pathway),

    median_pathway_NES =
      median(
        NES,
        na.rm = TRUE
      ),

    minimum_pathway_FDR =
      min(
        padj,
        na.rm = TRUE
      )
  ),
  by = .(
    contrast,
    theme,
    leading_edge_gene
  )
]

setorder(
  gene_recurrence,
  contrast,
  theme,
  -pathways_containing_gene,
  minimum_pathway_FDR
)

fwrite(
  gene_recurrence,
  file.path(
    outdir,
    "DS001_GSE38900_leading_edge_gene_recurrence.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# RSV cross-platform pathway direction audit
# ------------------------------------------------------------

a <- res[["RSV_GPL10558"]][
  ,
  .(
    pathway,
    NES_GPL10558 = NES,
    padj_GPL10558 = padj
  )
]

b <- res[["RSV_GPL6884"]][
  ,
  .(
    pathway,
    NES_GPL6884 = NES,
    padj_GPL6884 = padj
  )
]

cross <- merge(
  a,
  b,
  by = "pathway",
  all = TRUE
)

cross[
  ,
  same_direction :=
    !is.na(NES_GPL10558) &
    !is.na(NES_GPL6884) &
    sign(NES_GPL10558) ==
    sign(NES_GPL6884)
]

cross[
  ,
  significant_both :=
    !is.na(padj_GPL10558) &
    !is.na(padj_GPL6884) &
    padj_GPL10558 < 0.05 &
    padj_GPL6884 < 0.05
]

cross[
  ,
  replicated :=
    significant_both &
    same_direction
]

cross[
  ,
  discordant_significant :=
    significant_both &
    !same_direction
]

cross[
  ,
  abs_NES_difference :=
    abs(
      NES_GPL10558 -
      NES_GPL6884
    )
]

setorder(
  cross,
  -discordant_significant,
  -abs_NES_difference
)

fwrite(
  cross,
  file.path(
    outdir,
    "DS001_GSE38900_RSV_crossplatform_pathway_direction_audit.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# Theme-level summary
# ------------------------------------------------------------

theme_summary <- theme_table[
  ,
  .(
    pathways = uniqueN(pathway),

    significant_FDR05 =
      uniqueN(
        pathway[
          !is.na(padj) &
          padj < 0.05
        ]
      ),

    positive_NES =
      uniqueN(
        pathway[
          !is.na(NES) &
          NES > 0
        ]
      ),

    negative_NES =
      uniqueN(
        pathway[
          !is.na(NES) &
          NES < 0
        ]
      ),

    leading_edge_genes =
      uniqueN(
        unlist(
          leadingEdge
        )
      )
  ),
  by = .(
    contrast,
    theme
  )
]

fwrite(
  theme_summary,
  file.path(
    outdir,
    "DS001_GSE38900_major_theme_robustness_summary.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# Terminal reporting
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("MAJOR THEME SUMMARY\n")
cat("============================================================\n")

print(theme_summary)

cat("\n============================================================\n")
cat("RSV CROSS-PLATFORM DIRECTION SUMMARY\n")
cat("============================================================\n")

cat(
  "Complete pathway pairs:",
  sum(
    !is.na(cross$NES_GPL10558) &
    !is.na(cross$NES_GPL6884)
  ),
  "\n"
)

cat(
  "Significant on both RSV platforms:",
  sum(
    cross$significant_both,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Replicated same direction:",
  sum(
    cross$replicated,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant on both but opposite direction:",
  sum(
    cross$discordant_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n============================================================\n")
cat("TOP SIGNIFICANT DISCORDANT RSV PATHWAYS\n")
cat("============================================================\n")

print(
  cross[
    discordant_significant == TRUE
  ][
    order(
      -abs_NES_difference
    )
  ][
    1:min(30, .N),
    .(
      pathway,
      NES_GPL10558,
      padj_GPL10558,
      NES_GPL6884,
      padj_GPL6884,
      abs_NES_difference
    )
  ]
)

cat("\n============================================================\n")
cat("TOP RECURRENT LEADING-EDGE GENES\n")
cat("============================================================\n")

for (contrast_name in unique(gene_recurrence$contrast)) {

  cat("\n### ", contrast_name, "\n", sep = "")

  for (theme_name in unique(
    gene_recurrence[
      contrast == contrast_name
    ]$theme
  )) {

    cat("\n", theme_name, "\n")

    print(
      gene_recurrence[
        contrast == contrast_name &
        theme == theme_name
      ][
        1:min(10, .N)
      ]
    )
  }
}

cat("\nSaved robustness-audit tables to:\n")
cat(outdir, "\n")

cat("\n============================================================\n")
cat("SESSION 17 PATHWAY ROBUSTNESS AUDIT COMPLETE\n")
cat("============================================================\n")
