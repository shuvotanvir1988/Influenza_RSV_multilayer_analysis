#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

infile <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/robustness_audit",
  "DS001_GSE38900_major_theme_Reactome_pathways.tsv.gz"
)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/robustness_audit"
)

cat("============================================================\n")
cat("SESSION 17 — CORRECTED LEADING-EDGE GENE AUDIT\n")
cat("============================================================\n")

x <- fread(infile)

if (!"leadingEdge" %in% names(x)) {
  stop("leadingEdge column not found.")
}

# ------------------------------------------------------------
# Parse semicolon-delimited leading-edge genes
# ------------------------------------------------------------

rows <- list()

for (i in seq_len(nrow(x))) {

  le <- x$leadingEdge[i]

  if (
    is.na(le) ||
    le == ""
  ) {
    next
  }

  genes <- unlist(
    strsplit(
      le,
      ";",
      fixed = TRUE
    )
  )

  genes <- trimws(genes)

  genes <- genes[
    !is.na(genes) &
    genes != ""
  ]

  genes <- unique(genes)

  if (length(genes) == 0) {
    next
  }

  rows[[length(rows) + 1]] <- data.table(

    contrast = x$contrast[i],

    theme = x$theme[i],

    pathway = x$pathway[i],

    NES = x$NES[i],

    padj = x$padj[i],

    leading_edge_gene = genes
  )
}

leading <- rbindlist(
  rows,
  fill = TRUE
)

# ------------------------------------------------------------
# True gene recurrence
# ------------------------------------------------------------

recurrence <- leading[
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
  recurrence,
  contrast,
  theme,
  -pathways_containing_gene,
  minimum_pathway_FDR,
  leading_edge_gene
)

# ------------------------------------------------------------
# Theme summary
# ------------------------------------------------------------

theme_gene_summary <- recurrence[
  ,
  .(
    unique_leading_edge_genes =
      uniqueN(leading_edge_gene),

    genes_in_at_least_2_pathways =
      sum(
        pathways_containing_gene >= 2
      ),

    genes_in_at_least_3_pathways =
      sum(
        pathways_containing_gene >= 3
      ),

    maximum_pathway_recurrence =
      max(
        pathways_containing_gene
      )
  ),
  by = .(
    contrast,
    theme
  )
]

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fwrite(
  leading,
  file.path(
    outdir,
    "DS001_GSE38900_corrected_leading_edge_genes.tsv.gz"
  ),
  sep = "\t"
)

fwrite(
  recurrence,
  file.path(
    outdir,
    "DS001_GSE38900_corrected_leading_edge_gene_recurrence.tsv"
  ),
  sep = "\t"
)

fwrite(
  theme_gene_summary,
  file.path(
    outdir,
    "DS001_GSE38900_corrected_leading_edge_theme_summary.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CORRECTED THEME GENE SUMMARY\n")
cat("============================================================\n")

print(theme_gene_summary)

cat("\n============================================================\n")
cat("TOP RECURRENT LEADING-EDGE GENES\n")
cat("============================================================\n")

themes_to_show <- c(
  "Interferon",
  "Neutrophil_Myeloid",
  "Antigen_Presentation",
  "Translation",
  "Metabolism"
)

for (contrast_name in unique(recurrence$contrast)) {

  cat(
    "\n############################################\n"
  )

  cat(
    "CONTRAST: ",
    contrast_name,
    "\n",
    sep = ""
  )

  cat(
    "############################################\n"
  )

  for (theme_name in themes_to_show) {

    z <- recurrence[
      contrast == contrast_name &
      theme == theme_name
    ]

    if (nrow(z) == 0) {
      next
    }

    cat(
      "\n### ",
      theme_name,
      "\n",
      sep = ""
    )

    print(
      z[
        1:min(15, .N),
        .(
          leading_edge_gene,
          pathways_containing_gene,
          median_pathway_NES,
          minimum_pathway_FDR
        )
      ]
    )
  }
}

cat("\n============================================================\n")
cat("CORRECTED LEADING-EDGE AUDIT COMPLETE\n")
cat("============================================================\n")
