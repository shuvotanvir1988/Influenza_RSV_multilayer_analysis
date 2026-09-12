#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

hallmark_file <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/hallmark_comparison",
  "DS001_GSE38900_Hallmark_contrast_comparison.tsv"
)

reactome_file <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/reactome_comparison",
  "DS001_GSE38900_Reactome_contrast_comparison.tsv.gz"
)

domain_file <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/tables/preregistered_domains",
  "DS001_GSE38900_preregistered_domain_representative_pathways.tsv"
)

outdir <- file.path(
  root,
  "results/DS001_GSE38900/pathway_analysis/figures"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("============================================================\n")
cat("SESSION 17 — PATHWAY FIGURE GENERATION\n")
cat("============================================================\n")

# ============================================================
# Helper
# ============================================================

clean_name <- function(x) {

  x <- gsub(
    "^HALLMARK_",
    "",
    x
  )

  x <- gsub(
    "^REACTOME_",
    "",
    x
  )

  x <- gsub(
    "_",
    " ",
    x
  )

  x
}

# ============================================================
# FIGURE 1 — Hallmark NES heatmap
# ============================================================

h <- fread(
  hallmark_file
)

hallmark_long <- rbindlist(
  list(

    h[
      ,
      .(
        pathway,
        contrast = "RSV GPL10558",
        NES = NES_RSV_GPL10558,
        padj = padj_RSV_GPL10558
      )
    ],

    h[
      ,
      .(
        pathway,
        contrast = "RSV GPL6884",
        NES = NES_RSV_GPL6884,
        padj = padj_RSV_GPL6884
      )
    ],

    h[
      ,
      .(
        pathway,
        contrast = "Influenza A GPL6884",
        NES = NES_Influenza_GPL6884,
        padj = padj_Influenza_GPL6884
      )
    ],

    h[
      ,
      .(
        pathway,
        contrast = "Influenza A vs RSV",
        NES = NES_Influenza_vs_RSV_GPL6884,
        padj = padj_Influenza_vs_RSV_GPL6884
      )
    ]
  )
)

hallmark_long[
  ,
  label := clean_name(pathway)
]

hallmark_long[
  ,
  contrast := factor(
    contrast,
    levels = c(
      "RSV GPL10558",
      "RSV GPL6884",
      "Influenza A GPL6884",
      "Influenza A vs RSV"
    )
  )
]

# Order pathways by maximum absolute NES.
order_h <- hallmark_long[
  ,
  .(
    max_abs_NES = max(
      abs(NES),
      na.rm = TRUE
    )
  ),
  by = label
][
  order(max_abs_NES)
]$label

hallmark_long[
  ,
  label := factor(
    label,
    levels = order_h
  )
]

p1 <- ggplot(
  hallmark_long,
  aes(
    x = contrast,
    y = label,
    fill = NES
  )
) +
  geom_tile() +
  geom_point(
    data = hallmark_long[
      !is.na(padj) &
      padj < 0.05
    ],
    aes(
      x = contrast,
      y = label
    ),
    inherit.aes = FALSE,
    shape = 8,
    size = 1.7
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    name = "NES"
  ) +
  labs(
    title = "Hallmark pathway enrichment across influenza A and RSV contrasts",
    subtitle = "* indicates FDR < 0.05",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

ggsave(
  file.path(
    outdir,
    "Figure_S17_Hallmark_NES_heatmap.png"
  ),
  p1,
  width = 9,
  height = 12,
  dpi = 300
)

ggsave(
  file.path(
    outdir,
    "Figure_S17_Hallmark_NES_heatmap.pdf"
  ),
  p1,
  width = 9,
  height = 12
)

# ============================================================
# FIGURE 2 — Preregistered-domain representative pathways
# ============================================================

d <- fread(
  domain_file
)

# Keep up to five most pathogen-differential pathways per domain.
d[
  ,
  abs_direct_NES :=
    abs(
      NES_Influenza_vs_RSV_GPL6884
    )
]

selected <- d[
  order(
    preregistered_domain,
    padj_Influenza_vs_RSV_GPL6884,
    -abs_direct_NES
  ),
  head(.SD, 5),
  by = preregistered_domain
]

selected[
  ,
  display :=
    paste0(
      preregistered_domain,
      " | ",
      clean_name(pathway)
    )
]

domain_long <- rbindlist(
  list(

    selected[
      ,
      .(
        display,
        preregistered_domain,
        contrast = "RSV GPL10558",
        NES = NES_RSV_GPL10558
      )
    ],

    selected[
      ,
      .(
        display,
        preregistered_domain,
        contrast = "RSV GPL6884",
        NES = NES_RSV_GPL6884
      )
    ],

    selected[
      ,
      .(
        display,
        preregistered_domain,
        contrast = "Influenza A GPL6884",
        NES = NES_Influenza_GPL6884
      )
    ],

    selected[
      ,
      .(
        display,
        preregistered_domain,
        contrast = "Influenza A vs RSV",
        NES = NES_Influenza_vs_RSV_GPL6884
      )
    ]
  )
)

domain_long[
  ,
  contrast := factor(
    contrast,
    levels = c(
      "RSV GPL10558",
      "RSV GPL6884",
      "Influenza A GPL6884",
      "Influenza A vs RSV"
    )
  )
]

display_order <- unique(
  selected$display
)

domain_long[
  ,
  display := factor(
    display,
    levels = rev(display_order)
  )
]

p2 <- ggplot(
  domain_long,
  aes(
    x = contrast,
    y = display,
    fill = NES
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint = 0,
    name = "NES"
  ) +
  labs(
    title = "Representative pathways from preregistered biological domains",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

ggsave(
  file.path(
    outdir,
    "Figure_S17_Pregistered_domain_heatmap.png"
  ),
  p2,
  width = 11,
  height = 11,
  dpi = 300
)

ggsave(
  file.path(
    outdir,
    "Figure_S17_Pregistered_domain_heatmap.pdf"
  ),
  p2,
  width = 11,
  height = 11
)

# ============================================================
# FIGURE 3 — RSV cross-platform Reactome concordance
# ============================================================

r <- fread(
  reactome_file
)

r[
  ,
  replication_status :=
    fifelse(
      RSV_replicated,
      "Replicated FDR<0.05",
      "Other"
    )
]

rho <- cor(
  r$NES_RSV_GPL10558,
  r$NES_RSV_GPL6884,
  method = "spearman",
  use = "complete.obs"
)

p3 <- ggplot(
  r,
  aes(
    x = NES_RSV_GPL10558,
    y = NES_RSV_GPL6884
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(
      shape = replication_status
    ),
    alpha = 0.6,
    size = 1.6
  ) +
  labs(
    title = "Cross-platform replication of RSV Reactome pathway enrichment",
    subtitle = paste0(
      "Spearman rho = ",
      round(rho, 3)
    ),
    x = "RSV NES — GPL10558",
    y = "RSV NES — GPL6884",
    shape = NULL
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(
    outdir,
    "Figure_S17_RSV_Reactome_crossplatform_concordance.png"
  ),
  p3,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    outdir,
    "Figure_S17_RSV_Reactome_crossplatform_concordance.pdf"
  ),
  p3,
  width = 7,
  height = 6
)

cat("\nSaved figures:\n")

cat(
  file.path(
    outdir,
    "Figure_S17_Hallmark_NES_heatmap.png"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "Figure_S17_Pregistered_domain_heatmap.png"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "Figure_S17_RSV_Reactome_crossplatform_concordance.png"
  ),
  "\n"
)

cat("\n============================================================\n")
cat("SESSION 17 PATHWAY FIGURE GENERATION COMPLETE\n")
cat("============================================================\n")
