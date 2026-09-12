#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

infile <- paste0(
    "results/DS001_GSE38900/cell_composition/tables/",
    "DS001_GSE38900_MCPcounter_disease_comparisons_final.tsv"
)

outfile <- paste0(
    "results/DS001_GSE38900/cell_composition/tables/",
    "DS001_GSE38900_RSV_MCPcounter_crossplatform_replication.tsv"
)

x <- read.delim(
    infile,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

a <- x[
    x$contrast == "GPL10558_RSV_vs_control" &
    x$primary_immune_population,
    c(
        "cell_type",
        "hedges_g",
        "hedges_g_ci_low",
        "hedges_g_ci_high",
        "model_p",
        "model_FDR_primary8"
    )
]

b <- x[
    x$contrast == "GPL6884_RSV_vs_control" &
    x$primary_immune_population,
    c(
        "cell_type",
        "hedges_g",
        "hedges_g_ci_low",
        "hedges_g_ci_high",
        "model_p",
        "model_FDR_primary8"
    )
]

colnames(a)[-1] <- paste0(
    colnames(a)[-1],
    "_GPL10558"
)

colnames(b)[-1] <- paste0(
    colnames(b)[-1],
    "_GPL6884"
)

m <- merge(
    a,
    b,
    by = "cell_type",
    all = FALSE
)

m$direction_GPL10558 <- ifelse(
    m$hedges_g_GPL10558 > 0,
    "Higher_in_RSV",
    "Lower_in_RSV"
)

m$direction_GPL6884 <- ifelse(
    m$hedges_g_GPL6884 > 0,
    "Higher_in_RSV",
    "Lower_in_RSV"
)

m$direction_concordant <-
    m$direction_GPL10558 ==
    m$direction_GPL6884

m$significant_GPL10558 <-
    m$model_FDR_primary8_GPL10558 < 0.05

m$significant_GPL6884 <-
    m$model_FDR_primary8_GPL6884 < 0.05

m$replication_class <- ifelse(
    m$significant_GPL10558 &
    m$significant_GPL6884 &
    m$direction_concordant,
    "Significant_both",
    ifelse(
        m$direction_concordant &
        (m$significant_GPL10558 |
         m$significant_GPL6884),
        "Directionally_replicated_one_significant",
        ifelse(
            m$direction_concordant,
            "Directionally_concordant",
            "Directionally_discordant"
        )
    )
)

write.table(
    m,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("=== RSV CROSS-PLATFORM COMPOSITION REPLICATION ===\n\n")

print(
    m[, c(
        "cell_type",
        "hedges_g_GPL10558",
        "model_FDR_primary8_GPL10558",
        "hedges_g_GPL6884",
        "model_FDR_primary8_GPL6884",
        "direction_concordant",
        "replication_class"
    )],
    row.names = FALSE
)

cat("\n=== EFFECT-SIZE CORRELATION ===\n")

cat(
    "Pearson r = ",
    round(
        cor(
            m$hedges_g_GPL10558,
            m$hedges_g_GPL6884,
            method = "pearson"
        ),
        3
    ),
    "\n",
    sep = ""
)

cat(
    "Spearman rho = ",
    round(
        cor(
            m$hedges_g_GPL10558,
            m$hedges_g_GPL6884,
            method = "spearman"
        ),
        3
    ),
    "\n",
    sep = ""
)
