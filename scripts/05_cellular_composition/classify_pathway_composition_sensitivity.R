#!/usr/bin/env Rscript

x <- read.delim(
    "results/DS001_GSE38900/cell_composition/sensitivity/DS001_GSE38900_composition_adjusted_pathway_effects.tsv",
    check.names = FALSE,
    stringsAsFactors = FALSE
)

classify <- function(base_fdr, adj_fdr, attenuation, sign_preserved) {

    if (!is.finite(base_fdr) || base_fdr >= 0.05) {
        return("Unclassified_base_not_significant")
    }

    if (!isTRUE(sign_preserved)) {
        return("Composition_sensitive_sign_reversal")
    }

    if (is.finite(attenuation) && attenuation < 0 && adj_fdr < 0.05) {
        return("Composition_robust_amplified")
    }

    if (
        is.finite(attenuation) &&
        attenuation >= 70
    ) {
        return("Composition_sensitive")
    }

    if (
        adj_fdr >= 0.05
    ) {
        return("Composition_sensitive")
    }

    if (
        is.finite(attenuation) &&
        attenuation >= 25
    ) {
        return("Mixed")
    }

    if (
        is.finite(attenuation) &&
        attenuation < 25 &&
        adj_fdr < 0.05
    ) {
        return("Composition_robust")
    }

    return("Unclassified")
}

x$composition_class <- mapply(
    classify,
    x$FDR_base,
    x$FDR_adjusted,
    x$attenuation_percent,
    x$sign_preserved
)

outfile <- paste0(
    "results/DS001_GSE38900/cell_composition/sensitivity/",
    "DS001_GSE38900_pathway_composition_classification.tsv"
)

write.table(
    x,
    outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

show <- x[, c(
    "contrast",
    "preregistered_domain",
    "beta_base",
    "FDR_base",
    "beta_adjusted",
    "FDR_adjusted",
    "attenuation_percent",
    "composition_class"
)]

print(show, row.names = FALSE)
