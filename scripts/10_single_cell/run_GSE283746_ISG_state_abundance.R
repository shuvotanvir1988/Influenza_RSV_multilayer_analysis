suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(broom)
})

infile <- paste0(
    "results/scrnaseq_validation/GSE283746/isg_state/design/",
    "GSE283746_ISG_state_participant_counts_FROZEN_v1.0.tsv"
)

outfile <- paste0(
    "results/scrnaseq_validation/GSE283746/isg_state/",
    "GSE283746_ISG_state_abundance_primary.tsv"
)

d <- read_tsv(infile, show_col_types = FALSE)

# Recover participant-level age/sex from the frozen H5AD-derived
# pseudobulk metadata already generated for the primary analysis.
meta_file <- paste0(
    "results/scrnaseq_validation/GSE283746/pseudobulk/",
    "GSE283746_primary_pseudobulk_sample_metadata_FROZEN_v1.0.tsv"
)

if (!file.exists(meta_file)) {
    stop(
        "Expected pseudobulk metadata file not found: ",
        meta_file,
        "\nLocate the existing participant-level pseudobulk metadata ",
        "before running the abundance model. Do not reconstruct ",
        "age/sex manually."
    )
}

meta <- read_tsv(meta_file, show_col_types = FALSE)

cat("=== INPUT ===\n")
cat("Rows:", nrow(d), "\n")
cat("Pairs:", n_distinct(d$pair_id), "\n")

cat("\n=== METADATA COLUMNS ===\n")
print(names(meta))

# Identify one row per participant.
# Expected participant identifier is Sample.
meta <- meta %>%
    distinct(Sample, .keep_all = TRUE)

required <- c("Sample", "Age", "Sex")
missing <- setdiff(required, names(meta))

if (length(missing) > 0) {
    stop(
        "Missing required metadata columns: ",
        paste(missing, collapse = ", ")
    )
}

meta <- meta %>%
    mutate(
        Age = as.numeric(Age),
        Sex = factor(Sex)
    )

age_mean <- mean(meta$Age, na.rm = TRUE)
age_sd   <- sd(meta$Age, na.rm = TRUE)

meta <- meta %>%
    mutate(
        age_z = (Age - age_mean) / age_sd
    )

x <- d %>%
    left_join(
        meta %>% select(Sample, Age, age_z, Sex),
        by = "Sample"
    ) %>%
    mutate(
        RSV_status = factor(
            `Combined Condition`,
            levels = c("Healthy", "RSV")
        )
    )

results <- list()

for (p in unique(x$pair_id)) {

    z <- x %>% filter(pair_id == p)

    # Participants with no cells in either member of the pair
    # cannot contribute to a binomial denominator.
    z <- z %>%
        filter(
            ISG_cells + Reference_cells > 0,
            !is.na(age_z),
            !is.na(Sex)
        )

    cat("\n========================================\n")
    cat(p, "\n")
    cat("========================================\n")

    print(
        z %>%
            group_by(RSV_status) %>%
            summarise(
                participants = n(),
                ISG_cells = sum(ISG_cells),
                reference_cells = sum(Reference_cells),
                median_ISG_fraction =
                    median(ISG_fraction, na.rm = TRUE),
                mean_ISG_fraction =
                    mean(ISG_fraction, na.rm = TRUE),
                .groups = "drop"
            )
    )

    fit <- glm(
        cbind(ISG_cells, Reference_cells) ~
            age_z + Sex + RSV_status,
        family = quasibinomial(),
        data = z
    )

    tt <- broom::tidy(fit)

    r <- tt %>%
        filter(term == "RSV_statusRSV") %>%
        mutate(
            pair_id = p,
            participants = nrow(z),
            odds_ratio = exp(estimate),
            CI_low = exp(estimate - 1.96 * std.error),
            CI_high = exp(estimate + 1.96 * std.error)
        )

    results[[p]] <- r
}

res <- bind_rows(results) %>%
    mutate(
        FDR = p.adjust(p.value, method = "BH")
    ) %>%
    select(
        pair_id,
        participants,
        estimate,
        std.error,
        statistic,
        p.value,
        FDR,
        odds_ratio,
        CI_low,
        CI_high
    ) %>%
    arrange(FDR)

write_tsv(res, outfile)

cat("\n========================================\n")
cat("FINAL ABUNDANCE RESULTS\n")
cat("========================================\n")

print(res, n = Inf)

cat("\nWritten:\n")
cat(outfile, "\n")
