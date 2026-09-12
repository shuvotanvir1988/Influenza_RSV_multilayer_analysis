#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

root <- normalizePath(Sys.getenv("INFLUENZA_RSV_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

meta_file <- file.path(
    root,
    "data/processed/DS001_GSE38900/DS001_GSE38900_master_sample_metadata_harmonized.tsv"
)

outdir <- file.path(
    root,
    "results/DS001_GSE38900/cell_composition/tables"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(
    meta_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

to_bool <- function(x) {
    tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

meta$primary_eligible <-
    to_bool(meta$eligible_primary_infection_vs_control)

meta$direct_eligible <-
    to_bool(meta$eligible_direct_influenza_vs_rsv)

meta$rsv_eligible <-
    to_bool(meta$eligible_rsv_acute_vs_control)

meta$age_model_eligible <-
    to_bool(meta$eligible_age_adjusted_model)

keep_cols <- c(
    "canonical_sample_id",
    "platform_id",
    "geo_accession",
    "description",
    "harmonized_group",
    "pathogen",
    "infection_phase",
    "sex_standardized",
    "age_months_harmonized",
    "primary_eligible",
    "direct_eligible",
    "rsv_eligible",
    "age_model_eligible"
)

inventory <- meta[, keep_cols]

write.table(
    inventory,
    file.path(
        outdir,
        "DS001_GSE38900_cell_composition_sample_inventory.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("=== PLATFORM × HARMONIZED GROUP ===\n")
print(
    with(
        inventory,
        table(platform_id, harmonized_group)
    )
)

cat("\n=== PRIMARY-ELIGIBLE ===\n")
print(
    with(
        inventory[inventory$primary_eligible, ],
        table(platform_id, harmonized_group)
    )
)

cat("\n=== DIRECT INFLUENZA-vs-RSV ELIGIBLE ===\n")
print(
    with(
        inventory[inventory$direct_eligible, ],
        table(platform_id, harmonized_group)
    )
)

cat("\n=== RSV ACUTE-vs-CONTROL ELIGIBLE ===\n")
print(
    with(
        inventory[inventory$rsv_eligible, ],
        table(platform_id, harmonized_group)
    )
)

cat("\n=== INFECTION PHASE ===\n")
print(
    with(
        inventory,
        table(platform_id, infection_phase)
    )
)

cat("\n=== AGE-ADJUSTED MODEL ELIGIBILITY ===\n")
print(
    with(
        inventory,
        table(platform_id, age_model_eligible)
    )
)
