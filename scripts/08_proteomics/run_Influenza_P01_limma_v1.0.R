suppressPackageStartupMessages({
    library(limma)
})

ROOT <- "results/proteomics_validation/Influenza-P01"

expr_file <- file.path(
    ROOT,
    "analysis_inputs",
    "Influenza-P01_human_expression_matrix_v1.0.tsv"
)

ann_file <- file.path(
    ROOT,
    "analysis_inputs",
    "Influenza-P01_human_assay_annotation_v1.0.tsv"
)

meta_file <- file.path(
    ROOT,
    "analysis_inputs",
    "Influenza-P01_sample_metadata_v1.0.tsv"
)

outdir <- file.path(
    ROOT,
    "differential_proteomics"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

expr <- read.delim(
    expr_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

ann <- read.delim(
    ann_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

meta <- read.delim(
    meta_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

rownames(expr) <- as.character(expr$row_index)
expr$row_index <- NULL

stopifnot(
    identical(
        colnames(expr),
        meta$sample_ID
    )
)

Y <- as.matrix(expr)
storage.mode(Y) <- "double"

meta$infect_status <- factor(
    meta$infect_status,
    levels=c(
        "healthy_control",
        "infected"
    )
)

meta$sex <- factor(
    meta$sex,
    levels=c(
        "f",
        "m"
    )
)

meta$ICU <- factor(
    meta$ICU,
    levels=c(
        "n",
        "y"
    )
)

# ============================================================
# PRIMARY ADJUSTED MODEL
# protein ~ infection + sex
# ============================================================

design <- model.matrix(
    ~ infect_status + sex,
    data=meta
)

cat("=== PRIMARY DESIGN MATRIX ===\n")
print(colnames(design))
print(table(meta$infect_status, meta$sex))

fit0 <- lmFit(
    Y,
    design
)

coef_name <- "infect_statusinfected"

raw_sigma <- fit0$sigma

fit <- eBayes(
    fit0
)

tt <- topTable(
    fit,
    coef=coef_name,
    number=Inf,
    sort.by="none"
)

tt$row_index <- rownames(tt)

tt$residual_sigma <- raw_sigma[
    match(
        tt$row_index,
        rownames(Y)
    )
]

tt$standardized_effect <- (
    tt$logFC /
    tt$residual_sigma
)

primary <- merge(
    ann,
    tt,
    by="row_index",
    sort=FALSE
)

write.table(
    primary,
    file=file.path(
        outdir,
        "Influenza-P01_primary_adjusted_infection_limma_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# UNADJUSTED SENSITIVITY
# protein ~ infection
# ============================================================

design_u <- model.matrix(
    ~ infect_status,
    data=meta
)

fit_u0 <- lmFit(
    Y,
    design_u
)

sigma_u <- fit_u0$sigma

fit_u <- eBayes(
    fit_u0
)

tt_u <- topTable(
    fit_u,
    coef="infect_statusinfected",
    number=Inf,
    sort.by="none"
)

tt_u$row_index <- rownames(tt_u)

tt_u$residual_sigma <- sigma_u[
    match(
        tt_u$row_index,
        rownames(Y)
    )
]

tt_u$standardized_effect <- (
    tt_u$logFC /
    tt_u$residual_sigma
)

unadjusted <- merge(
    ann,
    tt_u,
    by="row_index",
    sort=FALSE
)

write.table(
    unadjusted,
    file=file.path(
        outdir,
        "Influenza-P01_unadjusted_infection_limma_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# NON-ICU SENSITIVITY
# controls + infected ICU=n
# ============================================================

keep_nonICU <- (
    meta$infect_status == "healthy_control" |
    (
        meta$infect_status == "infected" &
        meta$ICU == "n"
    )
)

meta_n <- droplevels(
    meta[
        keep_nonICU,
    ]
)

Y_n <- Y[
    ,
    keep_nonICU,
    drop=FALSE
]

stopifnot(
    sum(
        meta_n$infect_status ==
        "healthy_control"
    ) == 23
)

stopifnot(
    sum(
        meta_n$infect_status ==
        "infected"
    ) == 35
)

design_n <- model.matrix(
    ~ infect_status + sex,
    data=meta_n
)

fit_n0 <- lmFit(
    Y_n,
    design_n
)

sigma_n <- fit_n0$sigma

fit_n <- eBayes(
    fit_n0
)

tt_n <- topTable(
    fit_n,
    coef="infect_statusinfected",
    number=Inf,
    sort.by="none"
)

tt_n$row_index <- rownames(tt_n)

tt_n$residual_sigma <- sigma_n[
    match(
        tt_n$row_index,
        rownames(Y_n)
    )
]

tt_n$standardized_effect <- (
    tt_n$logFC /
    tt_n$residual_sigma
)

nonICU <- merge(
    ann,
    tt_n,
    by="row_index",
    sort=FALSE
)

write.table(
    nonICU,
    file=file.path(
        outdir,
        "Influenza-P01_nonICU_infection_limma_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# SECONDARY ICU SEVERITY MODEL
# infected only: protein ~ ICU + sex
# ============================================================

keep_inf <- (
    meta$infect_status ==
    "infected"
)

meta_i <- droplevels(
    meta[
        keep_inf,
    ]
)

Y_i <- Y[
    ,
    keep_inf,
    drop=FALSE
]

stopifnot(
    sum(
        meta_i$ICU ==
        "n"
    ) == 35
)

stopifnot(
    sum(
        meta_i$ICU ==
        "y"
    ) == 26
)

design_i <- model.matrix(
    ~ ICU + sex,
    data=meta_i
)

fit_i0 <- lmFit(
    Y_i,
    design_i
)

sigma_i <- fit_i0$sigma

fit_i <- eBayes(
    fit_i0
)

tt_i <- topTable(
    fit_i,
    coef="ICUy",
    number=Inf,
    sort.by="none"
)

tt_i$row_index <- rownames(tt_i)

tt_i$residual_sigma <- sigma_i[
    match(
        tt_i$row_index,
        rownames(Y_i)
    )
]

tt_i$standardized_effect <- (
    tt_i$logFC /
    tt_i$residual_sigma
)

icu <- merge(
    ann,
    tt_i,
    by="row_index",
    sort=FALSE
)

write.table(
    icu,
    file=file.path(
        outdir,
        "Influenza-P01_secondary_ICU_severity_limma_v1.0.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# SUMMARY
# ============================================================

cat("\n=== PRIMARY MODEL ===\n")
cat(
    "Assays:",
    nrow(primary),
    "\n"
)

cat(
    "FDR < 0.05:",
    sum(
        primary$adj.P.Val < 0.05
    ),
    "\n"
)

cat(
    "Positive FDR < 0.05:",
    sum(
        primary$adj.P.Val < 0.05 &
        primary$logFC > 0
    ),
    "\n"
)

cat(
    "Negative FDR < 0.05:",
    sum(
        primary$adj.P.Val < 0.05 &
        primary$logFC < 0
    ),
    "\n"
)

cat(
    "Nominal P < 0.05:",
    sum(
        primary$P.Value < 0.05
    ),
    "\n"
)

cat("\n=== UNADJUSTED MODEL ===\n")
cat(
    "FDR < 0.05:",
    sum(
        unadjusted$adj.P.Val < 0.05
    ),
    "\n"
)

cat("\n=== NON-ICU SENSITIVITY ===\n")

cat(
    "Controls:",
    sum(
        meta_n$infect_status ==
        "healthy_control"
    ),
    "\n"
)

cat(
    "Infected:",
    sum(
        meta_n$infect_status ==
        "infected"
    ),
    "\n"
)

cat(
    "FDR < 0.05:",
    sum(
        nonICU$adj.P.Val < 0.05
    ),
    "\n"
)

cat("\n=== ICU SEVERITY ===\n")

cat(
    "Non-ICU:",
    sum(
        meta_i$ICU == "n"
    ),
    "\n"
)

cat(
    "ICU:",
    sum(
        meta_i$ICU == "y"
    ),
    "\n"
)

cat(
    "FDR < 0.05:",
    sum(
        icu$adj.P.Val < 0.05
    ),
    "\n"
)

cat("\nSTATUS: LIMMA COMPLETE\n")
