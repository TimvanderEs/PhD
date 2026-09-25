#!/usr/bin/env Rscript

# Reconstruct the phenotype and covariate schemas required by the later
# refined Chinese-only HELIOS 10k REGENIE commands. The score input must
# already contain the historically selected column named g_raw; this script
# does not infer or rename a different cognitive-factor column.
#
# Usage:
#   Rscript prepare_regenie_10k_inputs.R \
#     cognitive_factor_scores.tsv sample_covariates.txt ancestry.eigenvec \
#     output_prefix

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: Rscript prepare_regenie_10k_inputs.R",
      "<score-file> <sample-covariates> <ancestry-eigenvec> <output-prefix>"
    ),
    call. = FALSE
  )
}

score_file <- args[[1L]]
sample_covariate_file <- args[[2L]]
eigenvector_file <- args[[3L]]
output_prefix <- args[[4L]]
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = ""
  )
}

first_column <- function(columns, candidates, label) {
  match <- candidates[candidates %in% columns]
  if (length(match) == 0L) {
    stop(
      sprintf("No %s column found; tried: %s", label, paste(candidates, collapse = ", ")),
      call. = FALSE
    )
  }
  match[[1L]]
}

scores <- read_table(score_file)
sample_covariates <- read_table(sample_covariate_file)
pcs <- read_table(eigenvector_file)

score_id <- first_column(names(scores), c("IID", "ID"), "score IID")
if (!"g_raw" %in% names(scores)) {
  stop(
    paste(
      "Score file must contain an explicit 'g_raw' column.",
      "The upstream score used by this historical branch was not retained,",
      "so another column must not be relabelled automatically."
    ),
    call. = FALSE
  )
}
covariate_id <- first_column(names(sample_covariates), c("IID", "ID"), "covariate IID")
age_column <- first_column(names(sample_covariates), c("Age", "FREG8_Age"), "age")
sex_column <- first_column(names(sample_covariates), c("Sex", "FREG7_Gender"), "sex")
pc_id <- first_column(names(pcs), c("IID", "ID"), "PC IID")

pc_columns <- paste0("PC", 1:20)
missing_pcs <- setdiff(pc_columns, names(pcs))
if (length(missing_pcs) > 0L) {
  stop(
    sprintf("Eigenvector file is missing: %s", paste(missing_pcs, collapse = ", ")),
    call. = FALSE
  )
}

scores$IID <- trimws(as.character(scores[[score_id]]))
scores$g_raw <- as.numeric(scores[["g_raw"]])
sample_covariates$IID <- trimws(as.character(sample_covariates[[covariate_id]]))
pcs$IID <- trimws(as.character(pcs[[pc_id]]))

scores <- scores[!duplicated(scores$IID) & nzchar(scores$IID), c("IID", "g_raw")]
sample_covariates <- sample_covariates[
  !duplicated(sample_covariates$IID) & nzchar(sample_covariates$IID),
  ,
  drop = FALSE
]
pcs <- pcs[!duplicated(pcs$IID) & nzchar(pcs$IID), , drop = FALSE]

covariate_rows <- match(scores$IID, sample_covariates$IID)
pc_rows <- match(scores$IID, pcs$IID)

combined <- data.frame(
  IID = scores$IID,
  g_raw = scores$g_raw,
  Age = as.numeric(sample_covariates[[age_column]][covariate_rows]),
  Sex = sample_covariates[[sex_column]][covariate_rows],
  stringsAsFactors = FALSE
)

for (column in pc_columns) {
  combined[[column]] <- as.numeric(pcs[[column]][pc_rows])
}

sex_numeric <- suppressWarnings(as.numeric(as.character(combined$Sex)))
if (any(!is.na(combined$Sex) & is.na(sex_numeric))) {
  stop("Sex must use the numeric coding applied in the historical analysis", call. = FALSE)
}

combined$Age2 <- combined$Age^2
combined$AgeSex <- combined$Age * sex_numeric
combined$Age2Sex <- combined$Age2 * sex_numeric

complete_columns <- c(
  "g_raw", "Age", "Sex", "Age2", "AgeSex", "Age2Sex", pc_columns
)
combined <- combined[complete.cases(combined[, complete_columns]), , drop = FALSE]
if (nrow(combined) == 0L) {
  stop("No participants have complete score, covariate and PC data", call. = FALSE)
}

pc_fid_column <- first_column(names(pcs), c("FID", "#FID", "IID"), "PC FID")
combined$FID <- pcs[[pc_fid_column]][match(combined$IID, pcs$IID)]
combined$FID[is.na(combined$FID) | !nzchar(as.character(combined$FID))] <-
  combined$IID[is.na(combined$FID) | !nzchar(as.character(combined$FID))]

phenotype_output <- combined[, c("FID", "IID", "g_raw")]
covariate_output <- combined[, c(
  "FID", "IID", "Age", "Sex", "Age2", "AgeSex", "Age2Sex", pc_columns
)]

write.table(
  phenotype_output,
  paste0(output_prefix, ".pheno"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
write.table(
  covariate_output,
  paste0(output_prefix, ".cov"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message(sprintf("Prepared REGENIE inputs for %d participants", nrow(combined)))
