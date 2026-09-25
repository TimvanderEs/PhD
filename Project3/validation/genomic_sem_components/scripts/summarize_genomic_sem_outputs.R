#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: summarize_genomic_sem_outputs.R <10k|22k> <result_dir>")
}

freeze <- args[[1L]]
result_dir <- normalizePath(args[[2L]], mustWork = TRUE)
trait_ids <- c(
  "reaction_time",
  "stroop_box",
  "stroop_ink",
  "quiz_score",
  "working_memory",
  "pairing_guesses"
)

ldsc <- readRDS(file.path(result_dir, "ldsc_output.rds"))
if (!identical(dim(ldsc$S), c(length(trait_ids), length(trait_ids)))) {
  stop("Unexpected genetic-covariance matrix dimensions")
}

dimnames(ldsc$S) <- list(trait_ids, trait_ids)
if (!is.null(ldsc$I) && identical(dim(ldsc$I), dim(ldsc$S))) {
  dimnames(ldsc$I) <- list(trait_ids, trait_ids)
}
saveRDS(ldsc, file.path(result_dir, "ldsc_output.rds"))
write.table(
  ldsc$S,
  file.path(result_dir, "genetic_covariance.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  stats::cov2cor(ldsc$S),
  file.path(result_dir, "genetic_correlation.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
if (!is.null(ldsc$I)) {
  write.table(
    ldsc$I,
    file.path(result_dir, "intercept_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

# GenomicSEM vectorises the lower triangle column by column. The genetic-
# variance entries therefore occupy 1, 7, 12, 16, 19 and 21 for six traits.
variance_indices <- c(1L, 1L + cumsum(length(trait_ids):2L))
h2 <- diag(ldsc$S)
h2_se <- sqrt(diag(ldsc$V)[variance_indices])
h2_z <- h2 / h2_se
h2_table <- data.frame(
  freeze = freeze,
  trait = trait_ids,
  h2 = h2,
  se = h2_se,
  z = h2_z,
  p = 2 * pnorm(abs(h2_z), lower.tail = FALSE),
  stringsAsFactors = FALSE
)
write.table(
  h2_table,
  file.path(result_dir, "heritability_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

fits <- read.delim(file.path(result_dir, "model_comparison.tsv"), check.names = FALSE)
one <- fits[fits$model == "one_factor", , drop = FALSE]
two <- fits[fits$model == "two_factor", , drop = FALSE]
if (nrow(one) != 1L || nrow(two) != 1L) {
  stop("Expected one row for each structural model")
}

delta_chisq <- one$chisq - two$chisq
delta_df <- one$df - two$df
contrast <- data.frame(
  freeze = freeze,
  comparison = "one_factor_vs_two_factor",
  delta_chisq = delta_chisq,
  delta_df = delta_df,
  p_chisq_difference_approx = if (delta_chisq >= 0 && delta_df > 0) {
    pchisq(delta_chisq, df = delta_df, lower.tail = FALSE)
  } else {
    NA_real_
  },
  delta_AIC_two_minus_one = two$AIC - one$AIC,
  delta_SRMR_two_minus_one = two$SRMR - one$SRMR,
  stringsAsFactors = FALSE
)
write.table(
  contrast,
  file.path(result_dir, "model_contrast.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

read_nonempty <- function(path) {
  file.exists(path) && file.info(path)$size > 0L
}

one_parameters <- read.delim(
  file.path(result_dir, "one_factor_parameters.tsv"),
  check.names = FALSE
)
two_parameters <- read.delim(
  file.path(result_dir, "two_factor_parameters.tsv"),
  check.names = FALSE
)

residual_is_negative <- function(parameters) {
  any(
    parameters$op == "~~" &
      parameters$lhs == parameters$rhs &
      is.finite(parameters$Unstand_Est) &
      parameters$Unstand_Est < 0
  )
}

eigenvalues <- eigen(ldsc$S, symmetric = TRUE, only.values = TRUE)$values
diagnostics <- data.frame(
  freeze = freeze,
  minimum_genetic_covariance_eigenvalue = min(eigenvalues),
  covariance_positive_semidefinite = min(eigenvalues) >= -sqrt(.Machine$double.eps),
  one_factor_warning_file_nonempty = read_nonempty(file.path(result_dir, "one_factor_warnings.txt")),
  two_factor_warning_file_nonempty = read_nonempty(file.path(result_dir, "two_factor_warnings.txt")),
  one_factor_negative_residual_variance = residual_is_negative(one_parameters),
  two_factor_negative_residual_variance = residual_is_negative(two_parameters),
  stringsAsFactors = FALSE
)
write.table(
  diagnostics,
  file.path(result_dir, "model_diagnostics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Summarized", freeze, "results in", result_dir, "\n")
