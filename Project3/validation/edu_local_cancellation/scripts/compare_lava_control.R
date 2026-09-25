#!/usr/bin/env Rscript

# Compare the isolated historical-control LAVA run with the archived output.
#
# LAVA 0.1.5 calculates bivariate confidence intervals and P values using
# unseeded Monte Carlo draws (ci.bivariate() and integral.p()). Exact equality
# is therefore expected for the row set, univariate estimates and point
# estimates, but not for the simulated interval bounds or P values. The gate
# below requires exact deterministic reproduction and close Monte Carlo
# agreement without any change in nominal or BH-FDR decisions.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: compare_lava_control.R archived.univ archived.bivar control.univ ",
    "control.bivar output.tsv"
  )
}

compare_table <- function(archived_path, control_path, keys, numeric_columns, table_name) {
  archived <- fread(archived_path, na.strings = c("NA", "NaN"))
  control <- fread(control_path, na.strings = c("NA", "NaN"))
  merged <- merge(
    archived,
    control,
    by = keys,
    all = TRUE,
    suffixes = c("_archived", "_control")
  )

  archived_marker <- paste0(setdiff(names(archived), keys)[1], "_archived")
  control_marker <- paste0(setdiff(names(control), keys)[1], "_control")
  archived_only <- sum(!is.na(merged[[archived_marker]]) & is.na(merged[[control_marker]]))
  control_only <- sum(is.na(merged[[archived_marker]]) & !is.na(merged[[control_marker]]))

  result <- data.table(
    table = table_name,
    metric = c("archived_rows", "control_rows", "archived_only_rows", "control_only_rows"),
    value = c(nrow(archived), nrow(control), archived_only, control_only)
  )

  for (column in numeric_columns) {
    first <- suppressWarnings(as.numeric(merged[[paste0(column, "_archived")]]))
    second <- suppressWarnings(as.numeric(merged[[paste0(column, "_control")]]))
    differences <- abs(first - second)
    maximum <- if (all(is.na(differences))) NA_real_ else max(differences, na.rm = TRUE)
    result <- rbind(
      result,
      data.table(
        table = table_name,
        metric = paste0("max_abs_difference_", column),
        value = maximum
      )
    )
  }
  result
}

monte_carlo_summary <- function(merged, label) {
  archived_p <- merged$p_archived
  control_p <- merged$p_control
  archived_fdr <- p.adjust(archived_p, method = "BH")
  control_fdr <- p.adjust(control_p, method = "BH")

  correlation <- function(first, second) {
    keep <- is.finite(first) & is.finite(second)
    if (sum(keep) < 2L) NA_real_ else cor(first[keep], second[keep])
  }
  maximum_difference <- function(first, second) {
    difference <- abs(first - second)
    if (all(is.na(difference))) NA_real_ else max(difference, na.rm = TRUE)
  }
  median_difference <- function(first, second) {
    difference <- abs(first - second)
    if (all(is.na(difference))) NA_real_ else median(difference, na.rm = TRUE)
  }

  metrics <- list(
    rows = nrow(merged),
    p_correlation = correlation(archived_p, control_p),
    p_median_abs_difference = median_difference(archived_p, control_p),
    p_max_abs_difference = maximum_difference(archived_p, control_p),
    rho_lower_correlation = correlation(merged$rho.lower_archived, merged$rho.lower_control),
    rho_upper_correlation = correlation(merged$rho.upper_archived, merged$rho.upper_control),
    r2_lower_correlation = correlation(merged$r2.lower_archived, merged$r2.lower_control),
    r2_upper_correlation = correlation(merged$r2.upper_archived, merged$r2.upper_control),
    nominal_significant_archived = sum(archived_p < 0.05, na.rm = TRUE),
    nominal_significant_control = sum(control_p < 0.05, na.rm = TRUE),
    nominal_decision_discordant = sum(
      (archived_p < 0.05) != (control_p < 0.05), na.rm = TRUE
    ),
    fdr_significant_archived = sum(archived_fdr < 0.05, na.rm = TRUE),
    fdr_significant_control = sum(control_fdr < 0.05, na.rm = TRUE),
    fdr_decision_discordant = sum(
      (archived_fdr < 0.05) != (control_fdr < 0.05), na.rm = TRUE
    )
  )
  data.table(
    table = paste0("bivariate_", label, "_monte_carlo"),
    metric = names(metrics),
    value = as.numeric(unlist(metrics, use.names = FALSE))
  )
}

comparison <- rbind(
  compare_table(
    args[1], args[3],
    c("locus", "phen"),
    c("chr", "start", "stop", "n.snps", "n.pcs", "h2.obs", "p"),
    "univariate"
  ),
  compare_table(
    args[2], args[4],
    c("locus", "phen1", "phen2"),
    c("chr", "start", "stop", "n.snps", "n.pcs", "rho", "r2"),
    "bivariate"
  )
)

archived_bivariate <- fread(args[2], na.strings = c("NA", "NaN"))
control_bivariate <- fread(args[4], na.strings = c("NA", "NaN"))
bivariate_merged <- merge(
  archived_bivariate,
  control_bivariate,
  by = c("locus", "phen1", "phen2"),
  all = FALSE,
  suffixes = c("_archived", "_control")
)
target_merged <- bivariate_merged[
  (phen1 == "HEL10k_EAS" & phen2 == "Edu_EAS") |
    (phen1 == "Edu_EAS" & phen2 == "HEL10k_EAS")
]
comparison <- rbind(
  comparison,
  monte_carlo_summary(bivariate_merged, "all"),
  monte_carlo_summary(target_merged, "helios_edu")
)

dir.create(dirname(args[5]), recursive = TRUE, showWarnings = FALSE)
fwrite(comparison, args[5], sep = "\t", na = "NA", quote = FALSE)
print(comparison)

row_mismatch <- comparison[
  metric %chin% c("archived_only_rows", "control_only_rows"),
  any(value != 0)
]
numeric_mismatch <- comparison[
  grepl("^max_abs_difference_", metric),
  any(!is.na(value) & value > 1e-10)
]
count_mismatch <- any(c(
  comparison[table == "univariate" & metric == "archived_rows", value] !=
    comparison[table == "univariate" & metric == "control_rows", value],
  comparison[table == "bivariate" & metric == "archived_rows", value] !=
    comparison[table == "bivariate" & metric == "control_rows", value]
))

monte_carlo_failure <- comparison[
  grepl("_monte_carlo$", table),
  any(c(
    value[metric == "p_correlation"] < 0.99,
    value[metric == "p_max_abs_difference"] > 0.01,
    value[metric == "rho_lower_correlation"] < 0.99,
    value[metric == "rho_upper_correlation"] < 0.99,
    value[metric == "r2_lower_correlation"] < 0.99,
    value[metric == "r2_upper_correlation"] < 0.99,
    value[metric == "nominal_decision_discordant"] != 0,
    value[metric == "fdr_decision_discordant"] != 0
  ), na.rm = TRUE),
  by = table
][, any(V1)]

if (row_mismatch || numeric_mismatch || count_mismatch || monte_carlo_failure) {
  stop("Historical-control LAVA output did not reproduce the archived output")
}
cat(
  "Historical-control reproduction passed: exact deterministic outputs and ",
  "concordant unseeded Monte Carlo inference.\n",
  sep = ""
)
