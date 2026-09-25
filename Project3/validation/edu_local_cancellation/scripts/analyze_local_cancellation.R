#!/usr/bin/env Rscript

# Analyse HELIOS cognition–KoGES EDU local covariance from corrected LAVA output.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: analyze_local_cancellation.R corrected.univ corrected.bivar ",
    "loci.txt global_summary.tsv output_root"
  )
}

HELIOS <- "HEL10k_EAS"
EDU <- "Edu_EAS"
UNIVARIATE_THRESHOLD <- 2e-5
MHC_CHR <- 6L
MHC_START <- 25000000L
MHC_STOP <- 34000000L

output_root <- normalizePath(args[5], mustWork = FALSE)
directories <- file.path(
  output_root,
  c("03_local_covariance", "04_cancellation", "05_sensitivity", "06_plots", "final")
)
invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))
local_directory <- directories[1]
cancellation_directory <- directories[2]
sensitivity_directory <- directories[3]
plot_directory <- directories[4]
final_directory <- directories[5]

univariate <- fread(args[1], na.strings = c("NA", "NaN"))
bivariate <- fread(args[2], na.strings = c("NA", "NaN"))
loci <- fread(args[3], na.strings = c("NA", "NaN"))
global_summary <- fread(args[4], na.strings = c("NA", "NaN"))

target <- bivariate[
  (phen1 == HELIOS & phen2 == EDU) | (phen1 == EDU & phen2 == HELIOS)
]
if (nrow(target) == 0L) stop("No HELIOS–EDU bivariate rows were found")
if (anyDuplicated(target$locus)) stop("Duplicate HELIOS–EDU bivariate rows were found")

helios_univ <- univariate[phen == HELIOS]
edu_univ <- univariate[phen == EDU]
helios_columns <- helios_univ[, .(
  locus,
  helios_local_h2 = h2.obs,
  helios_local_h2_p = p
)]
edu_columns <- edu_univ[, .(
  locus,
  edu_local_h2 = h2.obs,
  edu_local_h2_p = p
)]
target <- merge(target, helios_columns, by = "locus", all.x = TRUE)
target <- merge(target, edu_columns, by = "locus", all.x = TRUE)
if (target[, anyNA(.SD), .SDcols = c(
  "helios_local_h2", "helios_local_h2_p", "edu_local_h2", "edu_local_h2_p"
)]) {
  stop("Missing univariate result for an eligible bivariate row")
}

target[, local_rg_fdr := p.adjust(p, method = "BH")]
target[, local_cov_original := rho * sqrt(helios_local_h2 * edu_local_h2)]
target[, `:=`(
  local_cov_higher_performance = -local_cov_original,
  local_rg_higher_performance = -rho,
  local_rg_higher_performance_lower = -rho.upper,
  local_rg_higher_performance_upper = -rho.lower
)]
target[, sign_rg := fifelse(
  local_rg_higher_performance > 0, "positive",
  fifelse(local_rg_higher_performance < 0, "negative", "zero")
)]
target[, sign_covariance := fifelse(
  local_cov_higher_performance > 0, "positive",
  fifelse(local_cov_higher_performance < 0, "negative", "zero")
)]

regions <- target[, .(
  locus_id = as.integer(locus),
  chr = as.integer(chr),
  start = as.integer(start),
  stop = as.integer(stop),
  n_snps = as.integer(n.snps),
  n_pcs = as.integer(n.pcs),
  helios_local_h2,
  helios_local_h2_se = NA_real_,
  helios_local_h2_p,
  edu_local_h2,
  edu_local_h2_se = NA_real_,
  edu_local_h2_p,
  local_rg_original = rho,
  local_rg_original_lower = rho.lower,
  local_rg_original_upper = rho.upper,
  local_rg_higher_performance,
  local_rg_higher_performance_lower,
  local_rg_higher_performance_upper,
  local_rg_se = NA_real_,
  local_rg_p = p,
  local_rg_fdr,
  local_cov_original,
  local_cov_higher_performance,
  covariance_source = "derived_from_rho_and_observed_h2",
  eligible_for_bivariate = TRUE,
  sign_rg,
  sign_covariance,
  source_file = normalizePath(args[2])
)]
setorder(regions, chr, start, stop)
fwrite(
  regions,
  file.path(local_directory, "helios_koges_all_eligible_regions.tsv"),
  sep = "\t", na = "NA", quote = FALSE
)

aggregate_rg <- function(frame) {
  denominator <- sqrt(sum(frame$helios_local_h2) * sum(frame$edu_local_h2))
  if (is.finite(denominator) && denominator > 0) {
    sum(frame$local_cov_higher_performance) / denominator
  } else {
    NA_real_
  }
}

decomposition <- function(frame, label) {
  covariance <- frame$local_cov_higher_performance
  positive <- sum(pmax(covariance, 0))
  negative <- sum(abs(pmin(covariance, 0)))
  signed <- sum(covariance)
  absolute <- sum(abs(covariance))
  ratio <- if (absolute > 0) abs(signed) / absolute else NA_real_
  data.table(
    analysis_set = label,
    n_regions = nrow(frame),
    sum_positive_cov = positive,
    sum_negative_abs_cov = negative,
    sum_signed_cov = signed,
    sum_absolute_cov = absolute,
    net_to_total_ratio = ratio,
    cancellation_fraction = if (is.na(ratio)) NA_real_ else 1 - ratio
  )
}

describe_set <- function(frame, label) {
  data.table(
    analysis_set = label,
    n_regions = nrow(frame),
    n_positive_covariance = sum(frame$local_cov_higher_performance > 0),
    n_negative_covariance = sum(frame$local_cov_higher_performance < 0),
    n_zero_covariance = sum(frame$local_cov_higher_performance == 0),
    median_local_rg = median(frame$local_rg_higher_performance),
    local_rg_q1 = quantile(frame$local_rg_higher_performance, 0.25, names = FALSE),
    local_rg_q3 = quantile(frame$local_rg_higher_performance, 0.75, names = FALSE),
    median_local_covariance = median(frame$local_cov_higher_performance),
    local_covariance_q1 = quantile(frame$local_cov_higher_performance, 0.25, names = FALSE),
    local_covariance_q3 = quantile(frame$local_cov_higher_performance, 0.75, names = FALSE),
    aggregate_rg_eligible = aggregate_rg(frame)
  )
}

subsets <- list(
  all_eligible = regions,
  nominal_significant = regions[local_rg_p < 0.05],
  fdr_significant = regions[local_rg_fdr < 0.05]
)
descriptive <- rbindlist(Map(describe_set, subsets, names(subsets)))
fwrite(descriptive, file.path(local_directory, "local_covariance_summary.tsv"), sep = "\t", na = "NA", quote = FALSE)
decomposition_table <- rbindlist(Map(decomposition, subsets, names(subsets)))
fwrite(
  decomposition_table,
  file.path(cancellation_directory, "covariance_decomposition.tsv"),
  sep = "\t", na = "NA", quote = FALSE
)

chromosome_table <- rbindlist(lapply(sort(unique(regions$chr)), function(chromosome) {
  row <- decomposition(regions[chr == chromosome], paste0("chr", chromosome))
  row[, chr := chromosome]
  row
}), fill = TRUE)
setnames(
  chromosome_table,
  c("n_regions", "sum_positive_cov", "sum_negative_abs_cov", "sum_signed_cov", "sum_absolute_cov"),
  c("n_eligible_regions", "positive_covariance_sum", "negative_abs_covariance_sum", "signed_covariance_sum", "absolute_covariance_sum")
)
setcolorder(chromosome_table, c(
  "chr", "analysis_set", "n_eligible_regions", "positive_covariance_sum",
  "negative_abs_covariance_sum", "signed_covariance_sum", "absolute_covariance_sum",
  "net_to_total_ratio", "cancellation_fraction"
))
fwrite(
  chromosome_table,
  file.path(cancellation_directory, "chromosome_covariance_summary.tsv"),
  sep = "\t", na = "NA", quote = FALSE
)

leave_one_out <- rbindlist(lapply(1:22, function(chromosome) {
  frame <- regions[chr != chromosome]
  row <- decomposition(frame, paste0("exclude_chr", chromosome))
  row[, `:=`(excluded_chr = chromosome, rg_eligible = aggregate_rg(frame))]
  row
}))
fwrite(
  leave_one_out,
  file.path(sensitivity_directory, "leave_one_chromosome_out.tsv"),
  sep = "\t", na = "NA", quote = FALSE
)

ranked <- copy(regions)
ranked[, absolute_covariance := abs(local_cov_higher_performance)]
setorder(ranked, -absolute_covariance)
outlier_exclusion <- rbindlist(lapply(c(0L, 1L, 3L, 5L), function(count) {
  frame <- if (count == 0L) ranked else ranked[-seq_len(count)]
  row <- decomposition(frame, paste0("exclude_top_", count, "_absolute_covariance"))
  excluded <- if (count == 0L) "" else paste(ranked[seq_len(count), locus_id], collapse = ";")
  row[, `:=`(excluded_loci = excluded, rg_eligible = aggregate_rg(frame))]
  row
}))
fwrite(
  outlier_exclusion,
  file.path(sensitivity_directory, "outlier_exclusion.tsv"),
  sep = "\t", na = "NA", quote = FALSE
)

mhc <- regions$chr == MHC_CHR & regions$start <= MHC_STOP & regions$stop >= MHC_START
mhc_exclusion <- rbind(
  decomposition(regions, "mhc_included")[, `:=`(
    mhc_regions_removed = 0L,
    rg_eligible = aggregate_rg(regions)
  )],
  decomposition(regions[!mhc], "mhc_excluded")[, `:=`(
    mhc_regions_removed = sum(mhc),
    rg_eligible = aggregate_rg(regions[!mhc])
  )]
)
fwrite(
  mhc_exclusion,
  file.path(sensitivity_directory, "mhc_exclusion.tsv"),
  sep = "\t", na = "NA", quote = FALSE
)

rank_table <- function(frame, label) {
  answer <- copy(frame[seq_len(min(10L, nrow(frame)))])
  answer[, `:=`(ranking = label, rank = seq_len(.N))]
  setcolorder(answer, c("ranking", "rank", setdiff(names(answer), c("ranking", "rank"))))
  answer
}
positive_rank <- copy(regions); setorder(positive_rank, -local_cov_higher_performance)
negative_rank <- copy(regions); setorder(negative_rank, local_cov_higher_performance)
absolute_rg_rank <- copy(regions); absolute_rg_rank[, rank_value := abs(local_rg_higher_performance)]
setorder(absolute_rg_rank, -rank_value); absolute_rg_rank[, rank_value := NULL]
top_regions <- rbindlist(list(
  rank_table(positive_rank, "largest_positive_covariance"),
  rank_table(negative_rank, "largest_negative_covariance"),
  rank_table(absolute_rg_rank, "largest_absolute_local_rg")
), fill = TRUE)
fwrite(top_regions, file.path(local_directory, "top_regions.tsv"), sep = "\t", na = "NA", quote = FALSE)

helios_significant <- helios_univ[p < UNIVARIATE_THRESHOLD, .N]
edu_significant <- edu_univ[p < UNIVARIATE_THRESHOLD, .N]
helios_h2_all <- sum(helios_univ$h2.obs)
edu_h2_all <- sum(edu_univ$h2.obs)
helios_h2_eligible <- sum(regions$helios_local_h2)
edu_h2_eligible <- sum(regions$edu_local_h2)
helios_coverage <- if (helios_h2_all > 0) helios_h2_eligible / helios_h2_all else NA_real_
edu_coverage <- if (edu_h2_all > 0) edu_h2_eligible / edu_h2_all else NA_real_
eligibility <- data.table(
  metric = c(
    "total_lava_blocks", "helios_estimable_blocks", "edu_estimable_blocks",
    "helios_locally_heritable_blocks", "edu_locally_heritable_blocks",
    "bivariate_eligible_blocks", "percentage_all_blocks_eligible",
    "helios_h2_all_estimable", "helios_h2_bivariate_eligible", "helios_h2_coverage",
    "edu_h2_all_estimable", "edu_h2_bivariate_eligible", "edu_h2_coverage"
  ),
  value = c(
    nrow(loci), nrow(helios_univ), nrow(edu_univ), helios_significant,
    edu_significant, nrow(regions), 100 * nrow(regions) / nrow(loci),
    helios_h2_all, helios_h2_eligible, helios_coverage,
    edu_h2_all, edu_h2_eligible, edu_coverage
  )
)
fwrite(eligibility, file.path(local_directory, "eligibility_coverage.tsv"), sep = "\t", na = "NA", quote = FALSE)

metric_value <- function(target_metric, column = "value") {
  answer <- global_summary[get("metric") == target_metric, get(column)]
  if (length(answer) != 1L) stop("Global metric not found uniquely: ", target_metric, " / ", column)
  as.numeric(answer)
}
primary <- decomposition_table[analysis_set == "all_eligible"]
poster <- data.table(
  global_ldsc_rg_original_higher_worse = metric_value("genetic_correlation"),
  global_ldsc_rg_higher_performance = -metric_value("genetic_correlation"),
  global_ldsc_rg_se = metric_value("genetic_correlation", "se"),
  global_ldsc_p = metric_value("genetic_correlation", "p"),
  n_lava_blocks_total = nrow(loci),
  n_bivariate_eligible = nrow(regions),
  n_positive_cov = regions[local_cov_higher_performance > 0, .N],
  n_negative_cov = regions[local_cov_higher_performance < 0, .N],
  n_fdr_positive = regions[local_rg_fdr < 0.05 & local_cov_higher_performance > 0, .N],
  n_fdr_negative = regions[local_rg_fdr < 0.05 & local_cov_higher_performance < 0, .N],
  sum_positive_cov = primary$sum_positive_cov,
  sum_negative_abs_cov = primary$sum_negative_abs_cov,
  sum_signed_cov = primary$sum_signed_cov,
  sum_absolute_cov = primary$sum_absolute_cov,
  net_to_total_ratio = primary$net_to_total_ratio,
  cancellation_fraction = primary$cancellation_fraction,
  rg_eligible = aggregate_rg(regions),
  helios_h2_coverage = helios_coverage,
  edu_h2_coverage = edu_coverage
)
fwrite(poster, file.path(final_directory, "poster_key_results.tsv"), sep = "\t", na = "NA", quote = FALSE)

chromosome_boundaries <- loci[, .(
  first_locus = min(LOC),
  last_locus = max(LOC),
  centre = (min(LOC) + max(LOC)) / 2
), by = CHR]
setorder(chromosome_boundaries, CHR)

genomic_plot <- function(y_column, y_label, title, filename) {
  plot_data <- copy(regions)
  plot_data[, significant := local_rg_fdr < 0.05]
  plot_data[, direction := fifelse(get(y_column) >= 0, "Positive", "Negative")]
  p <- ggplot(plot_data, aes(x = locus_id, y = .data[[y_column]], colour = direction)) +
    geom_vline(
      xintercept = chromosome_boundaries$first_locus[-1],
      colour = "grey88", linewidth = 0.25
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
    geom_point(aes(shape = significant), alpha = 0.78, size = 1.8) +
    scale_colour_manual(values = c(Positive = "#2B6CB0", Negative = "#C53030")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), labels = c("FDR >= 0.05", "FDR < 0.05")) +
    scale_x_continuous(
      breaks = chromosome_boundaries$centre,
      labels = chromosome_boundaries$CHR,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(x = "Chromosome", y = y_label, title = title, colour = "Direction", shape = NULL) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
  ggsave(file.path(plot_directory, filename), p, width = 12, height = 5, dpi = 300)
}

genomic_plot(
  "local_rg_higher_performance", "Local genetic correlation",
  "HELIOS cognition and KoGES education: local genetic correlation",
  "helios_edu_local_rg_genome.png"
)
genomic_plot(
  "local_cov_higher_performance", "Local genetic covariance",
  "HELIOS cognition and KoGES education: local genetic covariance",
  "helios_edu_local_covariance_genome.png"
)

cumulative <- copy(regions)
cumulative[, cumulative_covariance := cumsum(local_cov_higher_performance)]
p_cumulative <- ggplot(cumulative, aes(locus_id, cumulative_covariance)) +
  geom_vline(
    xintercept = chromosome_boundaries$first_locus[-1],
    colour = "grey88", linewidth = 0.25
  ) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
  geom_line(colour = "#4A5568", linewidth = 0.75) +
  scale_x_continuous(
    breaks = chromosome_boundaries$centre,
    labels = chromosome_boundaries$CHR,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Chromosome", y = "Cumulative signed local covariance",
    title = "Cumulative local genetic covariance"
  ) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(plot_directory, "cumulative_signed_covariance.png"), p_cumulative, width = 12, height = 5, dpi = 300)

decomposition_plot <- data.table(
  contribution = factor(
    c("Positive", "Negative (absolute)", "Net (absolute)"),
    levels = c("Positive", "Negative (absolute)", "Net (absolute)")
  ),
  value = c(primary$sum_positive_cov, primary$sum_negative_abs_cov, abs(primary$sum_signed_cov))
)
p_decomposition <- ggplot(decomposition_plot, aes(contribution, value, fill = contribution)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = signif(value, 4)), vjust = -0.45, size = 3.3) +
  scale_fill_manual(values = c("#2B6CB0", "#C53030", "#4A5568")) +
  labs(
    x = NULL, y = "Summed local genetic covariance",
    title = "Opposing local covariance contributions",
    subtitle = sprintf("Cancellation fraction = %.3f", primary$cancellation_fraction)
  ) +
  expand_limits(y = max(decomposition_plot$value) * 1.12) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))
ggsave(
  file.path(plot_directory, "positive_negative_covariance_decomposition.png"),
  p_decomposition, width = 7.5, height = 5.2, dpi = 300
)

funnel <- data.table(
  stage = factor(
    c(
      "All LAVA blocks", "HELIOS locally heritable", "EDU locally heritable",
      "Both eligible", "Nominal bivariate", "FDR bivariate"
    ),
    levels = rev(c(
      "All LAVA blocks", "HELIOS locally heritable", "EDU locally heritable",
      "Both eligible", "Nominal bivariate", "FDR bivariate"
    ))
  ),
  regions = c(
    nrow(loci), helios_significant, edu_significant, nrow(regions),
    nrow(subsets$nominal_significant), nrow(subsets$fdr_significant)
  )
)
p_funnel <- ggplot(funnel, aes(stage, regions)) +
  geom_col(fill = "#2C7A7B", width = 0.7) +
  geom_text(aes(label = regions), hjust = -0.2, size = 3.2) +
  coord_flip() +
  expand_limits(y = max(funnel$regions) * 1.12) +
  labs(x = NULL, y = "Regions", title = "LAVA eligibility and evidence funnel") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(plot_directory, "lava_eligibility_funnel.png"), p_funnel, width = 9, height = 5.5, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(final_directory, "postprocessing_sessionInfo.txt"))
print(t(poster))
