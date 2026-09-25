#!/usr/bin/env Rscript

# Quantify the effect of correcting KoGES EDU from binary to quantitative in LAVA.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "Usage: compare_control_corrected.R control.univ control.bivar ",
    "corrected.univ corrected.bivar summary.tsv eligibility.tsv"
  )
}

HELIOS <- "HEL10k_EAS"
EDU <- "Edu_EAS"
UNIVARIATE_THRESHOLD <- 2e-5

control_univ <- fread(args[1], na.strings = c("NA", "NaN"))
control_bivar <- fread(args[2], na.strings = c("NA", "NaN"))
corrected_univ <- fread(args[3], na.strings = c("NA", "NaN"))
corrected_bivar <- fread(args[4], na.strings = c("NA", "NaN"))

target_pair <- function(frame, prefix) {
  answer <- frame[
    (phen1 == HELIOS & phen2 == EDU) | (phen1 == EDU & phen2 == HELIOS),
    .(locus, rho, p)
  ]
  setnames(answer, c("rho", "p"), paste0(prefix, c("_rho", "_p")))
  answer
}

control_target <- target_pair(control_bivar, "control")
corrected_target <- target_pair(corrected_bivar, "corrected")
eligibility <- merge(control_target, corrected_target, by = "locus", all = TRUE)
eligibility[, `:=`(
  in_control = !is.na(control_rho),
  in_corrected = !is.na(corrected_rho)
)]
setcolorder(eligibility, c(
  "locus", "in_control", "in_corrected", "control_rho", "corrected_rho",
  "control_p", "corrected_p"
))
setorder(eligibility, locus)

univ_comparison <- merge(
  control_univ[phen %chin% c(HELIOS, EDU), .(locus, phen, control_h2 = h2.obs, control_p = p)],
  corrected_univ[phen %chin% c(HELIOS, EDU), .(locus, phen, corrected_h2 = h2.obs, corrected_p = p)],
  by = c("locus", "phen"), all = TRUE
)

maximum_difference <- function(first, second) {
  difference <- abs(first - second)
  if (all(is.na(difference))) NA_real_ else max(difference, na.rm = TRUE)
}

common_target <- eligibility[in_control & in_corrected]
summary <- data.table(
  metric = c(
    "control_univariate_rows_all_traits",
    "corrected_univariate_rows_all_traits",
    "control_bivariate_rows_all_pairs",
    "corrected_bivariate_rows_all_pairs",
    "control_helios_locally_heritable_blocks",
    "corrected_helios_locally_heritable_blocks",
    "control_edu_locally_heritable_blocks",
    "corrected_edu_locally_heritable_blocks",
    "control_helios_edu_eligible_blocks",
    "corrected_helios_edu_eligible_blocks",
    "helios_edu_eligible_in_both",
    "helios_edu_control_only",
    "helios_edu_corrected_only",
    "helios_edu_eligibility_jaccard",
    "common_target_rho_correlation",
    "common_target_median_abs_rho_difference",
    "common_target_max_abs_rho_difference",
    "common_target_max_abs_p_difference",
    "common_target_sign_discordant",
    "control_target_nominal_significant",
    "corrected_target_nominal_significant",
    "common_target_nominal_decision_discordant",
    "control_target_fdr_significant",
    "corrected_target_fdr_significant",
    "common_target_fdr_decision_discordant",
    "helios_univariate_max_abs_h2_difference",
    "helios_univariate_max_abs_p_difference",
    "edu_univariate_max_abs_h2_difference",
    "edu_univariate_max_abs_p_difference"
  ),
  value = c(
    nrow(control_univ),
    nrow(corrected_univ),
    nrow(control_bivar),
    nrow(corrected_bivar),
    control_univ[phen == HELIOS & p < UNIVARIATE_THRESHOLD, .N],
    corrected_univ[phen == HELIOS & p < UNIVARIATE_THRESHOLD, .N],
    control_univ[phen == EDU & p < UNIVARIATE_THRESHOLD, .N],
    corrected_univ[phen == EDU & p < UNIVARIATE_THRESHOLD, .N],
    nrow(control_target),
    nrow(corrected_target),
    eligibility[in_control & in_corrected, .N],
    eligibility[in_control & !in_corrected, .N],
    eligibility[!in_control & in_corrected, .N],
    eligibility[in_control & in_corrected, .N] / nrow(eligibility),
    if (nrow(common_target) > 1L) cor(common_target$control_rho, common_target$corrected_rho) else NA_real_,
    median(abs(common_target$control_rho - common_target$corrected_rho)),
    maximum_difference(common_target$control_rho, common_target$corrected_rho),
    maximum_difference(common_target$control_p, common_target$corrected_p),
    sum(sign(common_target$control_rho) != sign(common_target$corrected_rho)),
    sum(control_target$control_p < 0.05),
    sum(corrected_target$corrected_p < 0.05),
    sum((common_target$control_p < 0.05) != (common_target$corrected_p < 0.05)),
    sum(p.adjust(control_target$control_p, method = "BH") < 0.05),
    sum(p.adjust(corrected_target$corrected_p, method = "BH") < 0.05),
    sum(
      (p.adjust(common_target$control_p, method = "BH") < 0.05) !=
        (p.adjust(common_target$corrected_p, method = "BH") < 0.05)
    ),
    univ_comparison[phen == HELIOS, maximum_difference(control_h2, corrected_h2)],
    univ_comparison[phen == HELIOS, maximum_difference(control_p, corrected_p)],
    univ_comparison[phen == EDU, maximum_difference(control_h2, corrected_h2)],
    univ_comparison[phen == EDU, maximum_difference(control_p, corrected_p)]
  )
)

dir.create(dirname(args[5]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args[6]), recursive = TRUE, showWarnings = FALSE)
fwrite(summary, args[5], sep = "\t", na = "NA", quote = FALSE)
fwrite(eligibility, args[6], sep = "\t", na = "NA", quote = FALSE)
print(summary)
