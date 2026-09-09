#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript tests/validate_synthetic_output.R <SYNTHETIC_OUTPUT_DIR>")
}

out_dir <- normalizePath(args[[1L]], mustWork = TRUE)
prefix <- "SYNTHETIC_FOURTRAIT"

master <- fread(file.path(out_dir, paste0(prefix, "_PLEIO_locus_master.tsv")))
qc <- fread(file.path(out_dir, paste0(prefix, "_classification_QC.tsv")))

expected_direction <- c(
  PLEIO_1 = "concordant",
  PLEIO_2 = "mixed",
  PLEIO_3 = "discordant",
  PLEIO_4 = "trait_dominant_or_unclear"
)
expected_prioritised <- c(
  PLEIO_1 = FALSE,
  PLEIO_2 = FALSE,
  PLEIO_3 = TRUE,
  PLEIO_4 = TRUE
)

if (!setequal(master$locus_key, names(expected_direction))) {
  stop("Synthetic master does not contain the expected locus keys.")
}

observed_direction <- setNames(master$locus_direction_class, master$locus_key)
observed_prioritised <- setNames(
  master$pleio_prioritised_vs_relevant_single_trait_FUMA,
  master$locus_key
)

if (!identical(observed_direction[names(expected_direction)], expected_direction)) {
  stop("Synthetic locus-direction labels differ from the expected labels.")
}
if (!identical(
  observed_prioritised[names(expected_prioritised)],
  expected_prioritised
)) {
  stop("Synthetic single-trait prioritisation flags differ from expected values.")
}

expected_qc <- list(
  n_pleio_loci = 4L,
  n_ind_sig_snps = 4L,
  n_ind_sig_snps_missing_raw_match = 0L,
  n_ind_sig_snps_unassigned_locus = 0L,
  n_loci_unassigned_corrected = 0L,
  n_corrected_prioritised = 2L
)

for (field in names(expected_qc)) {
  if (!field %in% names(qc) || nrow(qc) != 1L || qc[[field]][[1L]] != expected_qc[[field]]) {
    stop("Unexpected QC value for ", field, ".")
  }
}

message("Synthetic validation passed: 4 loci and all expected classifications confirmed.")
