#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_genomic_sem_models.R <10k|22k> <munged_dir> <output_dir>")
}

freeze <- args[[1L]]
munged_dir <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

trait_ids <- c(
  "reaction_time",
  "stroop_box",
  "stroop_ink",
  "quiz_score",
  "working_memory",
  "pairing_guesses"
)

files <- file.path(munged_dir, paste0(freeze, "_", trait_ids, ".sumstats.gz"))
if (!all(file.exists(files))) {
  stop("Missing munged inputs: ", paste(files[!file.exists(files)], collapse = ", "))
}

ld_path <- "/home/ec2-user/input_files/reference/eas_ldscores/"
log_path <- file.path(output_dir, paste0(freeze, "_multivariable_ldsc"))

sink(file.path(output_dir, "session_info.txt"), split = TRUE)
print(sessionInfo())
sink()

ldsc_output <- GenomicSEM::ldsc(
  traits = files,
  sample.prev = rep(NA_real_, length(files)),
  population.prev = rep(NA_real_, length(files)),
  ld = ld_path,
  wld = ld_path,
  trait.names = trait_ids,
  ldsc.log = log_path,
  stand = FALSE
)

saveRDS(ldsc_output, file.path(output_dir, "ldsc_output.rds"))
write.table(ldsc_output$S, file.path(output_dir, "genetic_covariance.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(stats::cov2cor(ldsc_output$S), file.path(output_dir, "genetic_correlation.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(ldsc_output$V, file.path(output_dir, "sampling_covariance.tsv"), sep = "\t", quote = FALSE, col.names = NA)
if (!is.null(ldsc_output$I)) {
  write.table(ldsc_output$I, file.path(output_dir, "intercept_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
}

one_factor_model <- '
g =~ NA*reaction_time + stroop_box + stroop_ink + quiz_score + working_memory + pairing_guesses
g ~~ 1*g
'

two_factor_model <- '
latency =~ NA*reaction_time + stroop_box + stroop_ink
latency ~~ 1*latency
non_latency =~ NA*quiz_score + working_memory + pairing_guesses
non_latency ~~ 1*non_latency
latency ~~ non_latency
'

fit_model <- function(model_name, syntax) {
  warnings_seen <- character()
  fitted <- withCallingHandlers(
    GenomicSEM::usermodel(
      covstruc = ldsc_output,
      estimation = "DWLS",
      model = syntax,
      CFIcalc = TRUE,
      std.lv = TRUE,
      imp_cov = TRUE,
      fix_resid = FALSE
    ),
    warning = function(warning_condition) {
      warnings_seen <<- c(warnings_seen, conditionMessage(warning_condition))
      invokeRestart("muffleWarning")
    }
  )
  saveRDS(fitted, file.path(output_dir, paste0(model_name, "_fit.rds")))
  if (!is.null(fitted$results)) {
    write.table(fitted$results, file.path(output_dir, paste0(model_name, "_parameters.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(fitted$modelfit)) {
    model_fit <- as.data.frame(fitted$modelfit)
    model_fit$model <- model_name
    write.table(model_fit, file.path(output_dir, paste0(model_name, "_fit_indices.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(fitted$imp_cov)) {
    write.table(fitted$imp_cov, file.path(output_dir, paste0(model_name, "_implied_covariance.tsv")), sep = "\t", quote = FALSE, col.names = NA)
  }
  writeLines(unique(warnings_seen), file.path(output_dir, paste0(model_name, "_warnings.txt")))
  fitted
}

one_factor <- fit_model("one_factor", one_factor_model)
two_factor <- fit_model("two_factor", two_factor_model)

fit_rows <- list()
if (!is.null(one_factor$modelfit)) {
  row <- as.data.frame(one_factor$modelfit)
  row$model <- "one_factor"
  fit_rows[[length(fit_rows) + 1L]] <- row
}
if (!is.null(two_factor$modelfit)) {
  row <- as.data.frame(two_factor$modelfit)
  row$model <- "two_factor"
  fit_rows[[length(fit_rows) + 1L]] <- row
}
if (length(fit_rows)) {
  all_names <- unique(unlist(lapply(fit_rows, names)))
  fit_rows <- lapply(fit_rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (name in missing) row[[name]] <- NA
    row[all_names]
  })
  write.table(do.call(rbind, fit_rows), file.path(output_dir, "model_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("Completed", freeze, "LDSC and Genomic SEM models in", output_dir, "\n")
