#!/usr/bin/env Rscript

# Prepare one ancestry-specific cognitive-factor phenotype for the historical
# HELIOS REGENIE workflow. Participant-level output is controlled data and
# must be written outside the public repository.
#
# Usage:
#   Rscript prepare_ancestry_g_for_regenie.R \
#     pooled_scores.tsv ancestry_ids.tsv ancestry_g.pheno [score_column] [sd_cutoff]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 5L) {
  stop(
    paste(
      "Usage: Rscript prepare_ancestry_g_for_regenie.R",
      "<pooled-score-file> <ancestry-ID-file> <output-file>",
      "[score-column; default: g] [SD-cutoff; default: 5]"
    ),
    call. = FALSE
  )
}

score_file <- args[[1L]]
ancestry_file <- args[[2L]]
output_file <- args[[3L]]
score_column <- if (length(args) >= 4L) args[[4L]] else "g"
sd_cutoff <- if (length(args) >= 5L) as.numeric(args[[5L]]) else 5

if (!is.finite(sd_cutoff) || sd_cutoff <= 0) {
  stop("SD cutoff must be a positive number", call. = FALSE)
}

read_input <- function(path) {
  read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = ""
  )
}

find_column <- function(columns, candidates, label) {
  found <- candidates[candidates %in% columns]
  if (length(found) == 0L) {
    stop(
      sprintf("No %s column found; tried: %s", label, paste(candidates, collapse = ", ")),
      call. = FALSE
    )
  }
  found[[1L]]
}

scores <- read_input(score_file)
ancestry <- read_input(ancestry_file)

score_id <- find_column(names(scores), c("IID", "ID"), "score IID")
ancestry_id <- find_column(names(ancestry), c("IID", "ID"), "ancestry IID")

if (!score_column %in% names(scores)) {
  stop(sprintf("Score column '%s' is absent", score_column), call. = FALSE)
}

scores$IID <- trimws(as.character(scores[[score_id]]))
ancestry$IID <- trimws(as.character(ancestry[[ancestry_id]]))

if (any(!nzchar(scores$IID)) || anyDuplicated(scores$IID)) {
  stop("Score IIDs must be non-empty and unique", call. = FALSE)
}
if (any(!nzchar(ancestry$IID)) || anyDuplicated(ancestry$IID)) {
  stop("Ancestry IIDs must be non-empty and unique", call. = FALSE)
}

score_values <- suppressWarnings(as.numeric(scores[[score_column]]))
if (any(!is.na(scores[[score_column]]) & is.na(score_values))) {
  stop(sprintf("Score column '%s' must be numeric", score_column), call. = FALSE)
}

keep <- scores$IID %in% ancestry$IID
analysis <- data.frame(
  IID = scores$IID[keep],
  score = score_values[keep],
  stringsAsFactors = FALSE
)

if (nrow(analysis) == 0L) {
  stop("No score IIDs overlap the ancestry ID file", call. = FALSE)
}

observed <- is.finite(analysis$score)
if (sum(observed) < 2L) {
  stop("Fewer than two finite scores remain after ancestry matching", call. = FALSE)
}

score_mean <- mean(analysis$score[observed])
score_sd <- sd(analysis$score[observed])
if (!is.finite(score_sd) || score_sd == 0) {
  stop("Score standard deviation is zero or undefined", call. = FALSE)
}

outlier <- observed & abs(analysis$score - score_mean) > sd_cutoff * score_sd
analysis$score[outlier] <- NA_real_

# Blom rank-based inverse-normal transformation used in the recovered
# ancestry-specific phenotype-preparation workflow.
remaining <- is.finite(analysis$score)
n_remaining <- sum(remaining)
analysis$g <- NA_real_
analysis$g[remaining] <- qnorm(
  (rank(analysis$score[remaining], ties.method = "average") - 3 / 8) /
    (n_remaining + 1 / 4)
)

fid_column <- find_column(names(ancestry), c("FID", "#FID", "IID", "ID"), "ancestry FID")
fid <- as.character(ancestry[[fid_column]][match(analysis$IID, ancestry$IID)])
fid[is.na(fid) | !nzchar(fid)] <- analysis$IID[is.na(fid) | !nzchar(fid)]

output <- data.frame(FID = fid, IID = analysis$IID, g = analysis$g)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(
  output,
  output_file,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA"
)

message(sprintf(
  "Matched %d participants; masked %d values beyond %.1f SD; transformed %d scores",
  nrow(output), sum(outlier), sd_cutoff, n_remaining
))
