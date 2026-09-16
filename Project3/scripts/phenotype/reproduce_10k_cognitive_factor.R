#!/usr/bin/env Rscript

# Reproduce the six-task cognitive factor used in the HELIOS 10k analysis.
#
# Usage:
#   Rscript reproduce_10k_cognitive_factor.R \
#     Cognition_phenotypes_10K_Pritesh.txt \
#     merged_pheno_cov_filtered_10k.txt \
#     output_directory \
#     [g_scores_10k.txt]
#
# The optional fourth argument validates the reproduced score against the
# historical score file. Participant-level outputs must not be committed.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L || length(args) > 4L) {
  stop(
    paste(
      "Usage: Rscript reproduce_10k_cognitive_factor.R",
      "<cognitive_file> <covariate_file> <output_directory>",
      "[historical_g_file]"
    ),
    call. = FALSE
  )
}

cognitive_file <- args[[1L]]
covariate_file <- args[[2L]]
output_directory <- args[[3L]]
historical_g_file <- if (length(args) == 4L) args[[4L]] else NULL

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

mode_value <- function(x) {
  observed <- x[!is.na(x)]
  if (length(observed) == 0L) {
    return(NA)
  }
  values <- unique(observed)
  values[[which.max(tabulate(match(observed, values)))]]
}

read_whitespace_table <- function(path) {
  read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = ""
  )
}

require_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        label,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

cognitive <- read_whitespace_table(cognitive_file)
covariates <- read.table(
  covariate_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

task_columns <- c(
  "DC7R4_Pairing7_Guesses",
  "DC7R6_React_Avg",
  "DC7R8_Stroopbox_Avg",
  "DC7R10_Stroopink_Avg",
  "DC7R12_Wm_Score",
  "DC7R5_Quiz_Score"
)

require_columns(cognitive, c("IID", task_columns), "Cognitive input")
require_columns(
  covariates,
  c("IID", "Age", "Sex", "Ancestry", "FSAQ21_Soc9"),
  "Covariate input"
)

cognitive$IID <- trimws(as.character(cognitive$IID))
covariates$IID <- trimws(as.character(covariates$IID))

# The historical merge retained the first record when an IID occurred more
# than once in the covariate extract.
cognitive <- cognitive[!duplicated(cognitive$IID) & nzchar(cognitive$IID), ]
covariates <- covariates[!duplicated(covariates$IID) & nzchar(covariates$IID), ]

task_data <- cognitive[, c("IID", task_columns)]
for (column in task_columns) {
  task_data[[column]] <- as.numeric(task_data[[column]])
}
rownames(task_data) <- task_data$IID
task_data$IID <- NULL

# Direction-align timed/error measures so that higher values indicate better
# performance before PCA.
task_data$DC7R6_React_Avg_inv <- 1 / task_data$DC7R6_React_Avg
task_data$DC7R8_Stroopbox_Avg_inv <- 1 / task_data$DC7R8_Stroopbox_Avg
task_data$DC7R10_Stroopink_Avg_inv <- 1 / task_data$DC7R10_Stroopink_Avg
task_data$DC7R4_Pairing7_Guesses_inv <-
  1 / task_data$DC7R4_Pairing7_Guesses

task_data <- task_data[, c(
  "DC7R12_Wm_Score",
  "DC7R5_Quiz_Score",
  "DC7R6_React_Avg_inv",
  "DC7R8_Stroopbox_Avg_inv",
  "DC7R10_Stroopink_Avg_inv",
  "DC7R4_Pairing7_Guesses_inv"
)]
task_data[!is.finite(as.matrix(task_data))] <- NA_real_
task_data <- task_data[complete.cases(task_data), , drop = FALSE]

# Apply the historical four-standard-deviation mask separately to each task.
for (column in names(task_data)) {
  standardized <- as.numeric(scale(task_data[[column]]))
  task_data[[column]][abs(standardized) > 4] <- NA_real_
}

covariates$age <- as.numeric(covariates$Age)
covariates$sex_numeric <- as.numeric(covariates$Sex)
covariates$age2 <- covariates$age^2
covariates$age_sex <- covariates$age * covariates$sex_numeric
covariates$age2_sex <- covariates$age2 * covariates$sex_numeric
covariates$HI <- as.numeric(covariates$FSAQ21_Soc9)
covariates$SES_dontknow <- as.integer(covariates$HI == -888)
covariates$SES_prefno <- as.integer(covariates$HI == -777)
covariates$indian <- as.integer(covariates$Ancestry == "Indian")
covariates$malay <- as.integer(covariates$Ancestry == "Malay")

factor_columns <- c(
  "Sex", "indian", "malay", "HI", "SES_dontknow", "SES_prefno"
)
for (column in factor_columns) {
  replacement <- mode_value(covariates[[column]])
  covariates[[column]][is.na(covariates[[column]])] <- replacement
  covariates[[column]] <- factor(covariates[[column]])
}
covariates$sex <- covariates$Sex

common_ids <- intersect(rownames(task_data), covariates$IID)
task_data <- task_data[common_ids, , drop = FALSE]
model_data <- covariates[match(common_ids, covariates$IID), , drop = FALSE]

residuals_by_task <- matrix(
  NA_real_,
  nrow = nrow(task_data),
  ncol = ncol(task_data),
  dimnames = list(common_ids, names(task_data))
)

model_formula <-
  response ~ age + sex + age_sex + age2 + age2_sex +
  indian + malay + HI + SES_prefno + SES_dontknow

for (column in names(task_data)) {
  analysis_data <- model_data
  analysis_data$response <- task_data[[column]]
  fitted_model <- lm(
    model_formula,
    data = analysis_data,
    na.action = na.exclude
  )
  residuals_by_task[, column] <- resid(fitted_model)
}

residuals_by_task <- residuals_by_task[
  complete.cases(residuals_by_task) &
    apply(residuals_by_task, 1L, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]

# The historical code standardized the residuals before calling prcomp with
# scaling enabled. This is retained to reproduce that calculation exactly.
scaled_residuals <- scale(residuals_by_task)
pca <- prcomp(scaled_residuals, center = TRUE, scale. = TRUE)
scores <- pca$x[, 1L]
loadings <- pca$rotation[, 1L]

# PCA signs are arbitrary. The stored historical score had a negative loading
# on direction-aligned Stroop-box performance, so enforce that orientation.
if (loadings[["DC7R8_Stroopbox_Avg_inv"]] > 0) {
  scores <- -scores
  loadings <- -loadings
}

score_output <- data.frame(
  IID = rownames(residuals_by_task),
  g = as.numeric(scores),
  stringsAsFactors = FALSE
)
write.table(
  score_output,
  file.path(output_directory, "cognitive_factor_scores.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

task_labels <- c(
  DC7R12_Wm_Score = "Working memory",
  DC7R5_Quiz_Score = "Quiz score",
  DC7R6_React_Avg_inv = "Reaction time",
  DC7R8_Stroopbox_Avg_inv = "Stroop box",
  DC7R10_Stroopink_Avg_inv = "Stroop ink",
  DC7R4_Pairing7_Guesses_inv = "Pairing guesses"
)

component_correlations <- vapply(
  names(loadings),
  function(column) cor(scaled_residuals[, column], scores),
  numeric(1L)
)

loading_output <- data.frame(
  task = unname(task_labels[names(loadings)]),
  variable = names(loadings),
  loading_on_stored_g = as.numeric(loadings),
  correlation_with_stored_g = as.numeric(component_correlations),
  loading_higher_performance = -as.numeric(loadings),
  correlation_higher_performance = -as.numeric(component_correlations),
  stringsAsFactors = FALSE
)
write.table(
  loading_output,
  file.path(output_directory, "pca_loadings.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

correlation_output <- cor(scaled_residuals)
write.table(
  cbind(variable = rownames(correlation_output), correlation_output),
  file.path(output_directory, "task_correlations.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

variance_explained <- (pca$sdev^2) / sum(pca$sdev^2)
summary_output <- data.frame(
  metric = c("pca_sample_n", "pc1_variance_explained"),
  value = c(nrow(residuals_by_task), variance_explained[[1L]]),
  stringsAsFactors = FALSE
)

if (!is.null(historical_g_file)) {
  historical <- read_whitespace_table(historical_g_file)
  require_columns(historical, c("IID", "g"), "Historical score input")
  historical$IID <- trimws(as.character(historical$IID))
  historical$g <- as.numeric(historical$g)
  validation_ids <- intersect(score_output$IID, historical$IID)
  reproduced <- score_output$g[match(validation_ids, score_output$IID)]
  observed <- historical$g[match(validation_ids, historical$IID)]
  validation_correlation <- cor(reproduced, observed)
  validation_slope <- coef(lm(observed ~ 0 + reproduced))[[1L]]
  aligned <- reproduced * validation_slope
  validation_rmse <- sqrt(mean((observed - aligned)^2))
  validation_max_difference <- max(abs(observed - aligned))

  summary_output <- rbind(
    summary_output,
    data.frame(
      metric = c(
        "validation_sample_n",
        "historical_score_correlation",
        "historical_score_slope",
        "historical_score_rmse",
        "historical_score_max_absolute_difference"
      ),
      value = c(
        length(validation_ids),
        validation_correlation,
        validation_slope,
        validation_rmse,
        validation_max_difference
      ),
      stringsAsFactors = FALSE
    )
  )
}

write.table(
  summary_output,
  file.path(output_directory, "pca_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message(sprintf("PCA complete: %d participants", nrow(residuals_by_task)))
message(
  "Participant-level output written to cognitive_factor_scores.tsv; ",
  "do not add this file to GitHub."
)
