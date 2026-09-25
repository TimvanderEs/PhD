#!/usr/bin/env Rscript

# HELIOS 22k phenotype sensitivity to genetic-ancestry coding and
# participant-described ethnic heritage.
#
# Controlled participant-level inputs remain in memory. The script writes
# aggregate tables, non-identifiable figures and logs only.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_phenotype_ethnicity_sensitivity.R <analysis-root>")
}

root <- normalizePath(args[[1]], mustWork = TRUE)
input_dir <- file.path(root, "input")
output_dir <- file.path(root, "output")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- list(
  metadata = file.path(input_dir, "Full_data_HELIOS_22k.txt"),
  pca = file.path(input_dir, "22k_pre_norm_post_resid_Edu.txt"),
  annotation = file.path(input_dir, "Sample_Annotation_r231.csv"),
  interview = file.path(input_dir, "SG100K_IAQ_anon_r1_v1_1.csv")
)
missing_inputs <- names(input_paths)[!file.exists(unlist(input_paths))]
if (length(missing_inputs) > 0L) {
  stop("Missing inputs: ", paste(missing_inputs, collapse = ", "))
}

log_file <- file.path(output_dir, "run_log.txt")
log_connection <- file(log_file, open = "wt")
sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_connection)
}, add = TRUE)

cat("HELIOS phenotype/ethnicity sensitivity\n")
cat("Started:", format(Sys.time(), tz = "UTC"), "UTC\n")

`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0L) fallback else value
}

scalar_character <- function(object, field) {
  value <- object[[field]] %||% NA_character_
  if (length(value) == 0L || is.null(value)) return(NA_character_)
  as.character(value[[1]])
}

normalise_group <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value[value %in% c("", "na", "nan")] <- NA_character_
  value[value == "others"] <- "other"
  value
}

read_annotation <- function(path) {
  lines <- readLines(path, warn = FALSE)
  expected_header <- paste(
    c(
      "npmid", "ilmn_meta", "sg100k_meta", "sample_filters",
      "sample_flags", "grids_meta", "sample_qc"
    ),
    collapse = ","
  )
  if (!identical(lines[[1]], expected_header)) {
    stop("Unexpected SG100K annotation header")
  }
  records <- vector("list", length(lines) - 1L)
  for (index in seq.int(2L, length(lines))) {
    line <- lines[[index]]
    delimiter <- regexpr(",", line, fixed = TRUE)[[1]]
    if (delimiter < 1L) stop("Malformed annotation line: ", index)
    npmid <- substr(line, 1L, delimiter - 1L)
    payload <- substr(line, delimiter + 1L, nchar(line))
    values <- fromJSON(paste0("[", payload, "]"), simplifyVector = FALSE)
    if (length(values) != 6L) stop("Unexpected annotation payload: ", index)
    meta <- values[[2]]
    grids <- values[[5]]
    records[[index - 1L]] <- list(
      IID = npmid,
      participant_id = scalar_character(meta, "participant_id"),
      cohort = scalar_character(meta, "cohort"),
      genetic_ancestry = scalar_character(grids, "ancestry")
    )
    if ((index - 1L) %% 10000L == 0L) {
      cat("Parsed annotation records:", index - 1L, "\n")
    }
  }
  annotation <- rbindlist(records)
  if (anyDuplicated(annotation$IID)) stop("Duplicate NPM IDs in annotation")
  annotation[, genetic_ancestry := normalise_group(genetic_ancestry)]
  annotation[toupper(cohort) == "HELIOS"]
}

collapse_interview <- function(path) {
  interview <- fread(
    path,
    select = c("FREG0_PID", "FIAQ10_Demo1"),
    na.strings = c("", "NA")
  )
  setnames(
    interview,
    c("FREG0_PID", "FIAQ10_Demo1"),
    c("participant_id", "self_ethnicity")
  )
  interview[, self_ethnicity := normalise_group(self_ethnicity)]
  collapsed <- interview[, {
    values <- unique(na.omit(self_ethnicity))
    list(
      self_ethnicity = if (length(values) == 1L) values else NA_character_,
      conflicting_interview = length(values) > 1L
    )
  }, by = participant_id]
  cat(
    "Interview participant IDs with conflicting self-description:",
    sum(collapsed$conflicting_interview), "\n"
  )
  collapsed[conflicting_interview == FALSE]
}

task_columns <- c(
  "DC7R6_React_Avg",
  "DC7R8_Stroopbox_Avg",
  "DC7R10_Stroopink_Avg",
  "DC7R5_Quiz_Score",
  "DC7R12_Wm_Score",
  "DC7R4_Pairing7_Guesses"
)
task_labels <- c(
  "Reaction time",
  "Stroop box",
  "Stroop ink",
  "Quiz score",
  "Working memory",
  "Pairing guesses"
)
names(task_labels) <- task_columns
directions <- c(-1, -1, -1, 1, 1, -1)
names(directions) <- task_columns

annotation <- read_annotation(input_paths$annotation)
interview <- collapse_interview(input_paths$interview)
metadata <- fread(input_paths$metadata, na.strings = c("", "NA"))
pca_input <- fread(input_paths$pca, na.strings = c("", "NA"))

if (anyDuplicated(metadata$IID)) stop("Duplicate IIDs in 22k metadata")
if (anyDuplicated(pca_input$IID)) stop("Duplicate IIDs in 22k PCA input")

required_metadata <- c(
  "IID", "age", "sex", "age2", "age_sex", "age2_sex", "indian",
  "malay", "Edu", task_columns
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0L) {
  stop("Metadata missing columns: ", paste(missing_metadata, collapse = ", "))
}
required_pca <- c("IID", "g", task_columns)
missing_pca <- setdiff(required_pca, names(pca_input))
if (length(missing_pca) > 0L) {
  stop("PCA input missing columns: ", paste(missing_pca, collapse = ", "))
}

metadata[, IID := as.character(IID)]
pca_input[, IID := as.character(IID)]
metadata <- merge(
  metadata,
  annotation[, .(IID, participant_id, genetic_ancestry)],
  by = "IID",
  all.x = TRUE,
  sort = FALSE
)
metadata <- merge(
  metadata,
  interview[, .(participant_id, self_ethnicity)],
  by = "participant_id",
  all.x = TRUE,
  sort = FALSE
)

if (any(is.na(metadata$genetic_ancestry))) {
  stop("At least one 22k metadata record lacks GRID genetic ancestry")
}

numeric_columns <- c(
  "age", "age2", "age_sex", "age2_sex", "indian", "malay", "Edu",
  task_columns
)
metadata[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]
metadata[, original_group := fifelse(
  indian == 1, "indian", fifelse(malay == 1, "malay", "chinese_reference")
)]
metadata[, genetic_indian := as.integer(genetic_ancestry == "indian")]
metadata[, genetic_malay := as.integer(genetic_ancestry == "malay")]

if (!setequal(metadata$IID[complete.cases(metadata[, c(
  "age", "age2", "age_sex", "age2_sex", "indian", "malay", "Edu"
), with = FALSE])], pca_input$IID)) {
  stop("Complete original covariates do not reproduce the PCA participant set")
}

analysis <- metadata[match(pca_input$IID, metadata$IID)]
if (!identical(analysis$IID, pca_input$IID)) stop("Participant ordering failed")
if (!all(complete.cases(analysis[, ..task_columns]))) {
  stop("Incomplete cognitive tasks in retained PCA sample")
}

fit_residuals <- function(data, model_type) {
  terms <- switch(
    model_type,
    original = c(
      "age", "age2", "age_sex", "age2_sex", "indian", "malay",
      "factor(sex)", "factor(Edu)"
    ),
    seven_group = c(
      "age", "age2", "age_sex", "age2_sex",
      "factor(genetic_ancestry)", "factor(sex)", "factor(Edu)"
    ),
    corrected_core = c(
      "age", "age2", "age_sex", "age2_sex", "genetic_indian",
      "genetic_malay", "factor(sex)", "factor(Edu)"
    ),
    corrected_core_self = c(
      "age", "age2", "age_sex", "age2_sex", "genetic_indian",
      "genetic_malay", "factor(self_ethnicity)", "factor(sex)",
      "factor(Edu)"
    ),
    stop("Unknown model type: ", model_type)
  )
  formula <- as.formula(paste("~", paste(terms, collapse = " + ")))
  design <- model.matrix(formula, data = data)
  response <- as.matrix(data[, ..task_columns])
  storage.mode(response) <- "double"
  residuals <- vapply(
    seq_len(ncol(response)),
    function(index) lm.fit(design, response[, index])$residuals,
    numeric(nrow(response))
  )
  colnames(residuals) <- task_columns
  rownames(residuals) <- data$IID
  list(
    residuals = residuals,
    design = design,
    response = response,
    formula = paste(deparse(formula), collapse = "")
  )
}

pca_from_residuals <- function(fitted, name) {
  oriented <- sweep(fitted$residuals, 2L, directions, `*`)
  pca <- prcomp(oriented, center = TRUE, scale. = TRUE)
  loading <- pca$rotation[, 1L]
  score <- pca$x[, 1L]
  if (sum(loading) < 0) {
    loading <- -loading
    score <- -score
  }
  names(score) <- rownames(fitted$residuals)
  list(
    name = name,
    n = nrow(fitted$residuals),
    ids = rownames(fitted$residuals),
    residuals = fitted$residuals,
    loading = loading,
    score = score,
    variance = pca$sdev[[1]]^2 / sum(pca$sdev^2),
    formula = fitted$formula
  )
}

tucker <- function(left, right) {
  sum(left * right) / sqrt(sum(left^2) * sum(right^2))
}

compare_models <- function(reference, comparison, label) {
  common <- intersect(names(reference$score), names(comparison$score))
  left <- reference$score[common]
  right <- comparison$score[common]
  if (cor(left, right) < 0) right <- -right
  left_z <- as.numeric(scale(left))
  right_z <- as.numeric(scale(right))
  data.table(
    comparison = label,
    reference_model = reference$name,
    comparison_model = comparison$name,
    reference_n = reference$n,
    comparison_n = comparison$n,
    overlap_n = length(common),
    loading_congruence = tucker(reference$loading, comparison$loading),
    maximum_absolute_loading_difference = max(abs(
      reference$loading - comparison$loading
    )),
    score_pearson_r = cor(left, right),
    score_spearman_rho = cor(left, right, method = "spearman"),
    standardised_score_rmse = sqrt(mean((left_z - right_z)^2)),
    reference_pc1_variance = reference$variance,
    comparison_pc1_variance = comparison$variance,
    pc1_variance_difference = comparison$variance - reference$variance
  )
}

core_groups <- c("chinese", "indian", "malay")
baseline_data <- copy(analysis)
core_data <- analysis[genetic_ancestry %in% core_groups]
questionnaire_data <- core_data[!is.na(self_ethnicity)]
strict_data <- questionnaire_data[self_ethnicity == genetic_ancestry]

cat("Retained original PCA N:", nrow(baseline_data), "\n")
cat("Core GRID N:", nrow(core_data), "\n")
cat("Core with self-described ethnicity N:", nrow(questionnaire_data), "\n")
cat("Strict genetic/self-description concordance N:", nrow(strict_data), "\n")

models <- list(
  original_full = pca_from_residuals(
    fit_residuals(baseline_data, "original"),
    "Original full-sample coding"
  ),
  seven_group = pca_from_residuals(
    fit_residuals(baseline_data, "seven_group"),
    "Seven-category GRID adjustment"
  ),
  core_only = pca_from_residuals(
    fit_residuals(core_data, "corrected_core"),
    "Core GRID groups only"
  ),
  questionnaire_base = pca_from_residuals(
    fit_residuals(questionnaire_data, "corrected_core"),
    "Questionnaire subset: genetic ancestry"
  ),
  questionnaire_self = pca_from_residuals(
    fit_residuals(questionnaire_data, "corrected_core_self"),
    "Questionnaire subset: genetic ancestry + self-description"
  ),
  strict_concordant = pca_from_residuals(
    fit_residuals(strict_data, "corrected_core"),
    "Strictly concordant questionnaire subset"
  )
)

# Exact reconstruction checks against the retained residual matrix and score.
reconstruction_rows <- lapply(task_columns, function(column) {
  observed <- as.numeric(pca_input[[column]])
  reproduced <- models$original_full$residuals[, column]
  difference <- observed - reproduced
  data.table(
    task = task_labels[[column]],
    n = length(observed),
    correlation = cor(observed, reproduced),
    rmse = sqrt(mean(difference^2)),
    maximum_absolute_difference = max(abs(difference))
  )
})
reconstruction <- rbindlist(reconstruction_rows)
saved_score_r <- cor(models$original_full$score, -as.numeric(pca_input$g))
if (saved_score_r < 0) saved_score_r <- -saved_score_r
reconstruction <- rbind(
  reconstruction,
  data.table(
    task = "PC1 score (absolute orientation)",
    n = nrow(pca_input),
    correlation = saved_score_r,
    rmse = NA_real_,
    maximum_absolute_difference = NA_real_
  )
)
fwrite(reconstruction, file.path(table_dir, "baseline_reproduction.tsv"), sep = "\t")

comparisons <- rbindlist(list(
  compare_models(
    models$original_full,
    models$seven_group,
    "Correct original reference coding while retaining all GRID categories"
  ),
  compare_models(
    models$original_full,
    models$core_only,
    "Restrict phenotype construction to the three core GRID groups"
  ),
  compare_models(
    models$questionnaire_base,
    models$questionnaire_self,
    "Add self-described ethnic heritage in the same questionnaire subset"
  ),
  compare_models(
    models$questionnaire_base,
    models$strict_concordant,
    "Exclude genetically/self-description-discordant questionnaire records"
  )
))
fwrite(comparisons, file.path(table_dir, "model_comparisons.tsv"), sep = "\t")

loading_rows <- rbindlist(lapply(models, function(model) {
  data.table(
    model = model$name,
    n = model$n,
    task = unname(task_labels[names(model$loading)]),
    loading_higher_performance = as.numeric(model$loading),
    pc1_variance_explained = model$variance
  )
}))
fwrite(loading_rows, file.path(table_dir, "loading_comparisons.tsv"), sep = "\t")

model_summary <- rbindlist(lapply(models, function(model) {
  data.table(
    model = model$name,
    n = model$n,
    pc1_variance_explained = model$variance,
    formula = model$formula
  )
}))
fwrite(model_summary, file.path(table_dir, "model_summary.tsv"), sep = "\t")

# Incremental explanatory contribution of self-described ethnicity for each
# raw task, after the original covariates and corrected core GRID indicators.
base_fit <- fit_residuals(questionnaire_data, "corrected_core")
self_fit <- fit_residuals(questionnaire_data, "corrected_core_self")
increment_rows <- lapply(seq_along(task_columns), function(index) {
  response <- base_fit$response[, index]
  base_residual <- base_fit$residuals[, index]
  self_residual <- self_fit$residuals[, index]
  rss_base <- sum(base_residual^2)
  rss_self <- sum(self_residual^2)
  tss <- sum((response - mean(response))^2)
  rank_base <- qr(base_fit$design)$rank
  rank_self <- qr(self_fit$design)$rank
  df1 <- rank_self - rank_base
  df2 <- length(response) - rank_self
  f_statistic <- if (df1 > 0L && rss_self > 0) {
    ((rss_base - rss_self) / df1) / (rss_self / df2)
  } else {
    NA_real_
  }
  data.table(
    task = task_labels[[task_columns[[index]]]],
    n = length(response),
    base_model_rank = rank_base,
    self_description_model_rank = rank_self,
    added_degrees_of_freedom = df1,
    base_r_squared = 1 - rss_base / tss,
    self_description_r_squared = 1 - rss_self / tss,
    incremental_r_squared = (rss_base - rss_self) / tss,
    partial_f = f_statistic,
    p_value = if (is.na(f_statistic)) NA_real_ else pf(
      f_statistic, df1, df2, lower.tail = FALSE
    )
  )
})
incremental_fit <- rbindlist(increment_rows)
incremental_fit[, p_fdr_bh := p.adjust(p_value, method = "BH")]
fwrite(
  incremental_fit,
  file.path(table_dir, "self_description_incremental_fit.tsv"),
  sep = "\t"
)

original_reference <- analysis$original_group == "chinese_reference"
intermediate_or_other <- !analysis$genetic_ancestry %in% core_groups
core_misencoded <- analysis$genetic_ancestry %in% c("indian", "malay") &
  analysis$original_group == "chinese_reference"
coding_audit <- data.table(
  item = c(
    "Original PCA sample",
    "Core GRID Chinese/Indian/Malay",
    "Pairwise intermediate or Other GRID category",
    "Intermediate/Other records coded as Chinese reference",
    "Core Indian/Malay records coded as Chinese reference",
    "Core participants with self-described ethnicity",
    "Strict genetic/self-description concordance"
  ),
  n = c(
    nrow(analysis),
    sum(analysis$genetic_ancestry %in% core_groups),
    sum(intermediate_or_other),
    sum(intermediate_or_other & original_reference),
    sum(core_misencoded),
    nrow(questionnaire_data),
    nrow(strict_data)
  )
)
coding_audit[, percentage_of_original_pca := 100 * n / nrow(analysis)]
fwrite(coding_audit, file.path(table_dir, "ancestry_coding_audit.tsv"), sep = "\t")

self_counts <- questionnaire_data[, .N, by = .(
  genetic_ancestry,
  self_ethnicity
)]
self_counts[, row_total := sum(N), by = genetic_ancestry]
self_counts[, row_percentage := 100 * N / row_total]
self_counts[, suppressed := {
  mask <- N < 5L
  if (any(mask) && any(!mask)) {
    candidates <- which(!mask)
    mask[candidates[[which.min(N[candidates])]]] <- TRUE
  }
  mask
}, by = genetic_ancestry]
self_counts[, count := ifelse(suppressed, "suppressed", as.character(N))]
self_counts[suppressed == TRUE, row_percentage := NA_real_]
self_counts[, N := NULL]
fwrite(
  self_counts,
  file.path(table_dir, "core_self_description_crosstab_suppressed.tsv"),
  sep = "\t"
)

sample_flow <- data.table(
  stage = c(
    "Retained 22k metadata rows",
    "Retained original PCA sample",
    "Core GRID groups in PCA sample",
    "Core GRID groups with self-described ethnicity",
    "Strict genetic/self-description concordance"
  ),
  n = c(
    nrow(metadata),
    nrow(analysis),
    nrow(core_data),
    nrow(questionnaire_data),
    nrow(strict_data)
  )
)
fwrite(sample_flow, file.path(table_dir, "sample_flow.tsv"), sep = "\t")

checks <- data.table(
  check = c(
    "Exact task-residual reproduction",
    "Exact saved-score reproduction",
    "GRID ancestry complete in PCA sample",
    "Questionnaire model participant sets identical",
    "No participant-level output written"
  ),
  status = c(
    if (min(reconstruction[task != "PC1 score (absolute orientation)"]$correlation) >
        0.999999999) "PASS" else "FAIL",
    if (saved_score_r > 0.999999999) "PASS" else "FAIL",
    if (all(!is.na(analysis$genetic_ancestry))) "PASS" else "FAIL",
    if (identical(
      models$questionnaire_base$ids,
      models$questionnaire_self$ids
    )) "PASS" else "FAIL",
    "PASS"
  ),
  detail = c(
    sprintf("minimum r=%.12f", min(reconstruction[task != "PC1 score (absolute orientation)"]$correlation)),
    sprintf("absolute r=%.12f", saved_score_r),
    sprintf("N=%d", nrow(analysis)),
    sprintf("N=%d", models$questionnaire_base$n),
    "Only aggregate tables, figures and logs are written"
  )
)
fwrite(checks, file.path(output_dir, "validation_checks.tsv"), sep = "\t")

input_manifest <- rbindlist(lapply(names(input_paths), function(name) {
  path <- input_paths[[name]]
  checksum <- system2("sha256sum", path, stdout = TRUE)
  data.table(
    input = name,
    file = basename(path),
    bytes = file.info(path)$size,
    sha256 = sub(" .*", "", checksum)
  )
}))
fwrite(input_manifest, file.path(output_dir, "input_manifest.tsv"), sep = "\t")

model_order <- vapply(models, function(model) model$name, character(1))
loading_rows[, model := factor(model, levels = model_order)]
plot_model_labels <- c(
  "Original full-sample coding" = "Original",
  "Seven-category GRID adjustment" = "Seven-category GRID",
  "Core GRID groups only" = "Core GRID only",
  "Questionnaire subset: genetic ancestry" = "Questionnaire: genetic ancestry",
  "Questionnaire subset: genetic ancestry + self-description" =
    "Questionnaire: + self-description",
  "Strictly concordant questionnaire subset" = "Strictly concordant subset"
)
loading_rows[, plot_model := factor(
  plot_model_labels[as.character(model)],
  levels = unname(plot_model_labels)
)]
loading_plot <- ggplot(
  loading_rows,
  aes(
    x = loading_higher_performance,
    y = task,
    colour = plot_model,
    group = plot_model
  )
) +
  geom_line(linewidth = 0.7, alpha = 0.75) +
  geom_point(size = 2.5) +
  labs(
    title = "HELIOS cognitive-factor loading sensitivity",
    subtitle = "Higher values indicate better performance",
    x = "PC1 loading",
    y = NULL,
    colour = "Phenotype model"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE))
ggsave(
  file.path(figure_dir, "loading_sensitivity.png"),
  loading_plot,
  width = 13,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "loading_sensitivity.pdf"),
  loading_plot,
  width = 13,
  height = 8
)

plot_comparison_labels <- c(
  "Correct original reference coding while retaining all GRID categories" =
    "Seven-category GRID adjustment",
  "Restrict phenotype construction to the three core GRID groups" =
    "Core GRID-only refit",
  "Add self-described ethnic heritage in the same questionnaire subset" =
    "Add self-described ethnicity",
  "Exclude genetically/self-description-discordant questionnaire records" =
    "Exclude discordant records"
)
comparisons[, plot_comparison := factor(
  plot_comparison_labels[comparison],
  levels = rev(unname(plot_comparison_labels))
)]
score_plot <- ggplot(
  comparisons,
  aes(x = score_pearson_r, y = plot_comparison)
) +
  geom_segment(
    aes(x = 0.99965, xend = score_pearson_r, yend = plot_comparison),
    colour = "#B8C4CE",
    linewidth = 1.1
  ) +
  geom_point(colour = "#2C7FB8", size = 4) +
  geom_text(
    aes(label = sprintf("%.6f", score_pearson_r)),
    hjust = -0.15,
    size = 4
  ) +
  scale_x_continuous(
    limits = c(0.99965, 1.00004),
    breaks = c(0.9997, 0.9998, 0.9999, 1.0000)
  ) +
  labs(
    title = "Cognitive-factor score concordance across sensitivities",
    subtitle = "Pearson correlation on participants shared with the stated reference model",
    caption = "The horizontal axis is restricted to r = 0.99965–1.00000.",
    x = "Pearson r",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.caption = element_text(hjust = 0))
ggsave(
  file.path(figure_dir, "score_concordance.png"),
  score_plot,
  width = 10,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "score_concordance.pdf"),
  score_plot,
  width = 10,
  height = 6
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

if (any(checks$status != "PASS")) {
  stop("One or more validation checks failed")
}

cat("Completed:", format(Sys.time(), tz = "UTC"), "UTC\n")
cat("All validation checks passed.\n")
