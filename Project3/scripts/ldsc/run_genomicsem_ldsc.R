#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

usage <- function(status = 0L) {
  cat(
    "Usage:\n",
    "  Rscript run_genomicsem_ldsc.R \\\n",
    "    --trait-config TRAITS.tsv --sumstats-dir DIR \\\n",
    "    --hm3 FILE --ld DIR --out-dir DIR --out-prefix NAME [options]\n\n",
    "Required trait-config columns:\n",
    "  file, trait, N, sample_prev, population_prev\n\n",
    "Options:\n",
    "  --wld DIR          Regression-weight LD scores (default: --ld)\n",
    "  --skip-munge       Reuse <out-dir>/<trait>.sumstats.gz files\n",
    "  --no-heatmap       Do not write the genetic-correlation heatmap\n",
    "  --validate-only    Validate paths/configuration without running LDSC\n",
    "  --help             Show this help\n",
    sep = ""
  )
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || "--help" %in% args) usage(0L)

flag_value <- function(flag, required = TRUE, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0L) {
    if (required) stop("Missing required option: ", flag, call. = FALSE)
    return(default)
  }
  if (length(hit) > 1L || hit == length(args) || startsWith(args[hit + 1L], "--")) {
    stop("Option must occur once and have a value: ", flag, call. = FALSE)
  }
  args[hit + 1L]
}

trait_config <- normalizePath(flag_value("--trait-config"), mustWork = TRUE)
sumstats_dir <- normalizePath(flag_value("--sumstats-dir"), mustWork = TRUE)
hm3 <- normalizePath(flag_value("--hm3"), mustWork = TRUE)
ld_folder <- normalizePath(flag_value("--ld"), mustWork = TRUE)
wld_folder <- normalizePath(flag_value("--wld", FALSE, ld_folder), mustWork = TRUE)
out_dir <- flag_value("--out-dir")
out_prefix <- flag_value("--out-prefix")
skip_munge <- "--skip-munge" %in% args
make_heatmap <- !("--no-heatmap" %in% args)
validate_only <- "--validate-only" %in% args

cfg <- read.delim(trait_config, sep = "\t", header = TRUE, check.names = FALSE)
required_columns <- c("file", "trait", "N", "sample_prev", "population_prev")
missing_columns <- setdiff(required_columns, names(cfg))
if (length(missing_columns)) {
  stop("Trait configuration is missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (nrow(cfg) < 2L) stop("At least two traits are required.", call. = FALSE)
if (anyDuplicated(cfg$trait)) stop("Trait names must be unique.", call. = FALSE)
if (any(!nzchar(cfg$file)) || any(!nzchar(cfg$trait))) {
  stop("The file and trait fields cannot be blank.", call. = FALSE)
}

resolve_input <- function(path) {
  if (grepl("^/", path)) path else file.path(sumstats_dir, path)
}
sumstat_files <- vapply(cfg$file, resolve_input, character(1))
missing_inputs <- sumstat_files[!file.exists(sumstat_files)]
if (length(missing_inputs)) {
  stop("Missing summary-statistics files:\n", paste(missing_inputs, collapse = "\n"), call. = FALSE)
}

as_numeric_na <- function(x, field) {
  x[x == "" | toupper(x) == "NA"] <- NA_character_
  y <- suppressWarnings(as.numeric(x))
  if (any(is.na(y) & !is.na(x))) stop("Non-numeric value in ", field, call. = FALSE)
  y
}

N_vec <- as_numeric_na(as.character(cfg$N), "N")
sample_prev <- as_numeric_na(as.character(cfg$sample_prev), "sample_prev")
population_prev <- as_numeric_na(as.character(cfg$population_prev), "population_prev")

if (any(sample_prev <= 0 | sample_prev >= 1, na.rm = TRUE) ||
    any(population_prev <= 0 | population_prev >= 1, na.rm = TRUE)) {
  stop("Prevalences must be between zero and one.", call. = FALSE)
}

cat("Validated", nrow(cfg), "traits:", paste(cfg$trait, collapse = ", "), "\n")
cat("HapMap3 list:", hm3, "\n")
cat("LD scores:", ld_folder, "\n")
cat("Weight LD scores:", wld_folder, "\n")
if (validate_only) quit(save = "no", status = 0L)

if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("R package GenomicSEM is required.", call. = FALSE)
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("R package Matrix is required.", call. = FALSE)
}
if (make_heatmap && !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("R package ggplot2 is required unless --no-heatmap is used.", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
trait_names <- cfg$trait
munged_files <- file.path(out_dir, paste0(trait_names, ".sumstats.gz"))

if (!skip_munge) {
  old_wd <- getwd()
  setwd(out_dir)
  message("Munging summary statistics in ", out_dir)
  tryCatch(
    GenomicSEM::munge(
      files = sumstat_files,
      hm3 = hm3,
      trait.names = trait_names,
      N = N_vec
    ),
    finally = setwd(old_wd)
  )
}

missing_munged <- munged_files[!file.exists(munged_files)]
if (length(missing_munged)) {
  stop("Missing munged files:\n", paste(missing_munged, collapse = "\n"), call. = FALSE)
}

message("Running GenomicSEM::ldsc()")
ldsc_result <- GenomicSEM::ldsc(
  traits = munged_files,
  sample.prev = sample_prev,
  population.prev = population_prev,
  ld = ld_folder,
  wld = wld_folder,
  trait.names = trait_names,
  stand = TRUE
)

saveRDS(ldsc_result, file.path(out_dir, paste0(out_prefix, "_LDSC_INT.rds")))
LDSC_INT_save <- ldsc_result
save(LDSC_INT_save, file = file.path(out_dir, paste0(out_prefix, "_LDSC_INT.RData")))

if (!("S" %in% names(ldsc_result))) stop("LDSC result has no S covariance matrix.", call. = FALSE)
S <- as.matrix(ldsc_result$S)
if (nrow(S) != length(trait_names)) stop("Unexpected covariance-matrix dimensions.", call. = FALSE)
write.csv(S, file.path(out_dir, paste0(out_prefix, "_genetic_covariance.csv")), quote = FALSE)

S_pd <- as.matrix(Matrix::nearPD(S)$mat)
rg <- stats::cov2cor(S_pd)
rg_file <- file.path(out_dir, paste0(out_prefix, "_genetic_correlations.csv"))
write.csv(rg, rg_file, quote = FALSE)

h2 <- data.frame(
  trait = trait_names,
  SNP_h2 = diag(S_pd),
  SNP_h2_SE = NA_real_,
  stringsAsFactors = FALSE
)
write.csv(
  h2,
  file.path(out_dir, paste0(out_prefix, "_heritability_summary.csv")),
  row.names = FALSE,
  quote = FALSE
)

if (make_heatmap) {
  heatmap_df <- as.data.frame(as.table(rg))
  names(heatmap_df) <- c("Trait1", "Trait2", "rg")
  heatmap_df$Trait1 <- factor(heatmap_df$Trait1, levels = trait_names)
  heatmap_df$Trait2 <- factor(heatmap_df$Trait2, levels = trait_names)
  p <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = Trait1, y = Trait2, fill = rg)) +
    ggplot2::geom_tile(color = "grey80") +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "rg") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = "Genetic correlations (LDSC / GenomicSEM)")
  ggplot2::ggsave(
    file.path(out_dir, paste0(out_prefix, "_rg_heatmap.png")),
    p,
    width = 2000 / 300,
    height = 2000 / 300,
    dpi = 300
  )
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0(out_prefix, "_sessionInfo.txt")))
message("Completed LDSC workflow: ", out_prefix)
