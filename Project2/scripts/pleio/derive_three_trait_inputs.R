#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

usage <- function(status = 0L) {
  cat(
    "Usage:\n",
    "  Rscript derive_three_trait_inputs.R --input-dir FOUR_TRAIT_DIR \\\n",
    "    --model-config MODELS.tsv --out-root DIR [--validate-only]\n\n",
    "The model config has two tab-separated columns: model and trait.\n",
    "Trait row order determines the PLEIO column/matrix order.\n",
    sep = ""
  )
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || "--help" %in% args) usage(0L)

flag_value <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) != 1L || hit == length(args) || startsWith(args[hit + 1L], "--")) {
    stop("Missing or invalid option: ", flag, call. = FALSE)
  }
  args[hit + 1L]
}

input_dir <- normalizePath(flag_value("--input-dir"), mustWork = TRUE)
model_config <- normalizePath(flag_value("--model-config"), mustWork = TRUE)
out_root <- flag_value("--out-root")
validate_only <- "--validate-only" %in% args

required_files <- file.path(input_dir, c("metain.txt.gz", "sg.txt.gz", "ce.txt.gz"))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing four-trait PLEIO inputs:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
}

models_long <- read.delim(model_config, sep = "\t", header = TRUE, check.names = FALSE)
if (!identical(names(models_long), c("model", "trait"))) {
  stop("Model config must have exactly two columns: model and trait.", call. = FALSE)
}
if (nrow(models_long) == 0L || any(!nzchar(models_long$model)) || any(!nzchar(models_long$trait))) {
  stop("Model config contains blank fields.", call. = FALSE)
}
model_order <- unique(models_long$model)
models <- lapply(model_order, function(model) models_long$trait[models_long$model == model])
names(models) <- model_order
if (any(vapply(models, function(x) length(x) != 3L || anyDuplicated(x), logical(1)))) {
  stop("Each model must contain exactly three unique traits.", call. = FALSE)
}

read_gz_table <- function(path, nrows = -1L) {
  read.delim(gzfile(path), header = TRUE, sep = "\t", check.names = FALSE, nrows = nrows)
}
matrix_from_pleio_file <- function(df, label) {
  m <- as.matrix(df)
  storage.mode(m) <- "numeric"
  rownames(m) <- colnames(df)
  if (nrow(m) != ncol(m) || !identical(rownames(m), colnames(m))) {
    stop(label, " must be a square matrix whose row order matches its columns.", call. = FALSE)
  }
  m
}

sg <- matrix_from_pleio_file(read_gz_table(file.path(input_dir, "sg.txt.gz")), "sg.txt.gz")
ce <- matrix_from_pleio_file(read_gz_table(file.path(input_dir, "ce.txt.gz")), "ce.txt.gz")
metain_header <- names(read_gz_table(file.path(input_dir, "metain.txt.gz"), nrows = 1L))

for (model in names(models)) {
  traits <- models[[model]]
  needed <- c("SNP", as.vector(rbind(paste0(traits, "_beta"), paste0(traits, "_se"))))
  missing_meta <- setdiff(needed, metain_header)
  missing_sg <- setdiff(traits, colnames(sg))
  missing_ce <- setdiff(traits, colnames(ce))
  if (length(missing_meta) || length(missing_sg) || length(missing_ce)) {
    stop(
      "Missing fields for ", model, ": ",
      paste(unique(c(missing_meta, missing_sg, missing_ce)), collapse = ", "),
      call. = FALSE
    )
  }
  cat("Validated", model, ":", paste(traits, collapse = ", "), "\n")
}
if (validate_only) quit(save = "no", status = 0L)

cat("Reading four-trait metain file; this may require several GB of RAM.\n")
metain <- read_gz_table(file.path(input_dir, "metain.txt.gz"))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write_gz_table <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

for (model in names(models)) {
  traits <- models[[model]]
  out_dir <- file.path(out_root, model)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  metain_cols <- c("SNP", as.vector(rbind(paste0(traits, "_beta"), paste0(traits, "_se"))))

  write_gz_table(metain[, metain_cols, drop = FALSE], file.path(out_dir, "metain.txt.gz"))
  write_gz_table(
    as.data.frame(sg[traits, traits, drop = FALSE], check.names = FALSE),
    file.path(out_dir, "sg.txt.gz")
  )
  write_gz_table(
    as.data.frame(ce[traits, traits, drop = FALSE], check.names = FALSE),
    file.path(out_dir, "ce.txt.gz")
  )
  write.table(
    data.frame(model = model, trait = traits),
    file.path(out_dir, "model_traits.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat("Wrote", model, "to", out_dir, "(", nrow(metain), "variants )\n")
}
