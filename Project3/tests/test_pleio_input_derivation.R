#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript test_pleio_input_derivation.R /path/to/Project3", call. = FALSE)
}
project_dir <- normalizePath(args[1L], mustWork = TRUE)
derive_script <- file.path(project_dir, "scripts", "pleio", "derive_three_trait_inputs.R")

tmp <- tempfile("project3_pleio_test_")
dir.create(tmp, recursive = TRUE)
input_dir <- file.path(tmp, "four_trait")
output_dir <- file.path(tmp, "three_trait")
dir.create(input_dir)

write_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

traits <- c("T1", "T2", "T3", "T4")
metain <- data.frame(SNP = c("rs1", "rs2"), check.names = FALSE)
for (i in seq_along(traits)) {
  metain[[paste0(traits[i], "_beta")]] <- c(i, -i) / 10
  metain[[paste0(traits[i], "_se")]] <- c(0.01, 0.02)
}
sg <- matrix(seq_len(16) / 100, nrow = 4, dimnames = list(traits, traits))
sg <- (sg + t(sg)) / 2
ce <- diag(4)
colnames(ce) <- rownames(ce) <- traits

write_gz(metain, file.path(input_dir, "metain.txt.gz"))
write_gz(as.data.frame(sg, check.names = FALSE), file.path(input_dir, "sg.txt.gz"))
write_gz(as.data.frame(ce, check.names = FALSE), file.path(input_dir, "ce.txt.gz"))

model_config <- file.path(tmp, "models.tsv")
write.table(
  data.frame(model = rep("TEST_MODEL", 3), trait = c("T1", "T3", "T4")),
  model_config,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

run_output <- system2(
  "Rscript",
  c(
    derive_script,
    "--input-dir", input_dir,
    "--model-config", model_config,
    "--out-root", output_dir
  ),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(run_output, "status")) && attr(run_output, "status") != 0L) {
  stop(paste(run_output, collapse = "\n"), call. = FALSE)
}

model_dir <- file.path(output_dir, "TEST_MODEL")
derived_meta <- read.delim(gzfile(file.path(model_dir, "metain.txt.gz")), check.names = FALSE)
derived_sg <- read.delim(gzfile(file.path(model_dir, "sg.txt.gz")), check.names = FALSE)
derived_ce <- read.delim(gzfile(file.path(model_dir, "ce.txt.gz")), check.names = FALSE)

expected_columns <- c(
  "SNP",
  "T1_beta", "T1_se",
  "T3_beta", "T3_se",
  "T4_beta", "T4_se"
)
stopifnot(
  identical(names(derived_meta), expected_columns),
  identical(derived_meta$SNP, c("rs1", "rs2")),
  identical(names(derived_sg), c("T1", "T3", "T4")),
  identical(names(derived_ce), c("T1", "T3", "T4")),
  nrow(derived_sg) == 3L,
  nrow(derived_ce) == 3L
)

unlink(tmp, recursive = TRUE)
cat("PLEIO three-trait input derivation test passed.\n")
