#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(LAVA))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: run_four_trait_lava.R",
      "<input.info.tsv> <sample_overlap.tsv> <reference_prefix>",
      "<output_prefix>"
    ),
    call. = FALSE
  )
}

input_info_file <- normalizePath(args[[1L]], mustWork = TRUE)
sample_overlap_file <- normalizePath(args[[2L]], mustWork = TRUE)
reference_prefix <- args[[3L]]
output_prefix <- args[[4L]]

locus_file <- "/home/ec2-user/software/COLOC-reporter/EAS_LAVA.locfile"
phenotypes <- c("MDD_EAS", "Edu_EAS", "SCZ_EAS", "HEL10k_EAS")
univariate_threshold <- 2e-5

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

cat("LAVA version:", as.character(packageVersion("LAVA")), "\n")
cat("Input info:", input_info_file, "\n")
cat("Sample overlap:", sample_overlap_file, "\n")
cat("Reference prefix:", reference_prefix, "\n")
cat("Locus file:", locus_file, "\n")
cat("Univariate threshold:", univariate_threshold, "\n")

loci <- read.loci(locus_file)
cat("Number of loci:", nrow(loci), "\n")

input <- process.input(
  input.info.file = input_info_file,
  sample.overlap.file = sample_overlap_file,
  ref.prefix = reference_prefix,
  phenos = phenotypes
)

univariate_results <- vector("list", nrow(loci))
bivariate_results <- vector("list", nrow(loci))
progress <- unique(ceiling(quantile(seq_len(nrow(loci)), seq(0.05, 1, 0.05))))

for (index in seq_len(nrow(loci))) {
  if (index %in% progress) {
    cat("Progress:", index, "of", nrow(loci), "loci\n")
    flush.console()
  }

  locus <- process.locus(loci[index, ], input)
  if (is.null(locus)) {
    next
  }

  locus_information <- data.frame(
    locus = locus$id,
    chr = locus$chr,
    start = locus$start,
    stop = locus$stop,
    n.snps = locus$n.snps,
    n.pcs = locus$K
  )
  result <- run.univ.bivar(locus, univ.thresh = univariate_threshold)
  univariate_results[[index]] <- cbind(locus_information, result$univ)
  if (!is.null(result$bivar)) {
    bivariate_results[[index]] <- cbind(locus_information, result$bivar)
  }
}

bind_complete <- function(results) {
  results <- results[!vapply(results, is.null, logical(1L))]
  if (length(results) == 0L) {
    return(NULL)
  }
  all_columns <- Reduce(union, lapply(results, colnames))
  results <- lapply(results, function(result) {
    missing_columns <- setdiff(all_columns, colnames(result))
    for (column in missing_columns) {
      result[[column]] <- NA
    }
    result[all_columns]
  })
  do.call(rbind, results)
}

univariate_output <- bind_complete(univariate_results)
bivariate_output <- bind_complete(bivariate_results)

write.table(
  univariate_output,
  paste0(output_prefix, ".univ.lava"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  bivariate_output,
  paste0(output_prefix, ".bivar.lava"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
writeLines(capture.output(sessionInfo()), paste0(output_prefix, ".sessionInfo.txt"))

cat("Univariate rows:", nrow(univariate_output), "\n")
cat("Bivariate rows:", nrow(bivariate_output), "\n")
cat("Completed output prefix:", output_prefix, "\n")
