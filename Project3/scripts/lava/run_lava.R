#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

usage <- function(status = 0L) {
  cat(
    "Usage:\n",
    "  Rscript run_lava.R --input-info FILE --sample-overlap FILE \\\n",
    "    --locus-file FILE --ref-prefix PREFIX --phenos ID1,ID2,... \\\n",
    "    --sumstats-dir DIR --out-dir DIR --out-prefix NAME [options]\n\n",
    "Options:\n",
    "  --univ-threshold P  Univariate screening P threshold (default: 2e-5)\n",
    "  --validate-only     Validate all inputs without running LAVA\n",
    "  --help              Show this help\n",
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

input_info_file <- normalizePath(flag_value("--input-info"), mustWork = TRUE)
sample_overlap_file <- normalizePath(flag_value("--sample-overlap"), mustWork = TRUE)
locus_file <- normalizePath(flag_value("--locus-file"), mustWork = TRUE)
ref_prefix <- flag_value("--ref-prefix")
phenos <- trimws(strsplit(flag_value("--phenos"), ",", fixed = TRUE)[[1L]])
sumstats_dir <- normalizePath(flag_value("--sumstats-dir"), mustWork = TRUE)
out_dir <- flag_value("--out-dir")
out_prefix <- flag_value("--out-prefix")
univ_threshold <- suppressWarnings(as.numeric(flag_value("--univ-threshold", FALSE, "2e-5")))
validate_only <- "--validate-only" %in% args

if (!length(phenos) || any(!nzchar(phenos)) || anyDuplicated(phenos)) {
  stop("--phenos must contain unique comma-separated IDs.", call. = FALSE)
}
if (is.na(univ_threshold) || univ_threshold <= 0 || univ_threshold > 1) {
  stop("--univ-threshold must be in (0,1].", call. = FALSE)
}

input_info <- read.delim(input_info_file, sep = "\t", header = TRUE, check.names = FALSE)
required_info <- c("phenotype", "cases", "controls", "filename")
missing_info <- setdiff(required_info, names(input_info))
if (length(missing_info)) {
  stop("Input-info file is missing: ", paste(missing_info, collapse = ", "), call. = FALSE)
}
input_info <- input_info[nzchar(input_info$phenotype), required_info, drop = FALSE]
if (anyDuplicated(input_info$phenotype)) stop("Input-info phenotype IDs must be unique.", call. = FALSE)
if (!all(phenos %in% input_info$phenotype)) {
  stop("Phenotypes absent from input-info: ", paste(setdiff(phenos, input_info$phenotype), collapse = ", "), call. = FALSE)
}

resolve_input <- function(path) {
  if (grepl("^/", path)) path else file.path(sumstats_dir, path)
}
input_info$filename <- vapply(input_info$filename, resolve_input, character(1))
missing_sumstats <- input_info$filename[!file.exists(input_info$filename)]
if (length(missing_sumstats)) {
  stop("Missing LAVA summary-statistics files:\n", paste(missing_sumstats, collapse = "\n"), call. = FALSE)
}

overlap <- read.table(
  sample_overlap_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
if (!all(phenos %in% rownames(overlap)) || !all(phenos %in% colnames(overlap))) {
  stop("Sample-overlap row and column names must include every --phenos ID.", call. = FALSE)
}
overlap_sub <- as.matrix(overlap[phenos, phenos, drop = FALSE])
storage.mode(overlap_sub) <- "numeric"
if (any(!is.finite(overlap_sub)) || !isTRUE(all.equal(overlap_sub, t(overlap_sub)))) {
  stop("Selected sample-overlap matrix must be finite and symmetric.", call. = FALSE)
}

locus_header <- names(read.table(locus_file, header = TRUE, nrows = 1, check.names = FALSE))
if (!all(c("LOC", "CHR", "START", "STOP") %in% locus_header)) {
  stop("Locus file must contain LOC, CHR, START, and STOP.", call. = FALSE)
}
missing_ref <- paste0(ref_prefix, c(".bed", ".bim", ".fam"))
missing_ref <- missing_ref[!file.exists(missing_ref)]
if (length(missing_ref)) {
  stop("Missing PLINK reference files:\n", paste(missing_ref, collapse = "\n"), call. = FALSE)
}

cat("Validated phenotypes:", paste(phenos, collapse = ", "), "\n")
cat("Univariate threshold:", format(univ_threshold, scientific = TRUE), "\n")
cat("Reference prefix:", ref_prefix, "\n")
if (validate_only) quit(save = "no", status = 0L)

if (!requireNamespace("LAVA", quietly = TRUE)) {
  stop("R package LAVA is required.", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
resolved_input_file <- file.path(out_dir, paste0(out_prefix, "_resolved_input.info.tsv"))
write.table(
  input_info,
  resolved_input_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

cat("Reading loci from:", locus_file, "\n")
loci <- LAVA::read.loci(locus_file)
n_loci <- nrow(loci)
cat("Number of loci:", n_loci, "\n")

input <- LAVA::process.input(
  input.info = resolved_input_file,
  sample.overlap = sample_overlap_file,
  ref.prefix = ref_prefix,
  phenos = phenos
)

progress <- unique(ceiling(stats::quantile(seq_len(n_loci), seq(0.05, 1, 0.05))))
univ <- vector("list", n_loci)
bivar <- vector("list", n_loci)

for (i in seq_len(n_loci)) {
  if (i %in% progress) {
    cat("Locus", i, "of", n_loci, "\n")
    flush.console()
  }
  locus <- LAVA::process.locus(loci[i, ], input)
  if (is.null(locus)) next

  loc_info <- data.frame(
    locus = locus$id,
    chr = locus$chr,
    start = locus$start,
    stop = locus$stop,
    n.snps = locus$n.snps,
    n.pcs = locus$K
  )
  loc_out <- LAVA::run.univ.bivar(locus, univ.thresh = univ_threshold)
  univ[[i]] <- cbind(loc_info, loc_out$univ)
  if (!is.null(loc_out$bivar)) bivar[[i]] <- cbind(loc_info, loc_out$bivar)
}

bind_unequal <- function(items) {
  items <- items[!vapply(items, is.null, logical(1))]
  if (!length(items)) return(NULL)
  all_names <- Reduce(union, lapply(items, names))
  items <- lapply(items, function(x) {
    for (field in setdiff(all_names, names(x))) x[[field]] <- NA
    x[all_names]
  })
  do.call(rbind, items)
}

univ_out <- bind_unequal(univ)
bivar_out <- bind_unequal(bivar)
univ_file <- file.path(out_dir, paste0(out_prefix, ".univ.lava"))
bivar_file <- file.path(out_dir, paste0(out_prefix, ".bivar.lava"))

if (!is.null(univ_out)) {
  write.table(univ_out, univ_file, sep = "\t", row.names = FALSE, quote = FALSE)
}
if (!is.null(bivar_out)) {
  write.table(bivar_out, bivar_file, sep = "\t", row.names = FALSE, quote = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0(out_prefix, "_sessionInfo.txt")))

cat("Univariate output:", univ_file, "\n")
cat("Bivariate output:", bivar_file, "\n")
