#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# Post-hoc directional classification used for the manuscript.
#
# Independent significant PLEIO SNPs are classified from allele-aligned
# marginal Z scores. Opposing psychiatric and cognitive/educational directions
# are Concordant, matching directions are Discordant, and within-domain sign
# heterogeneity is Mixed. SNP classes are collapsed to FUMA loci with the
# hierarchy Mixed > Dual > Concordant/Discordant > Unassigned.

parse_args <- function(x) {
  if (length(x) == 0L || any(x %in% c("-h", "--help"))) {
    return(list(help = TRUE))
  }
  if (length(x) %% 2L != 0L) {
    stop("Arguments must be supplied as --name value pairs. Use --help for usage.")
  }
  keys <- sub("^--", "", x[seq(1L, length(x), by = 2L)])
  vals <- x[seq(2L, length(x), by = 2L)]
  as.list(setNames(vals, keys))
}

usage <- function() {
  cat(
    "\nUsage:\n",
    "Rscript classify_pleio_locus_direction.R \\\n",
    "  --ancestry EUR \\\n",
    "  --model FOURTRAIT \\\n",
    "  --pleio-dir /path/to/model/FUMA/PLEIO \\\n",
    "  --single-trait-root /path/to/ancestry/FUMA \\\n",
    "  --raw-file /path/to/raw_merged_PLEIO.tsv[.gz] \\\n",
    "  --out-dir /path/to/output \\\n",
    "  [--ea-folder EA] [--cf-folder CF] [--zthr 1.96] \\\n",
    "  [--prefix EUR_FOURTRAIT]\n\n",
    "Accepted models: FOURTRAIT, MDD_3TRAIT, SCZ_3TRAIT\n",
    "The single-trait root must contain the model-relevant EA, CF, MDD, ",
    "and/or SCZ FUMA directories. --ea-folder and --cf-folder only map ",
    "historical source-directory names to the canonical EA and CF labels.\n",
    sep = ""
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0L)
}

required_args <- c(
  "ancestry", "model", "pleio-dir", "single-trait-root", "raw-file", "out-dir"
)
missing_args <- required_args[!required_args %in% names(args)]
if (length(missing_args) > 0L) {
  usage()
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

ancestry <- toupper(args$ancestry)
model_raw <- toupper(gsub("-", "_", args$model))
model <- switch(
  model_raw,
  "FOUR" = "FOURTRAIT",
  "4TRAIT" = "FOURTRAIT",
  "FOUR_TRAIT" = "FOURTRAIT",
  "FOURTRAIT" = "FOURTRAIT",
  "MDD" = "MDD_3TRAIT",
  "MDD3" = "MDD_3TRAIT",
  "MDD_3TRAIT" = "MDD_3TRAIT",
  "SCZ" = "SCZ_3TRAIT",
  "SCZ3" = "SCZ_3TRAIT",
  "SCZ_3TRAIT" = "SCZ_3TRAIT",
  stop(
    "Unknown --model: ", args$model,
    ". Accepted: FOURTRAIT, MDD_3TRAIT, SCZ_3TRAIT"
  )
)

pleio_dir <- normalizePath(args[["pleio-dir"]], mustWork = TRUE)
single_trait_root <- normalizePath(args[["single-trait-root"]], mustWork = TRUE)
raw_file <- normalizePath(args[["raw-file"]], mustWork = TRUE)
out_dir <- args[["out-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

ea_folder <- if ("ea-folder" %in% names(args)) args[["ea-folder"]] else "EA"
cf_folder <- if ("cf-folder" %in% names(args)) args[["cf-folder"]] else "CF"
zthr <- if ("zthr" %in% names(args)) as.numeric(args$zthr) else 1.96
if (!is.finite(zthr) || zthr < 0) stop("--zthr must be a non-negative number.")
prefix <- if ("prefix" %in% names(args)) args$prefix else paste(ancestry, model, sep = "_")

all_traits <- c("EA", "CF", "MDD", "SCZ")
model_traits <- switch(
  model,
  "FOURTRAIT" = c("EA", "CF", "MDD", "SCZ"),
  "MDD_3TRAIT" = c("EA", "CF", "MDD"),
  "SCZ_3TRAIT" = c("EA", "CF", "SCZ")
)

message("Ancestry: ", ancestry)
message("Model: ", model)
message("PLEIO FUMA directory: ", pleio_dir)
message("Single-trait FUMA root: ", single_trait_root)
message("Raw merged PLEIO file: ", raw_file)
message("Output directory: ", out_dir)
message("Direction threshold |Z| >= ", zthr)
message("Model comparators: ", paste(model_traits, collapse = ", "))

find_unique_file <- function(dir, pattern, label, allow_zero = FALSE) {
  if (!dir.exists(dir)) {
    if (allow_zero) {
      message("Directory absent for ", label, ": ", dir, "; treating as 0 loci.")
      return(NA_character_)
    }
    stop("Directory not found for ", label, ": ", dir)
  }
  files <- list.files(dir, recursive = TRUE, full.names = TRUE)
  hits <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
  if (length(hits) == 0L) {
    if (allow_zero) {
      message("No file matched for ", label, " in ", dir, "; treating as 0 loci.")
      return(NA_character_)
    }
    stop("No file matched for ", label, " in ", dir, " using pattern ", pattern)
  }
  if (length(hits) > 1L) {
    stop(
      "Multiple files matched for ", label, ":\n", paste(hits, collapse = "\n"),
      "\nPass a directory containing only the selected FUMA run."
    )
  }
  normalizePath(hits[1L], mustWork = TRUE)
}

find_col <- function(dt, candidates, required = TRUE) {
  nms <- names(dt)
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(nms[idx[1L]])
  if (required) {
    stop(
      "Could not find any of these columns: ", paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(nms, collapse = ", ")
    )
  }
  NA_character_
}

normalise_chr <- function(x) {
  x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
  x[x %chin% c("23", "X")] <- "X"
  x[x %chin% c("24", "Y")] <- "Y"
  x
}

read_fuma_loci <- function(trait_dir, trait_name, allow_zero = FALSE) {
  f <- find_unique_file(
    trait_dir,
    "^GenomicRiskLoci.*\\.txt$",
    paste0(trait_name, " GenomicRiskLoci"),
    allow_zero = allow_zero
  )
  if (is.na(f)) return(data.table())

  message("Reading ", trait_name, " loci: ", f)
  dt <- fread(f)
  if (nrow(dt) == 0L) return(data.table())

  chr_col <- find_col(dt, c("chr", "CHR", "chromosome"))
  start_col <- find_col(dt, c("start", "Start", "START", "posStart", "locusStart"))
  end_col <- find_col(dt, c("end", "End", "END", "stop", "STOP", "posEnd", "locusEnd"))
  locus_col <- find_col(
    dt, c("GenomicLocus", "genomicLocus", "locus", "Locus"), required = FALSE
  )

  dt[, trait := trait_name]
  dt[, chr_norm := normalise_chr(get(chr_col))]
  dt[, start_int := suppressWarnings(as.integer(get(start_col)))]
  dt[, end_int := suppressWarnings(as.integer(get(end_col)))]
  if (anyNA(dt$chr_norm) || anyNA(dt$start_int) || anyNA(dt$end_int)) {
    stop("Missing or invalid coordinates in ", f)
  }
  if (any(dt$start_int > dt$end_int)) stop("At least one locus has start > end in ", f)

  if (!is.na(locus_col)) {
    dt[, fuma_locus_id := as.character(get(locus_col))]
  } else {
    dt[, fuma_locus_id := paste0(chr_norm, ":", start_int, "-", end_int)]
  }
  dt[, locus_key := paste0(trait_name, "_", fuma_locus_id)]
  dt
}

overlap_with_trait <- function(pleio_loci, trait_loci, trait_name) {
  overlap_col <- paste0("overlap_", trait_name)
  overlap_n_col <- paste0("n_overlap_", trait_name)
  pleio_loci[, (overlap_col) := FALSE]
  pleio_loci[, (overlap_n_col) := 0L]
  if (nrow(trait_loci) == 0L) return(pleio_loci)

  for (i in seq_len(nrow(pleio_loci))) {
    hits <- trait_loci[
      chr_norm == pleio_loci$chr_norm[i] &
        start_int <= pleio_loci$end_int[i] &
        end_int >= pleio_loci$start_int[i]
    ]
    if (nrow(hits) > 0L) {
      pleio_loci[i, (overlap_col) := TRUE]
      pleio_loci[i, (overlap_n_col) := as.integer(nrow(hits))]
    }
  }
  pleio_loci
}

domain_sign <- function(z, threshold = 1.96) {
  z <- suppressWarnings(as.numeric(z))
  z <- z[!is.na(z) & abs(z) >= threshold]
  if (length(z) == 0L) return("unassigned")
  signs <- unique(sign(z))
  signs <- signs[signs != 0]
  if (length(signs) == 0L) return("unassigned")
  if (length(signs) > 1L) return("mixed")
  if (signs[1L] > 0) return("positive")
  "negative"
}

classify_snp <- function(ea_z, cf_z, mdd_z, scz_z, model, threshold = 1.96) {
  cognitive_domain <- domain_sign(c(ea_z, cf_z), threshold)
  psychiatric_z <- switch(
    model,
    "FOURTRAIT" = c(mdd_z, scz_z),
    "MDD_3TRAIT" = c(mdd_z),
    "SCZ_3TRAIT" = c(scz_z)
  )
  psychiatric_domain <- domain_sign(psychiatric_z, threshold)

  if (cognitive_domain == "mixed" || psychiatric_domain == "mixed") {
    return("mixed_within_domain")
  }
  if (cognitive_domain == "unassigned" || psychiatric_domain == "unassigned") {
    return("Unassigned")
  }
  if (!identical(cognitive_domain, psychiatric_domain)) return("concordant_raw")
  "discordant_raw"
}

collapse_locus <- function(x) {
  x <- unique(x[!is.na(x)])
  if (length(x) == 0L) return("Unassigned")
  if ("mixed_within_domain" %chin% x) return("Mixed")
  has_concordant <- "concordant_raw" %chin% x
  has_discordant <- "discordant_raw" %chin% x
  if (has_concordant && has_discordant) return("Dual")
  if (has_concordant) return("Concordant")
  if (has_discordant) return("Discordant")
  "Unassigned"
}

safe_sum <- function(x) sum(x, na.rm = TRUE)

# Read PLEIO and model-relevant single-trait FUMA loci.
pleio_loci <- read_fuma_loci(pleio_dir, "PLEIO", allow_zero = FALSE)
if (nrow(pleio_loci) == 0L) stop("PLEIO GenomicRiskLoci file contains zero rows.")
if (anyDuplicated(pleio_loci$locus_key)) stop("Duplicate PLEIO locus_key values detected.")

trait_dir_map <- list(
  EA = file.path(single_trait_root, ea_folder),
  CF = file.path(single_trait_root, cf_folder),
  MDD = file.path(single_trait_root, "MDD"),
  SCZ = file.path(single_trait_root, "SCZ")
)
trait_loci <- list()
for (trait in all_traits) {
  trait_loci[[trait]] <- read_fuma_loci(
    trait_dir_map[[trait]], trait, allow_zero = TRUE
  )
  pleio_loci <- overlap_with_trait(pleio_loci, trait_loci[[trait]], trait)
}

overlap_cols_model <- paste0("overlap_", model_traits)
pleio_loci[, overlaps_any_single_trait := Reduce(`|`, .SD), .SDcols = overlap_cols_model]
pleio_loci[, pleio_prioritised := !overlaps_any_single_trait]
pleio_loci[, prioritisation_status := ifelse(
  pleio_prioritised, "PLEIO_PRIORITISED", "OVERLAPS_SINGLE_TRAIT"
)]

# Read independent significant SNPs.
indsig_file <- find_unique_file(
  pleio_dir, "^IndSigSNPs.*\\.txt$", "PLEIO IndSigSNPs", allow_zero = FALSE
)
message("Reading PLEIO independent significant SNPs: ", indsig_file)
indsig <- fread(indsig_file)
if (nrow(indsig) == 0L) stop("PLEIO IndSigSNPs file contains zero rows.")

snp_col <- find_col(indsig, c("SNP", "rsID", "rsid", "IndSigSNP", "uniqID"))
chr_col <- find_col(indsig, c("chr", "CHR", "chromosome"), required = FALSE)
pos_col <- find_col(indsig, c("pos", "BP", "bp", "position"), required = FALSE)
locus_col <- find_col(
  indsig, c("GenomicLocus", "genomicLocus", "locus", "Locus"), required = FALSE
)
setnames(indsig, snp_col, "SNP")
indsig[, SNP := as.character(SNP)]
if (anyNA(indsig$SNP) || any(indsig$SNP == "")) stop("Missing SNP identifiers in IndSigSNPs.")
if (anyDuplicated(indsig$SNP)) {
  duplicated_snps <- unique(indsig$SNP[duplicated(indsig$SNP)])
  stop("Duplicate SNP identifiers in IndSigSNPs: ", paste(head(duplicated_snps, 20L), collapse = ", "))
}
if (!is.na(chr_col)) indsig[, chr_norm_fuma := normalise_chr(get(chr_col))]
if (!is.na(pos_col)) indsig[, bp_fuma := suppressWarnings(as.integer(get(pos_col)))]
if (!is.na(locus_col)) indsig[, fuma_locus_id := as.character(get(locus_col))]

# Read merged PLEIO statistics and normalise source aliases to EA and CF.
message("Reading raw merged PLEIO file...")
raw <- fread(raw_file)
raw_snp_col <- find_col(raw, c("SNP", "rsID", "rsid", "uniqID"))
ea_z_col <- find_col(raw, c("EA_Z", "EDU_Z", "Edu_Z"))
cf_z_col <- find_col(
  raw,
  c("CF_Z", "HELIOS10K_Z", "COG_Z", "G_COGENT_Z", "Cognition_Z", "G_Z")
)
mdd_z_col <- find_col(raw, c("MDD_Z"), required = model != "SCZ_3TRAIT")
scz_z_col <- find_col(raw, c("SCZ_Z"), required = model != "MDD_3TRAIT")

optional_cols <- c("CHR", "BP", "pleio_p", "LS_p", "min_single_trait_p")
keep_cols <- unique(c(raw_snp_col, optional_cols, ea_z_col, cf_z_col, mdd_z_col, scz_z_col))
keep_cols <- keep_cols[!is.na(keep_cols) & keep_cols %in% names(raw)]
raw_sub <- raw[, ..keep_cols]
setnames(raw_sub, raw_snp_col, "SNP")
setnames(raw_sub, ea_z_col, "EA_Z")
setnames(raw_sub, cf_z_col, "CF_Z")
if (!is.na(mdd_z_col)) setnames(raw_sub, mdd_z_col, "MDD_Z")
if (!is.na(scz_z_col)) setnames(raw_sub, scz_z_col, "SCZ_Z")
if (!("MDD_Z" %in% names(raw_sub))) raw_sub[, MDD_Z := NA_real_]
if (!("SCZ_Z" %in% names(raw_sub))) raw_sub[, SCZ_Z := NA_real_]
raw_sub[, SNP := as.character(SNP)]
if ("CHR" %in% names(raw_sub)) raw_sub[, CHR := normalise_chr(CHR)]
if ("BP" %in% names(raw_sub)) raw_sub[, BP := suppressWarnings(as.integer(BP))]

if (anyDuplicated(raw_sub$SNP)) {
  duplicated_snps <- unique(raw_sub$SNP[duplicated(raw_sub$SNP)])
  stop(
    "Duplicate SNP identifiers in raw merged file: ",
    paste(head(duplicated_snps, 20L), collapse = ", "),
    ". Resolve duplicates before classification."
  )
}

# Classify SNPs and assign them to PLEIO loci.
indsig_z <- merge(indsig, raw_sub, by = "SNP", all.x = TRUE, sort = FALSE)
indsig_z[, snp_direction_class := mapply(
  classify_snp,
  EA_Z, CF_Z, MDD_Z, SCZ_Z,
  MoreArgs = list(model = model, threshold = zthr)
)]

if ("fuma_locus_id" %in% names(indsig_z)) {
  indsig_z[, locus_key := paste0("PLEIO_", fuma_locus_id)]
} else {
  indsig_z[, locus_key := NA_character_]
  for (i in seq_len(nrow(indsig_z))) {
    chr_i <- if ("CHR" %in% names(indsig_z) && !is.na(indsig_z$CHR[i])) {
      indsig_z$CHR[i]
    } else if ("chr_norm_fuma" %in% names(indsig_z)) {
      indsig_z$chr_norm_fuma[i]
    } else {
      NA_character_
    }
    bp_i <- if ("BP" %in% names(indsig_z) && !is.na(indsig_z$BP[i])) {
      indsig_z$BP[i]
    } else if ("bp_fuma" %in% names(indsig_z)) {
      indsig_z$bp_fuma[i]
    } else {
      NA_integer_
    }
    if (is.na(chr_i) || is.na(bp_i)) next

    hits <- pleio_loci[chr_norm == chr_i & start_int <= bp_i & end_int >= bp_i]
    if (nrow(hits) == 1L) {
      indsig_z$locus_key[i] <- hits$locus_key[1L]
    } else if (nrow(hits) > 1L) {
      stop("SNP ", indsig_z$SNP[i], " maps to more than one PLEIO locus by coordinates.")
    }
  }
}

unknown_locus_keys <- setdiff(unique(na.omit(indsig_z$locus_key)), pleio_loci$locus_key)
if (length(unknown_locus_keys) > 0L) {
  stop(
    "Independent significant SNPs reference locus IDs absent from GenomicRiskLoci: ",
    paste(head(unknown_locus_keys, 20L), collapse = ", ")
  )
}

locus_classes <- indsig_z[
  !is.na(locus_key),
  .(
    locus_direction_class = collapse_locus(snp_direction_class),
    n_ind_sig_snps_classified = .N,
    n_concordant_raw_snps = safe_sum(snp_direction_class == "concordant_raw"),
    n_discordant_raw_snps = safe_sum(snp_direction_class == "discordant_raw"),
    n_mixed_within_domain_snps = safe_sum(snp_direction_class == "mixed_within_domain"),
    n_unassigned_snps = safe_sum(snp_direction_class == "Unassigned")
  ),
  by = locus_key
]

pleio_loci <- merge(pleio_loci, locus_classes, by = "locus_key", all.x = TRUE, sort = FALSE)
pleio_loci[is.na(locus_direction_class), locus_direction_class := "Unassigned"]
count_cols <- c(
  "n_ind_sig_snps_classified", "n_concordant_raw_snps",
  "n_discordant_raw_snps", "n_mixed_within_domain_snps", "n_unassigned_snps"
)
for (field in count_cols) pleio_loci[is.na(get(field)), (field) := 0L]

pleio_loci[, `:=`(ancestry = ancestry, PLEIO_model = model)]
indsig_z[, `:=`(ancestry = ancestry, PLEIO_model = model)]
setcolorder(pleio_loci, c("ancestry", "PLEIO_model", setdiff(names(pleio_loci), c("ancestry", "PLEIO_model"))))
setcolorder(indsig_z, c("ancestry", "PLEIO_model", setdiff(names(indsig_z), c("ancestry", "PLEIO_model"))))

# Publication-facing summaries and QC.
locus_summary <- pleio_loci[, .(n_loci = .N), by = .(
  ancestry, PLEIO_model, locus_direction_class
)]
locus_summary[, percent_loci := n_loci / sum(n_loci), by = .(ancestry, PLEIO_model)]
class_order <- c("Concordant", "Discordant", "Dual", "Mixed", "Unassigned")
locus_summary[, class_order := match(locus_direction_class, class_order)]
setorder(locus_summary, ancestry, PLEIO_model, class_order)
locus_summary[, class_order := NULL]

prioritisation_summary <- data.table(
  ancestry = ancestry,
  PLEIO_model = model,
  n_original_pleio_loci = nrow(pleio_loci),
  n_overlapping_single_trait = safe_sum(pleio_loci$overlaps_any_single_trait),
  n_pleio_prioritised = safe_sum(pleio_loci$pleio_prioritised),
  pct_pleio_prioritised = safe_sum(pleio_loci$pleio_prioritised) / nrow(pleio_loci),
  n_loci_overlapping_EA = safe_sum(pleio_loci$overlap_EA),
  n_loci_overlapping_CF = safe_sum(pleio_loci$overlap_CF),
  n_loci_overlapping_MDD = safe_sum(pleio_loci$overlap_MDD),
  n_loci_overlapping_SCZ = safe_sum(pleio_loci$overlap_SCZ)
)

required_model_z <- switch(
  model,
  "FOURTRAIT" = c("EA_Z", "CF_Z", "MDD_Z", "SCZ_Z"),
  "MDD_3TRAIT" = c("EA_Z", "CF_Z", "MDD_Z"),
  "SCZ_3TRAIT" = c("EA_Z", "CF_Z", "SCZ_Z")
)
qc <- data.table(
  ancestry = ancestry,
  PLEIO_model = model,
  z_threshold = zthr,
  model_comparator_traits = paste(model_traits, collapse = ";"),
  n_pleio_loci = nrow(pleio_loci),
  n_ind_sig_snps = nrow(indsig_z),
  n_ind_sig_snps_missing_raw_match = safe_sum(
    rowSums(!is.na(indsig_z[, ..required_model_z])) == 0L
  ),
  n_ind_sig_snps_unassigned_locus = safe_sum(is.na(indsig_z$locus_key)),
  n_loci_unassigned = safe_sum(pleio_loci$locus_direction_class == "Unassigned"),
  n_pleio_prioritised = safe_sum(pleio_loci$pleio_prioritised)
)

path_out <- function(suffix) file.path(out_dir, paste0(prefix, suffix))
fwrite(pleio_loci, path_out("_PLEIO_locus_master.tsv"), sep = "\t", na = "NA")
fwrite(indsig_z, path_out("_PLEIO_IndSigSNPs_with_direction.tsv"), sep = "\t", na = "NA")
fwrite(locus_summary, path_out("_PLEIO_locus_direction_summary.tsv"), sep = "\t", na = "NA")
fwrite(
  prioritisation_summary,
  path_out("_PLEIO_locus_prioritisation_summary.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(qc, path_out("_classification_QC.tsv"), sep = "\t", na = "NA")

cat("\n================ LOCUS DIRECTION SUMMARY ================\n")
print(locus_summary)
cat("\n================ LOCUS PRIORITISATION SUMMARY ================\n")
print(prioritisation_summary)
cat("\n================ QC ================\n")
print(qc)

message("Done. Outputs written under: ", out_dir)
message("Locus-level output: ", path_out("_PLEIO_locus_master.tsv"))
message("SNP-level output: ", path_out("_PLEIO_IndSigSNPs_with_direction.tsv"))
message("Direction summary: ", path_out("_PLEIO_locus_direction_summary.tsv"))
message("QC: ", path_out("_classification_QC.tsv"))
