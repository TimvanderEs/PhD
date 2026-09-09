#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# ============================================================
# PLEIO/FUMA locus classification audit and corrected rerun
#
# This script performs two classifications on the same data:
#   1) legacy: same cognitive/psychiatric sign = concordant
#              and all four single-trait FUMA sets are comparators
#   2) corrected: opposite cognitive/psychiatric sign = concordant
#                 and only model-relevant traits are comparators
#
# It writes corrected outputs plus legacy-vs-corrected transition/QC files.
# ============================================================

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
  cat("\nUsage:\n",
      "Rscript classify_PLEIO_FUMA_loci_v2.R \\\n",
      "  --ancestry EUR \\\n",
      "  --model FOURTRAIT \\\n",
      "  --pleio-dir /path/to/model/FUMA/PLEIO \\\n",
      "  --single-trait-root /path/to/ancestry/FUMA \\\n",
      "  --raw-file /path/to/raw_merged_PLEIO.tsv[.gz] \\\n",
      "  --out-dir /path/to/new_output_dir \\\n",
      "  [--cf-folder G] [--zthr 1.96] [--prefix EUR_FOURTRAIT]\n\n",
      "Accepted models: FOURTRAIT, MDD_3TRAIT, SCZ_3TRAIT\n",
      "single-trait-root must contain EDU, <cf-folder>, MDD, and SCZ subdirectories.\n",
      sep = "")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  usage()
  quit(save = "no", status = 0L)
}

required_args <- c("ancestry", "model", "pleio-dir", "single-trait-root", "raw-file", "out-dir")
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
  stop("Unknown --model: ", args$model,
       ". Accepted: FOURTRAIT, MDD_3TRAIT, SCZ_3TRAIT")
)

pleio_dir <- normalizePath(args[["pleio-dir"]], mustWork = TRUE)
single_trait_root <- normalizePath(args[["single-trait-root"]], mustWork = TRUE)
raw_file <- normalizePath(args[["raw-file"]], mustWork = TRUE)
out_dir <- args[["out-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

cf_folder <- if ("cf-folder" %in% names(args)) args[["cf-folder"]] else "G"
zthr <- if ("zthr" %in% names(args)) as.numeric(args$zthr) else 1.96
if (!is.finite(zthr) || zthr < 0) stop("--zthr must be a non-negative number.")
prefix <- if ("prefix" %in% names(args)) args$prefix else paste(ancestry, model, sep = "_")

all_traits <- c("EDU", "G", "MDD", "SCZ")
model_traits <- switch(
  model,
  "FOURTRAIT" = c("EDU", "G", "MDD", "SCZ"),
  "MDD_3TRAIT" = c("EDU", "G", "MDD"),
  "SCZ_3TRAIT" = c("EDU", "G", "SCZ")
)

message("Ancestry: ", ancestry)
message("Model: ", model)
message("PLEIO FUMA directory: ", pleio_dir)
message("Single-trait FUMA root: ", single_trait_root)
message("Raw merged PLEIO file: ", raw_file)
message("Output directory: ", out_dir)
message("Direction threshold |Z| >= ", zthr)
message("Corrected model comparators: ", paste(model_traits, collapse = ", "))

find_unique_file <- function(dir, pattern, label, allow_zero = FALSE) {
  if (!dir.exists(dir)) {
    if (allow_zero) {
      message("Directory absent for ", label, ": ", dir, " ; treating as 0 loci.")
      return(NA_character_)
    }
    stop("Directory not found for ", label, ": ", dir)
  }
  files <- list.files(dir, recursive = TRUE, full.names = TRUE)
  hits <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
  if (length(hits) == 0L) {
    if (allow_zero) {
      message("No file matched for ", label, " in ", dir, " ; treating as 0 loci.")
      return(NA_character_)
    }
    stop("No file matched for ", label, " in ", dir,
         " using pattern ", pattern)
  }
  if (length(hits) > 1L) {
    stop("Multiple files matched for ", label, ":\n", paste(hits, collapse = "\n"),
         "\nMove old files or pass a clean canonical directory.")
  }
  normalizePath(hits[1L], mustWork = TRUE)
}

find_col <- function(dt, candidates, required = TRUE) {
  nms <- names(dt)
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(nms[idx[1L]])
  if (required) {
    stop("Could not find any of these columns: ", paste(candidates, collapse = ", "),
         "\nAvailable columns: ", paste(nms, collapse = ", "))
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
  locus_col <- find_col(dt, c("GenomicLocus", "genomicLocus", "locus", "Locus"), required = FALSE)

  dt[, trait := trait_name]
  dt[, chr_norm := normalise_chr(get(chr_col))]
  dt[, start_int := suppressWarnings(as.integer(get(start_col)))]
  dt[, end_int := suppressWarnings(as.integer(get(end_col)))]

  if (anyNA(dt$chr_norm) || anyNA(dt$start_int) || anyNA(dt$end_int)) {
    stop("Missing/invalid coordinates in ", f)
  }
  if (any(dt$start_int > dt$end_int)) {
    stop("At least one locus has start > end in ", f)
  }

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

domain_sign <- function(z, zthr = 1.96) {
  z <- suppressWarnings(as.numeric(z))
  z <- z[!is.na(z) & abs(z) >= zthr]

  if (length(z) == 0L) return("unclear")
  s <- unique(sign(z))
  s <- s[s != 0]
  if (length(s) == 0L) return("unclear")
  if (length(s) > 1L) return("mixed")
  if (s[1L] > 0) return("positive")
  "negative"
}

classify_snp <- function(edu_z, cog_z, mdd_z, scz_z,
                         model, zthr = 1.96,
                         concordant_relationship = c("inverse", "same")) {
  concordant_relationship <- match.arg(concordant_relationship)

  cog_domain <- domain_sign(c(edu_z, cog_z), zthr)
  psych_z <- switch(
    model,
    "FOURTRAIT" = c(mdd_z, scz_z),
    "MDD_3TRAIT" = c(mdd_z),
    "SCZ_3TRAIT" = c(scz_z)
  )
  psych_domain <- domain_sign(psych_z, zthr)

  if (cog_domain == "mixed" || psych_domain == "mixed") {
    return("mixed_within_domain")
  }
  if (cog_domain == "unclear" || psych_domain == "unclear") {
    return("trait_dominant_or_unclear")
  }

  same_sign <- identical(cog_domain, psych_domain)
  is_concordant <- if (concordant_relationship == "inverse") !same_sign else same_sign
  if (is_concordant) return("concordant_raw")
  "discordant_raw"
}

collapse_locus <- function(x) {
  x <- unique(x[!is.na(x)])
  if (length(x) == 0L) return("unassigned")
  if ("mixed_within_domain" %chin% x) return("mixed")

  has_conc <- "concordant_raw" %chin% x
  has_disc <- "discordant_raw" %chin% x

  if (has_conc && has_disc) return("mixed_dual")
  if (has_conc) return("concordant")
  if (has_disc) return("discordant")
  "trait_dominant_or_unclear"
}

safe_sum <- function(x) sum(x, na.rm = TRUE)

# ------------------------------------------------------------
# 1. Read PLEIO and single-trait FUMA loci
# ------------------------------------------------------------
pleio_loci <- read_fuma_loci(pleio_dir, "PLEIO", allow_zero = FALSE)
if (nrow(pleio_loci) == 0L) stop("PLEIO GenomicRiskLoci file contains zero rows.")
if (anyDuplicated(pleio_loci$locus_key)) stop("Duplicate PLEIO locus_key values detected.")

trait_dir_map <- list(
  EDU = file.path(single_trait_root, "EDU"),
  G = file.path(single_trait_root, cf_folder),
  MDD = file.path(single_trait_root, "MDD"),
  SCZ = file.path(single_trait_root, "SCZ")
)

trait_loci <- list()
for (tr in all_traits) {
  trait_loci[[tr]] <- read_fuma_loci(trait_dir_map[[tr]], tr, allow_zero = TRUE)
  pleio_loci <- overlap_with_trait(pleio_loci, trait_loci[[tr]], tr)
}

overlap_cols_all <- paste0("overlap_", all_traits)
overlap_cols_model <- paste0("overlap_", model_traits)

# Legacy prioritisation reproduced from the supplied script.
pleio_loci[, legacy_overlaps_any_single_trait := Reduce(`|`, .SD), .SDcols = overlap_cols_all]
pleio_loci[, legacy_novel_vs_single_trait_FUMA := !legacy_overlaps_any_single_trait]

# Corrected, model-specific prioritisation.
pleio_loci[, overlaps_any_relevant_single_trait := Reduce(`|`, .SD), .SDcols = overlap_cols_model]
pleio_loci[, pleio_prioritised_vs_relevant_single_trait_FUMA := !overlaps_any_relevant_single_trait]

# ------------------------------------------------------------
# 2. Read FUMA independent significant SNPs
# ------------------------------------------------------------
indsig_file <- find_unique_file(
  pleio_dir,
  "^IndSigSNPs.*\\.txt$",
  "PLEIO IndSigSNPs",
  allow_zero = FALSE
)
message("Reading PLEIO independent significant SNPs: ", indsig_file)
indsig <- fread(indsig_file)
if (nrow(indsig) == 0L) stop("PLEIO IndSigSNPs file contains zero rows.")

snp_col <- find_col(indsig, c("SNP", "rsID", "rsid", "IndSigSNP", "uniqID"))
chr_col <- find_col(indsig, c("chr", "CHR", "chromosome"), required = FALSE)
pos_col <- find_col(indsig, c("pos", "BP", "bp", "position"), required = FALSE)
locus_col <- find_col(indsig, c("GenomicLocus", "genomicLocus", "locus", "Locus"), required = FALSE)

setnames(indsig, snp_col, "SNP")
indsig[, SNP := as.character(SNP)]
if (anyNA(indsig$SNP) || any(indsig$SNP == "")) stop("Missing SNP identifiers in IndSigSNPs.")
if (anyDuplicated(indsig$SNP)) {
  dup <- unique(indsig$SNP[duplicated(indsig$SNP)])
  stop("Duplicate SNP identifiers in IndSigSNPs: ", paste(head(dup, 20L), collapse = ", "))
}

if (!is.na(chr_col)) indsig[, chr_norm_fuma := normalise_chr(get(chr_col))]
if (!is.na(pos_col)) indsig[, bp_fuma := suppressWarnings(as.integer(get(pos_col)))]
if (!is.na(locus_col)) indsig[, fuma_locus_id := as.character(get(locus_col))]

# ------------------------------------------------------------
# 3. Read raw merged PLEIO file and harmonise Z columns
# ------------------------------------------------------------
message("Reading raw merged PLEIO file...")
raw <- fread(raw_file)
raw_snp_col <- find_col(raw, c("SNP", "rsID", "rsid", "uniqID"))
edu_z_col <- find_col(raw, c("EDU_Z", "Edu_Z", "EA_Z"))
cog_z_col <- find_col(raw, c("G_Z", "HELIOS10K_Z", "COG_Z", "G_COGENT_Z", "Cognition_Z", "CF_Z"))
mdd_z_col <- find_col(raw, c("MDD_Z"), required = model != "SCZ_3TRAIT")
scz_z_col <- find_col(raw, c("SCZ_Z"), required = model != "MDD_3TRAIT")

# For a faithful legacy comparison in a 3-trait run, both psychiatric columns are needed.
if (is.na(mdd_z_col)) mdd_z_col <- find_col(raw, c("MDD_Z"), required = FALSE)
if (is.na(scz_z_col)) scz_z_col <- find_col(raw, c("SCZ_Z"), required = FALSE)
legacy_direction_available <- !is.na(mdd_z_col) && !is.na(scz_z_col)

optional_cols <- c("CHR", "BP", "pleio_p", "LS_p", "min_single_trait_p")
keep_cols <- unique(c(raw_snp_col, optional_cols, edu_z_col, cog_z_col, mdd_z_col, scz_z_col))
keep_cols <- keep_cols[!is.na(keep_cols) & keep_cols %in% names(raw)]
raw_sub <- raw[, ..keep_cols]
setnames(raw_sub, raw_snp_col, "SNP")
setnames(raw_sub, edu_z_col, "EDU_Z")
setnames(raw_sub, cog_z_col, "COG_Z")
if (!is.na(mdd_z_col)) setnames(raw_sub, mdd_z_col, "MDD_Z")
if (!is.na(scz_z_col)) setnames(raw_sub, scz_z_col, "SCZ_Z")
if (!("MDD_Z" %in% names(raw_sub))) raw_sub[, MDD_Z := NA_real_]
if (!("SCZ_Z" %in% names(raw_sub))) raw_sub[, SCZ_Z := NA_real_]

raw_sub[, SNP := as.character(SNP)]
if ("CHR" %in% names(raw_sub)) raw_sub[, CHR := normalise_chr(CHR)]
if ("BP" %in% names(raw_sub)) raw_sub[, BP := suppressWarnings(as.integer(BP))]

if (anyDuplicated(raw_sub$SNP)) {
  dup <- unique(raw_sub$SNP[duplicated(raw_sub$SNP)])
  stop("Duplicate SNP identifiers in raw merged file: ",
       paste(head(dup, 20L), collapse = ", "),
       ". Resolve duplicates before classification to avoid a many-to-many merge.")
}

# ------------------------------------------------------------
# 4. Merge, classify SNPs, and assign loci
# ------------------------------------------------------------
indsig_z <- merge(indsig, raw_sub, by = "SNP", all.x = TRUE, sort = FALSE)

indsig_z[, corrected_snp_direction_class := mapply(
  classify_snp,
  EDU_Z, COG_Z, MDD_Z, SCZ_Z,
  MoreArgs = list(
    model = model,
    zthr = zthr,
    concordant_relationship = "inverse"
  )
)]

if (legacy_direction_available) {
  # Reproduce the supplied script: always use all four traits and same-sign = concordant.
  indsig_z[, legacy_snp_direction_class := mapply(
    classify_snp,
    EDU_Z, COG_Z, MDD_Z, SCZ_Z,
    MoreArgs = list(
      model = "FOURTRAIT",
      zthr = zthr,
      concordant_relationship = "same"
    )
  )]
} else {
  warning("Both MDD_Z and SCZ_Z are not available; legacy direction comparison will be NA.")
  indsig_z[, legacy_snp_direction_class := NA_character_]
}

# Prefer the explicit FUMA locus ID; fall back to coordinate overlap.
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

    hits <- pleio_loci[
      chr_norm == chr_i & start_int <= bp_i & end_int >= bp_i
    ]
    if (nrow(hits) == 1L) {
      indsig_z$locus_key[i] <- hits$locus_key[1L]
    } else if (nrow(hits) > 1L) {
      stop("SNP ", indsig_z$SNP[i], " maps to more than one PLEIO locus by coordinates.")
    }
  }
}

unknown_locus_keys <- setdiff(unique(na.omit(indsig_z$locus_key)), pleio_loci$locus_key)
if (length(unknown_locus_keys) > 0L) {
  stop("Independent significant SNPs reference locus IDs absent from GenomicRiskLoci: ",
       paste(head(unknown_locus_keys, 20L), collapse = ", "))
}

corrected_locus_class <- indsig_z[
  !is.na(locus_key),
  .(
    corrected_locus_direction_class = collapse_locus(corrected_snp_direction_class),
    n_ind_sig_snps_classified = .N,
    n_corrected_concordant_raw_snps = safe_sum(corrected_snp_direction_class == "concordant_raw"),
    n_corrected_discordant_raw_snps = safe_sum(corrected_snp_direction_class == "discordant_raw"),
    n_corrected_mixed_within_domain_snps = safe_sum(corrected_snp_direction_class == "mixed_within_domain"),
    n_corrected_trait_dominant_or_unclear_snps = safe_sum(corrected_snp_direction_class == "trait_dominant_or_unclear")
  ),
  by = locus_key
]

if (legacy_direction_available) {
  legacy_locus_class <- indsig_z[
    !is.na(locus_key),
    .(
      legacy_locus_direction_class = collapse_locus(legacy_snp_direction_class),
      n_legacy_concordant_raw_snps = safe_sum(legacy_snp_direction_class == "concordant_raw"),
      n_legacy_discordant_raw_snps = safe_sum(legacy_snp_direction_class == "discordant_raw"),
      n_legacy_mixed_within_domain_snps = safe_sum(legacy_snp_direction_class == "mixed_within_domain"),
      n_legacy_trait_dominant_or_unclear_snps = safe_sum(legacy_snp_direction_class == "trait_dominant_or_unclear")
    ),
    by = locus_key
  ]
} else {
  legacy_locus_class <- unique(indsig_z[!is.na(locus_key), .(locus_key)])
  legacy_locus_class[, `:=`(
    legacy_locus_direction_class = NA_character_,
    n_legacy_concordant_raw_snps = NA_integer_,
    n_legacy_discordant_raw_snps = NA_integer_,
    n_legacy_mixed_within_domain_snps = NA_integer_,
    n_legacy_trait_dominant_or_unclear_snps = NA_integer_
  )]
}

pleio_loci <- merge(pleio_loci, corrected_locus_class, by = "locus_key", all.x = TRUE, sort = FALSE)
pleio_loci <- merge(pleio_loci, legacy_locus_class, by = "locus_key", all.x = TRUE, sort = FALSE)
pleio_loci[is.na(corrected_locus_direction_class), corrected_locus_direction_class := "unassigned"]
if (legacy_direction_available) {
  pleio_loci[is.na(legacy_locus_direction_class), legacy_locus_direction_class := "unassigned"]
}

# Canonical corrected aliases for downstream use.
pleio_loci[, locus_direction_class := corrected_locus_direction_class]
pleio_loci[, novel_vs_single_trait_FUMA := pleio_prioritised_vs_relevant_single_trait_FUMA]
indsig_z[, snp_direction_class := corrected_snp_direction_class]

# ------------------------------------------------------------
# 5. Summaries, transitions, and QC
# ------------------------------------------------------------
# data.table does not recycle length-1 external values inside by=.().
# Group using row-level columns first, then append run metadata.
ancestry_value <- ancestry
model_value <- model

corrected_summary <- pleio_loci[
  , .N,
  by = .(
    locus_direction_class = corrected_locus_direction_class,
    pleio_prioritised =
      pleio_prioritised_vs_relevant_single_trait_FUMA
  )
][order(locus_direction_class, pleio_prioritised)]

corrected_summary[, `:=`(
  ancestry = ancestry_value,
  model = model_value
)]

setcolorder(
  corrected_summary,
  c(
    "ancestry",
    "model",
    "locus_direction_class",
    "pleio_prioritised",
    "N"
  )
)

legacy_summary <- pleio_loci[
  , .N,
  by = .(
    locus_direction_class = legacy_locus_direction_class,
    legacy_novel = legacy_novel_vs_single_trait_FUMA
  )
][order(locus_direction_class, legacy_novel)]

legacy_summary[, `:=`(
  ancestry = ancestry_value,
  model = model_value
)]

setcolorder(
  legacy_summary,
  c(
    "ancestry",
    "model",
    "locus_direction_class",
    "legacy_novel",
    "N"
  )
)

locus_transition <- pleio_loci[
  , .N,
  by = .(
    legacy_locus_direction_class,
    corrected_locus_direction_class
  )
][order(legacy_locus_direction_class, corrected_locus_direction_class)]

prioritisation_transition <- pleio_loci[
  , .N,
  by = .(
    legacy_novel_vs_single_trait_FUMA,
    pleio_prioritised_vs_relevant_single_trait_FUMA
  )
][order(legacy_novel_vs_single_trait_FUMA,
        pleio_prioritised_vs_relevant_single_trait_FUMA)]

snp_transition <- indsig_z[
  , .N,
  by = .(
    legacy_snp_direction_class,
    corrected_snp_direction_class
  )
][order(legacy_snp_direction_class, corrected_snp_direction_class)]

required_model_z <- switch(
  model,
  "FOURTRAIT" = c("EDU_Z", "COG_Z", "MDD_Z", "SCZ_Z"),
  "MDD_3TRAIT" = c("EDU_Z", "COG_Z", "MDD_Z"),
  "SCZ_3TRAIT" = c("EDU_Z", "COG_Z", "SCZ_Z")
)

qc <- data.table(
  ancestry = ancestry,
  model = model,
  z_threshold = zthr,
  n_pleio_loci = nrow(pleio_loci),
  n_ind_sig_snps = nrow(indsig_z),
  n_ind_sig_snps_missing_raw_match = safe_sum(
    rowSums(!is.na(indsig_z[, ..required_model_z])) == 0L
  ),
  n_ind_sig_snps_unassigned_locus = safe_sum(is.na(indsig_z$locus_key)),
  n_loci_unassigned_corrected = safe_sum(pleio_loci$corrected_locus_direction_class == "unassigned"),
  n_loci_direction_changed = if (legacy_direction_available) {
    safe_sum(pleio_loci$legacy_locus_direction_class != pleio_loci$corrected_locus_direction_class)
  } else NA_integer_,
  n_loci_prioritisation_changed = safe_sum(
    pleio_loci$legacy_novel_vs_single_trait_FUMA !=
      pleio_loci$pleio_prioritised_vs_relevant_single_trait_FUMA
  ),
  n_legacy_prioritised = safe_sum(pleio_loci$legacy_novel_vs_single_trait_FUMA),
  n_corrected_prioritised = safe_sum(pleio_loci$pleio_prioritised_vs_relevant_single_trait_FUMA),
  legacy_direction_comparison_available = legacy_direction_available
)

# ------------------------------------------------------------
# 6. Write outputs
# ------------------------------------------------------------
path_out <- function(suffix) file.path(out_dir, paste0(prefix, suffix))

# Corrected canonical outputs, with legacy/audit columns retained in the master.
fwrite(pleio_loci, path_out("_PLEIO_locus_master.tsv"), sep = "\t", na = "NA")
fwrite(indsig_z, path_out("_PLEIO_IndSigSNPs_with_direction.tsv"), sep = "\t", na = "NA")
fwrite(corrected_summary, path_out("_PLEIO_locus_direction_summary.tsv"), sep = "\t", na = "NA")

# Explicit audit/comparison outputs.
fwrite(legacy_summary, path_out("_legacy_locus_direction_summary.tsv"), sep = "\t", na = "NA")
fwrite(locus_transition, path_out("_locus_direction_transition.tsv"), sep = "\t", na = "NA")
fwrite(snp_transition, path_out("_snp_direction_transition.tsv"), sep = "\t", na = "NA")
fwrite(prioritisation_transition, path_out("_prioritisation_transition.tsv"), sep = "\t", na = "NA")
fwrite(qc, path_out("_classification_QC.tsv"), sep = "\t", na = "NA")

cat("\n================ CORRECTED LOCUS SUMMARY ================\n")
print(corrected_summary)
cat("\n================ LEGACY -> CORRECTED LOCUS TRANSITIONS ================\n")
print(locus_transition)
cat("\n================ PRIORITISATION TRANSITIONS ================\n")
print(prioritisation_transition)
cat("\n================ QC ================\n")
print(qc)

message("Done. Outputs written under: ", out_dir)
message("Primary corrected master: ", path_out("_PLEIO_locus_master.tsv"))
message("Primary corrected summary: ", path_out("_PLEIO_locus_direction_summary.tsv"))
message("Audit transition table: ", path_out("_locus_direction_transition.tsv"))
message("QC: ", path_out("_classification_QC.tsv"))
