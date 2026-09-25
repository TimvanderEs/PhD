#!/usr/bin/env Rscript

library(LAVA)

## ------------ PATHS (EDIT IF NEEDED) ----------------

BASE <- "/home/ec2-user/software/COLOC-reporter"

input.info.file     <- file.path(BASE, "input.info_EAS.txt")
sample.overlap.file <- file.path(BASE, "sample.overlap.txt")
locus.file          <- file.path(BASE, "EAS_LAVA.locfile")

#input.info.file     <- file.path(BASE, "input.info_EAS_noHEL10k.txt")i
#sample.overlap.file <- file.path(BASE, "sample_overlap_EAS_noHEL10k.txt")
#locus.file          <- file.path(BASE, "EAS_LAVA_40loci.locfile")

# Trait IDs (must match ID column + sample.overlap row/colnames)
#phenos <- c("MDD_EAS", "Edu_EAS", "SCZ_EAS")


# reference prefix (no .bed/.bim/.fam extension)
ref.prefix          <- file.path("/home/ec2-user/input_files/reference/g1000_eas")

output.path         <- file.path(BASE, "results_EAS")
dir.create(output.path, showWarnings = FALSE, recursive = TRUE)

# Trait IDs (must match ID column + sample.overlap row/colnames)
phenos <- c("MDD_EAS", "Edu_EAS", "SCZ_EAS", "HEL10k_EAS")

univ.p.thresh <- 2e-5

## ------------ LAVA PIPELINE -------------------------

cat("Reading loci from:", locus.file, "\n")
loci  <- read.loci(locus.file)
n.loc <- nrow(loci)
cat("Number of loci:", n.loc, "\n")

cat("Processing input...\n")
input <- process.input(
  input.info = input.info.file,
  sample.overlap = sample.overlap.file,
  ref.prefix = ref.prefix,
  phenos = phenos
)

cat("Starting LAVA analysis for", n.loc, "loci\n")

progress <- ceiling(quantile(1:n.loc, seq(.05, 1, .05)))
u <- vector("list", n.loc)
b <- vector("list", n.loc)

for (i in seq_len(n.loc)) {
  if (i %in% progress) {
    cat("..", names(progress)[which(progress == i)], " (locus", i, "of", n.loc, ")\n")
    flush.console()
  }

  locus <- process.locus(loci[i, ], input)

  if (!is.null(locus)) {
    loc.info <- data.frame(
      locus  = locus$id,
      chr    = locus$chr,
      start  = locus$start,
      stop   = locus$stop,
      n.snps = locus$n.snps,
      n.pcs  = locus$K
    )

    loc.out <- run.univ.bivar(locus, univ.thresh = univ.p.thresh)

    u[[i]] <- cbind(loc.info, loc.out$univ)
    if (!is.null(loc.out$bivar)) {
      b[[i]] <- cbind(loc.info, loc.out$bivar)
    }
  }
}

cat("Collapsing output and writing files...\n")

## Remove NULL entries (loci with no results)
u <- u[!sapply(u, is.null)]
b <- b[!sapply(b, is.null)]

## Harmonise columns across loci before rbind

if (length(u) > 0L) {
  all_u_cols <- Reduce(union, lapply(u, colnames))
  u <- lapply(u, function(df) {
    missing <- setdiff(all_u_cols, colnames(df))
    for (m in missing) df[[m]] <- NA
    df[all_u_cols]
  })
  u.out <- do.call(rbind, u)
} else {
  u.out <- NULL
}

if (length(b) > 0L) {
  all_b_cols <- Reduce(union, lapply(b, colnames))
  b <- lapply(b, function(df) {
    missing <- setdiff(all_b_cols, colnames(df))
    for (m in missing) df[[m]] <- NA
    df[all_b_cols]
  })
  b.out <- do.call(rbind, b)
} else {
  b.out <- NULL
}

## Output filenames (now no HEL10k in the name)

univ.out.file  <- file.path(output.path, "EAS_MDD_Edu_SCZ_HELIOS.univ.lava")
bivar.out.file <- file.path(output.path, "EAS_MDD_Edu_SCZ_HELIOS.bivar.lava")

if (!is.null(u.out)) {
  write.table(
    u.out,
    file      = univ.out.file,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
}

if (!is.null(b.out)) {
  write.table(
    b.out,
    file      = bivar.out.file,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
}

cat("Done!\n")
cat("Univariate output:", univ.out.file,  "\n")
cat("Bivariate output: ", bivar.out.file, "\n")

