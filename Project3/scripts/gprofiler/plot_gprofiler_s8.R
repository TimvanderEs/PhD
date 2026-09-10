#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || any(args %in% c("-h", "--help"))) {
  cat(
    paste0(
      "Usage: Rscript plot_gprofiler_s8.R SIGNIFICANT_TERMS.tsv OUTPUT_DIR ",
      "[TOP_PER_QUERY=8] [PLOT_ALL=false]\n"
    )
  )
  quit(status = if (any(args %in% c("-h", "--help"))) 0L else 2L)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

input_file <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
top_per_query <- if (length(args) >= 3L) as.integer(args[[3]]) else 8L
plot_all <- if (length(args) >= 4L) tolower(args[[4]]) == "true" else FALSE
if (!is.finite(top_per_query) || top_per_query < 1L) {
  stop("TOP_PER_QUERY must be a positive integer", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gp <- fread(input_file, na.strings = c("", "NA"))
required <- c(
  "ancestry", "model", "direction_class", "source_label", "term_id",
  "term_name", "adjusted_p_value", "negative_log10_adjP",
  "intersection_size", "broad_category"
)
missing <- setdiff(required, names(gp))
if (length(missing)) {
  stop("Missing input columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

model_labels <- c(
  FOURTRAIT = "Four-trait",
  MDD_3TRAIT = "MDD three-trait",
  SCZ_3TRAIT = "SCZ three-trait",
  EAS_4TRAIT = "Four-trait",
  EAS_MDD_EA_CF = "MDD three-trait",
  EAS_SCZ_EA_CF = "SCZ three-trait"
)
direction_labels <- c(
  concordant = "Concordant",
  discordant = "Discordant",
  dual = "Dual",
  all_mapped_genes = "All mapped genes"
)

wrap_text <- function(x, width = 46L) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

gp <- gp[is.finite(adjusted_p_value) & adjusted_p_value < 0.05]
gp[, model_plot := unname(model_labels[model])]
gp[is.na(model_plot), model_plot := model]
gp[, direction_plot := unname(direction_labels[direction_class])]
gp[is.na(direction_plot), direction_plot := direction_class]
gp[, query := ifelse(
  ancestry == "EAS",
  model_plot,
  paste0(model_plot, " · ", direction_plot)
)]
gp[, term_wrapped := wrap_text(term_name)]

select_terms <- function(x) {
  if (plot_all) return(copy(x))
  x[order(adjusted_p_value), head(.SD, top_per_query), by = query]
}

eur <- select_terms(gp[ancestry == "EUR"])
eas <- select_terms(gp[ancestry == "EAS"])

make_ora_plot <- function(x, title, subtitle) {
  if (!nrow(x)) {
    return(ggplot() + theme_void() + labs(title = title, subtitle = "No significant terms"))
  }
  query_order <- unique(x[order(model_plot, direction_plot), query])
  x[, query := factor(query, levels = query_order)]
  setorder(x, query, adjusted_p_value)
  x[, term_key := paste0(term_wrapped, "___", query)]
  x[, term_key := factor(term_key, levels = rev(unique(term_key)))]

  ggplot(x, aes(negative_log10_adjP, term_key)) +
    geom_segment(
      aes(x = 0, xend = negative_log10_adjP, yend = term_key),
      colour = "grey78", linewidth = 0.45
    ) +
    geom_point(aes(size = intersection_size, colour = broad_category), alpha = 0.92) +
    facet_grid(
      rows = vars(query), scales = "free_y", space = "free_y", switch = "y"
    ) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_colour_viridis_d(option = "D", end = 0.9) +
    scale_size_continuous(range = c(2, 7.5)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(-log[10]("g:SCS-adjusted " * italic(P))),
      y = NULL,
      colour = "Term family",
      size = "Intersection size"
    ) +
    theme_bw(base_size = 8.5) +
    theme(
      legend.position = "right",
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      panel.grid.major.y = element_blank()
    )
}

p_eur <- make_ora_plot(
  eur,
  "EUR directional over-representation results",
  "Class-specific mapped genes; model-wide mapped-gene custom backgrounds"
)
p_eas <- make_ora_plot(
  eas,
  "EAS model-wide over-representation results",
  "Model-wide mapped-gene queries; original background option was not recorded"
)
figure_s8 <- p_eur / p_eas + plot_annotation(tag_levels = "A")
height <- max(15, 5 + 0.23 * (nrow(eur) + nrow(eas)))
stem <- file.path(output_dir, "Supplementary_Figure_S8_gProfiler_overrepresentation")
ggsave(paste0(stem, ".pdf"), figure_s8, width = 14, height = height, units = "in")
ggsave(paste0(stem, ".png"), figure_s8, width = 14, height = height, units = "in", dpi = 600)
ggsave(
  paste0(stem, ".tiff"), figure_s8, width = 14, height = height,
  units = "in", dpi = 600, compression = "lzw"
)

setorder(gp, ancestry, model_plot, direction_plot, adjusted_p_value)
setorder(eur, query, adjusted_p_value)
setorder(eas, query, adjusted_p_value)
fwrite(gp, file.path(output_dir, "S8_all_significant_gProfiler_terms.csv"))
fwrite(rbindlist(list(eur, eas), use.names = TRUE), file.path(output_dir, "S8_terms_plotted.csv"))

message("Wrote S8 files to: ", normalizePath(output_dir))
