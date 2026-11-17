#!/usr/bin/env Rscript
# Restyle frozen Cook-derived EMT scores; no expression/scoring pipeline is run.
# Layout adapted from FigureYa175quadrant (top/right distributions + quadrants):
# https://github.com/ying-ge/FigureYa/blob/main/FigureYa175quadrant/FigureYa175quadrant.Rmd
# Histogram marginals are replaced by geom_density, also used in FigureYa227boxdensity.
# Reference style: /Users/sodajay/Desktop/EMT_Resutls/code/02_figure1_arial.R
# FigureYa citation: Lu et al. (2025), iMetaMed, doi:10.1002/imm3.70005.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})
options(stringsAsFactors = FALSE)
cowplot::set_null_device("agg")
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
project <- normalizePath(file.path(dirname(script), "../.."))
out <- file.path(project, "EMT_Resutls/figures/Cook_EMT_BCD")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source_file <- file.path(project, "results_v4/tables/S5_Cook2020_two_axis_EMT_scores_by_cell.csv")
summary_file <- file.path(project, "results_v4/tables/S5_Cook2020_two_axis_EMT_quadrant_by_cohort.csv")
data <- fread(source_file)
reference <- fread(summary_file)
cohorts <- c("China-HCC", "USA-HCC", "CCA-combined")
keys <- c("ChinaHCC", "USAHCC", "CCAcombined")
states <- c("Epithelial-high", "E/M-high", "EMT-high", "Low-both")
# Reuse four colours from the completed UMAP; these encode score quadrants here,
# not cell-type identities. A separate legend makes that distinction explicit.
cols <- c("Epithelial-high" = "#7895B2", "E/M-high" = "#E89B86",
          "EMT-high" = "#9CCDBD", "Low-both" = "#A3ADCA")
ink <- "#222222"
font <- "Arial"
stopifnot(nrow(data) == 29648L, uniqueN(data$barcode) == 29648L,
          all(is.finite(data$epithelial_score)), all(is.finite(data$hallmark_EMT_score)),
          identical(vapply(cohorts, function(x) nrow(data[cohort == x]), integer(1)),
                    setNames(c(13430L, 10914L, 5304L), cohorts)))
xcut <- median(data$epithelial_score)
ycut <- median(data$hallmark_EMT_score)
xrange <- as.numeric(quantile(data$epithelial_score, c(.005, .995)))
yrange <- as.numeric(quantile(data$hallmark_EMT_score, c(.005, .995)))
check_quadrant <- with(data, ifelse(epithelial_score > xcut,
  ifelse(hallmark_EMT_score > ycut, "E/M-high", "Epithelial-high"),
  ifelse(hallmark_EMT_score > ycut, "EMT-high", "Low-both")))
stopifnot(identical(check_quadrant, data$quadrant))
summary <- data[, .(cells = .N), by = .(cohort, quadrant)]
summary[, percentage := cells / sum(cells) * 100, by = cohort]
comparison <- merge(summary, reference, by = c("cohort", "quadrant"), suffixes = c("", "_reference"))
stopifnot(nrow(comparison) == 12L,
          all(comparison$cells == comparison$cells_reference),
          max(abs(comparison$percentage - comparison$percentage_reference)) < 1e-10)
data[, quadrant := factor(quadrant, levels = states)]

arial_pdf <- function(filename, width, height, bg = "white", ...) {
  grDevices::quartz(type = "pdf", file = filename, width = width, height = height,
                   family = font, bg = bg)
}

make_panel <- function(cohort_name, letter, standalone = FALSE) {
  z <- data[cohort == cohort_name]
  zs <- summary[cohort == cohort_name]
  set.seed(20260827L + match(cohort_name, cohorts))
  points <- z[sample.int(.N, min(.N, 6000L))]
  # Same sampled cells and all-cell density inputs as the original Figure 3.
  dx <- .025 * diff(xrange); dy <- .035 * diff(yrange)
  labels <- data.table(quadrant = states,
    x = c(xrange[2] - dx, xrange[2] - dx, xrange[1] + dx, xrange[1] + dx),
    y = c(yrange[1] + dy, yrange[2] - dy, yrange[2] - dy, yrange[1] + dy),
    hj = c(1, 1, 0, 0), vj = c(0, 1, 1, 0))
  labels <- merge(labels, zs[, .(quadrant, percentage)], by = "quadrant", sort = FALSE)
  labels[, label := sprintf("%s\n%.1f%%", quadrant, percentage)]
  label_pt <- if (standalone) 10 else 8
  tick_pt <- if (standalone) 10 else 9
  axis_pt <- if (standalone) 11 else 10
  p <- ggplot(points, aes(epithelial_score, hallmark_EMT_score)) +
    stat_density_2d(data = z, bins = 6, colour = "#777777", alpha = .4,
                    linewidth = .22, show.legend = FALSE) +
    geom_point(aes(colour = quadrant), size = if (standalone) .37 else .24,
               alpha = .45, shape = 16, stroke = 0) +
    geom_vline(xintercept = xcut, colour = "#555555", linetype = "dashed", linewidth = .32) +
    geom_hline(yintercept = ycut, colour = "#555555", linetype = "dashed", linewidth = .32) +
    geom_label(data = labels, aes(x, y, label = label, hjust = hj, vjust = vj),
               inherit.aes = FALSE, family = font, size = label_pt / .pt,
               fontface = "bold", lineheight = 1.05, colour = ink,
               fill = scales::alpha("white", .88), linewidth = 0,
               label.padding = unit(.08, "lines"), label.r = unit(.04, "lines")) +
    scale_colour_manual(values = cols, breaks = states, drop = FALSE) +
    scale_x_continuous(breaks = c(0, .3, .6, .9)) +
    scale_y_continuous(breaks = c(0, .1, .2), labels = c("0.0", "0.1", "0.2")) +
    coord_cartesian(xlim = xrange, ylim = yrange, expand = FALSE) +
    labs(x = "Epithelial-associated score", y = "Hallmark EMT score") +
    theme_classic(base_family = font, base_size = tick_pt) +
    theme(legend.position = "none",
          panel.border = element_rect(fill = NA, colour = ink, linewidth = .35),
          axis.line = element_blank(),
          axis.ticks = element_line(colour = ink, linewidth = .3),
          axis.ticks.length = unit(2, "pt"),
          axis.text = element_text(colour = ink, size = tick_pt),
          axis.title = element_text(colour = ink, size = axis_pt),
          axis.title.x = element_text(margin = margin(t = 4)),
          axis.title.y = element_text(margin = margin(r = 4)),
          plot.margin = margin(1, 1, 1, 1))
  top <- ggplot(z, aes(epithelial_score)) +
    geom_density(fill = cols[["Epithelial-high"]], colour = cols[["Epithelial-high"]],
                 alpha = .65, linewidth = .4, adjust = .85) +
    geom_vline(xintercept = xcut, colour = "#555555", linetype = "dashed", linewidth = .32) +
    coord_cartesian(xlim = xrange, expand = FALSE) +
    theme_void(base_family = font) + theme(plot.margin = margin(0, 1, 1, 1))
  right <- ggplot(z, aes(hallmark_EMT_score)) +
    geom_density(fill = cols[["EMT-high"]], colour = cols[["EMT-high"]],
                 alpha = .8, linewidth = .4, adjust = .85) +
    geom_vline(xintercept = ycut, colour = "#555555", linetype = "dashed", linewidth = .32) +
    coord_flip(xlim = yrange, expand = FALSE) +
    theme_void(base_family = font) + theme(plot.margin = margin(1, 0, 1, 1))
  # FigureYa175's marginal/main ratio is 1:6. Patchwork aligns the data panels.
  panel <- wrap_plots(top, plot_spacer(), p, right, ncol = 2,
                     widths = c(6, 1), heights = c(1, 6))
  ggdraw() + draw_plot(panel, x = 0, y = 0, width = 1, height = .86) +
    draw_label(cohort_name, x = .54, y = .967, hjust = .5, vjust = 1,
               fontfamily = font, fontface = "bold", size = 14, colour = ink) +
    draw_label(sprintf("%s malignant epithelial cells", format(nrow(z), big.mark = ",")),
               x = .54, y = .887, hjust = .5, vjust = 0,
               fontfamily = font, size = if (standalone) 10 else 9, colour = "#555555") +
    draw_label(letter, x = .015, y = .98, hjust = 0, vjust = 1,
               fontfamily = font, fontface = "bold", size = 16, colour = ink)
}

legend_data <- data.table(state = factor(states, levels = states), x = seq_along(states), y = 1)
legend_plot <- ggplot(legend_data, aes(x, y, colour = state)) + geom_point() +
  scale_colour_manual(values = cols, name = NULL, breaks = states, drop = FALSE) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.6, alpha = 1))) +
  theme_void(base_family = font) +
  theme(legend.position = "bottom", legend.text = element_text(size = 10, colour = ink),
        legend.key.width = unit(10, "pt"), legend.spacing.x = unit(4, "pt"),
        legend.margin = margin(0, 0, 0, 0))
legend <- get_legend(legend_plot)

save_pair <- function(p, stem, width, height) {
  ggsave(file.path(out, paste0(stem, ".pdf")), p, device = arial_pdf,
         width = width, height = height, bg = "white")
  ggsave(file.path(out, paste0(stem, ".png")), p, device = ragg::agg_png,
         width = width, height = height, dpi = 400, bg = "white")
}

panels <- lapply(seq_along(cohorts), function(i) make_panel(cohorts[i], LETTERS[i + 1L]))
row <- plot_grid(plotlist = panels, nrow = 1)
combined <- ggdraw() + draw_plot(row, x = 0, y = .15, width = 1, height = .85) +
  draw_grob(legend, x = 0, y = .045, width = 1, height = .09) +
  draw_label("Dashed lines: pooled medians; quadrants summarize continuous scores.",
               x = .5, y = .026, fontfamily = font, size = 8.5, colour = "#555555")
save_pair(combined, "Fig1_BCD_Cook_EMT_FigureYa175_Arial", 7.5, 3.35)
for (i in seq_along(cohorts)) {
  p <- make_panel(cohorts[i], LETTERS[i + 1L], standalone = TRUE)
  save_pair(p, paste0("Fig1", LETTERS[i + 1L], "_Cook_EMT_", keys[i], "_Arial"), 3.75, 3.65)
}

fwrite(summary, file.path(out, "quadrant_proportions_verified.csv"))
fwrite(data.table(state = states, colour = unname(cols), font = font), file.path(out, "style_palette.csv"))
audit <- data[, .(cells = .N, outside_display_range = sum(
  epithelial_score < xrange[1] | epithelial_score > xrange[2] |
  hallmark_EMT_score < yrange[1] | hallmark_EMT_score > yrange[2])), by = cohort]
audit[, `:=`(epithelial_cutoff = xcut, EMT_cutoff = ycut,
             x_min = xrange[1], x_max = xrange[2], y_min = yrange[1], y_max = yrange[2],
             scores_recomputed = FALSE, source_md5 = unname(tools::md5sum(source_file)))]
fwrite(audit, file.path(out, "plot_audit.csv"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
message("Completed B-D: ", out)
