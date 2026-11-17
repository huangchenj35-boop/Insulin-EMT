#!/usr/bin/env Rscript
# Frozen MFB-focused PAGA ratios; presentation changes only.
# Observed/expected dot-plot design adapted from FigureYa317RoeDot:
# https://github.com/ying-ge/FigureYa/blob/main/FigureYa317RoeDot/FigureYa317RoeDot.Rmd
# Lu et al. (2025), iMetaMed, doi:10.1002/imm3.70005.
# Native R grid keeps type sizes identical in combined and individual exports.
suppressPackageStartupMessages({ library(data.table); library(grid) })
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
project <- normalizePath(file.path(dirname(script), "../.."))
out <- file.path(project, "EMT_Resutls/figures/MFB_connectivity_EFG")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tables <- file.path(project, "emt_invivo_bodyfig_20260901/tables")
inputs <- file.path(tables, c("fig1a_MFB_focused_ChinaHCC_values.csv",
  "fig5_USAHCC_current_joint_PAGA_edges.csv", "fig1a_MFB_focused_CCAcombined_values.csv"))
cohorts <- c("China-HCC", "USA-HCC", "CCA-combined")
keys <- c("ChinaHCC", "USAHCC", "CCAcombined")
target <- "fibrotic_MFB_terminal"
roles <- c("hepatocyte", "cholangiocyte", "malignant", "qHSC", "activated_HSC", "portal")
source_map <- c("Normal hepatocyte" = "hepatocyte", "Normal cholangiocyte" = "cholangiocyte",
  "Malignant epithelial" = "malignant", "qHSC" = "qHSC", "Activated HSC" = "activated_HSC",
  "Portal fibroblast" = "portal", "Portal fibroblast-like" = "portal")
usa_map <- c("normal_hepatocyte_root" = "hepatocyte", "malignant_epithelial_EMT_root" = "malignant",
  "qHSC_root" = "qHSC", "portal_fibroblast_like_root" = "portal")
china <- fread(inputs[1]); cca <- fread(inputs[3]); usa <- fread(inputs[2])
china[, role := unname(source_map[candidate])]
cca[, role := unname(source_map[candidate])]
usa <- usa[role_1 == target | role_2 == target]
usa[, candidate := fifelse(role_1 == target, role_2, role_1)]
usa[, `:=`(role = unname(usa_map[candidate]), ratio = paga_connectivity_k30_uncapped)]
keep <- c("cohort", "candidate", "role", "ratio")
dat <- rbindlist(list(china[, ..keep], usa[, ..keep], cca[, ..keep]))
stopifnot(nrow(dat) == 14L, !anyNA(dat), !anyDuplicated(dat[, .(cohort, role)]),
  all(dat$ratio > 0), all(dat$ratio < 3),
  identical(as.integer(table(factor(dat$cohort, levels = cohorts))), c(6L, 4L, 4L)))
dat[, `:=`(metric = "uncapped observed/expected PAGA", graph_k = 30L)]
dat[, source_file := inputs[match(cohort, cohorts)]]
dat[, source_md5 := unname(tools::md5sum(inputs))[match(cohort, cohorts)]]
dat[, display_value := sprintf("%.3f", ratio)]
setorder(dat, cohort, role)
fwrite(dat, file.path(out, "MFB_connectivity_frozen_values.csv"))

# Exact hues from the completed UMAP; related populations use distinct hues.
# Portal fibroblast-like and portal fibroblast are aligned for display only.
palette <- c(hepatocyte = "#E89B86", cholangiocyte = "#7895B2", malignant = "#F2B39A",
  qHSC = "#83C4D8", activated_HSC = "#72B9A8", portal = "#9CCDBD")
labels <- c(hepatocyte = "Non-malignant\nhepatocyte", cholangiocyte = "Non-malignant\ncholangiocyte",
  malignant = "Malignant\nepithelial", qHSC = "qHSC", activated_HSC = "Activated HSC",
  portal = "Portal fibroblast\n/ fibroblast-like")
font <- "Arial"; ink <- "#222222"; secondary <- "#555555"
style <- data.table(role = roles, label = unname(labels[roles]), colour = unname(palette[roles]),
  font = font, cohort_title_pt = 14, cell_label_pt = 11, ratio_pt = 10, axis_tick_pt = 10,
  panel_letter_pt = 16)
fwrite(style, file.path(out, "style_palette_and_font_sizes.csv"))

draw_chart <- function(indices, width, height, standalone = FALSE) {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, width), yscale = c(0, height)))
  txt <- function(label, x, y, size = 10, colour = ink, face = "plain", just = "left", ...) {
    grid.text(label, x, y, default.units = "native", just = just,
      gp = gpar(fontfamily = font, fontsize = size, col = colour, fontface = face, lineheight = .97), ...)
  }
  line <- function(x0, y0, x1, y1, colour = "#DDDFE2", lwd = .6, lty = "solid") {
    grid.segments(x0, y0, x1, y1, default.units = "native", gp = gpar(col = colour, lwd = lwd, lty = lty))
  }
  shown <- if (standalone) roles[roles %in% dat[cohort == cohorts[indices], role]] else roles
  n <- length(shown)
  row_top <- height - 1.07
  row_bottom <- .99
  yvals <- setNames(seq(row_top, row_bottom, length.out = n), shown)
  # Reserve a separate ratio column to the right of the actual log-scaled plot.
  label_left <- .10
  label_width <- if (standalone) 1.56 else 1.68
  gap <- if (standalone) 0 else .12
  panel_width <- (width - label_width - .10 - gap * (length(indices) - 1L)) / length(indices)
  header_y <- height - .71
  axis_y <- .71
  for (j in seq_along(shown)) {
    if (j %% 2 == 1L) grid.rect(x = width / 2, y = yvals[j], width = width - .13, height = .37,
      default.units = "native", gp = gpar(fill = "#F7F8FA", col = NA))
    label <- labels[shown[j]]
    if (standalone && shown[j] == "portal") label <- if (indices == 1L) "Portal fibroblast" else "Portal\nfibroblast-like"
    txt(label, label_left, yvals[j], size = 11)
  }
  txt("Cell population", label_left, header_y, size = 10, colour = secondary, face = "bold")
  for (j in seq_along(indices)) {
    i <- indices[j]
    start <- label_width + (j - 1L) * (panel_width + gap)
    plot_left <- start + .06
    plot_right <- start + panel_width - .46
    number_x <- start + panel_width - .025
    # Uniform domain, including the full uncapped range, across all cohorts.
    xpos <- function(value) plot_left + (log10(value) - log10(.001)) / (log10(3) - log10(.001)) * (plot_right - plot_left)
    txt(LETTERS[i + 4L], if (standalone) .10 else start - .025, height - .15, size = 16, face = "bold")
    title_x <- if (standalone) width / 2 else start + panel_width / 2 + .035
    txt(cohorts[i], title_x, height - .16, size = 14, face = "bold", just = "centre")
    txt("MFB connectivity", title_x, height - .43, size = 10, colour = secondary, just = "centre")
    txt("Ratio", number_x, header_y, size = 10, colour = secondary, face = "bold", just = "right")
    line(xpos(1), axis_y, xpos(1), row_top + .22, colour = "#8B9299", lwd = .8, lty = "dashed")
    line(plot_left, axis_y, plot_right, axis_y, colour = ink, lwd = .7)
    ticks <- c(.001, .01, .1, 1)
    for (k in seq_along(ticks)) {
      line(xpos(ticks[k]), axis_y, xpos(ticks[k]), axis_y - .045, colour = ink, lwd = .7)
      txt(c("0.001", "0.01", "0.1", "1")[k], xpos(ticks[k]), axis_y - .14, size = 10, just = "centre")
    }
    z <- dat[cohort == cohorts[i]]
    for (role_name in shown) {
      y <- yvals[role_name]
      value <- z[role == role_name, ratio]
      if (length(value) == 0L) {
        txt("-", number_x, y, size = 10, colour = "#8B9299", just = "right")
        next
      }
      # Thin horizontal segments compare a ratio with 1; they are not CIs or trajectories.
      line(xpos(value), y, xpos(1), y, colour = palette[role_name], lwd = 1.4)
      grid.points(xpos(value), y, default.units = "native", pch = 21, size = unit(3.1, "mm"),
        gp = gpar(fill = palette[role_name], col = "white", lwd = .75))
      txt(sprintf("%.3f", value), number_x, y, size = 10, just = "right")
    }
  }
  if (standalone) {
    txt("Uncapped PAGA observed/expected ratio (log scale)", width / 2, .32, size = 10, just = "centre")
    txt("Vertical reference: observed = expected", width / 2, .10, size = 9, colour = secondary, just = "centre")
  } else {
    txt("Uncapped PAGA observed/expected ratio (log scale)", width / 2, .32, size = 11, just = "centre")
    txt("Vertical reference: ratio = 1; dash: population not included in this comparison.",
      width / 2, .10, size = 9, colour = secondary, just = "centre")
  }
  popViewport()
}

save_chart <- function(stem, indices, width, height, standalone = FALSE) {
  grDevices::quartz(type = "pdf", file = file.path(out, paste0(stem, ".pdf")),
    width = width, height = height, family = font, bg = "white")
  draw_chart(indices, width, height, standalone)
  invisible(dev.off())
  ragg::agg_png(file.path(out, paste0(stem, ".png")), width = width, height = height, units = "in", res = 400,
    background = "white")
  draw_chart(indices, width, height, standalone)
  invisible(dev.off())
}
save_chart("Fig1_EFG_MFB_connectivity_FigureYa317_Arial", 1:3, 7.5, 4)
for (i in 1:3) save_chart(paste0("Fig1", LETTERS[i + 4L], "_MFB_connectivity_", keys[i], "_Arial"),
  i, 4.15, if (i == 1L) 4 else 3.5, standalone = TRUE)

audit <- dat[, .(candidate_populations = .N, min_ratio = min(ratio), max_ratio = max(ratio),
  source_file = source_file[1], source_md5 = source_md5[1]), by = cohort]
audit[, `:=`(ratios_recomputed = FALSE, roles_dropped = FALSE, ratios_capped = FALSE)]
fwrite(audit, file.path(out, "plot_audit.csv"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
message("Completed E-G: ", out)
