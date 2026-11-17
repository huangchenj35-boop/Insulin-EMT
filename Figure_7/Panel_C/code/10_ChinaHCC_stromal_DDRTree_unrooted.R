#!/usr/bin/env Rscript
# Refit the unrooted Monocle-2 DDRTree for the stromal compartment only.
# This deliberately excludes the malignant epithelial root: the panel is the
# Figure-4 stromal topology (qHSC, activated HSC, portal fibroblast, fibrotic
# MFB) and is not a lineage-direction assignment.

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(monocle)
  library(ggplot2)
  library(data.table)
  library(grid)
})
set.seed(20260901)

ROOT <- "emt_invivo_bodyfig_20260901"
OUTF <- file.path(ROOT, "figures")
OUTP <- file.path("emt_invivo_bodyfig_standalone_20260901", "panels")
OUTT <- file.path("emt_invivo_bodyfig_standalone_20260901", "tables")
dir.create(OUTF, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTT, recursive = TRUE, showWarnings = FALSE)

tag <- Sys.getenv("EMT_DDRTREE_COHORT", "ChinaHCC")
stopifnot(tag %in% c("ChinaHCC", "USAHCC", "CCAcombined"))
source_cds <- if (tag == "ChinaHCC") {
  "results_v4/S6_ChinaHCC_joint_epithelial_sources_MFB_monocle2_topology.rds"
} else if (tag == "USAHCC") {
  "results_v4/S6_USAHCC_joint_epithelial_sources_MFB_monocle2_topology.rds"
} else {
  "results_v4/S6_CCAcombined_joint_sources_MFB_monocle2_topology.rds"
}
cds <- readRDS(source_cds)
stromal_roles <- if (tag == "ChinaHCC") {
  c("qHSC_root", "activated_HSC_intermediate", "portal_fibroblast_root", "fibrotic_MFB_terminal")
} else if (tag == "USAHCC") {
  c("qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal")
} else {
  c("qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal")
}
keep <- pData(cds)$frozen_role %in% stromal_roles
expected_n <- if (tag == "ChinaHCC") 377L else if (tag == "USAHCC") 223L else 157L
stopifnot(sum(keep) == expected_n)
stromal <- cds[, keep]

# Refit DDRTree on the 377 stromal cells so the displayed coordinates and MST
# are internally consistent after epithelial cells are removed.
stromal <- reduceDimension(stromal, max_components = 2, method = "DDRTree",
                           norm_method = "log", verbose = FALSE)
coords <- as.data.table(t(reducedDimS(stromal)), keep.rownames = "barcode")
setnames(coords, c("barcode", "DDRTree_1", "DDRTree_2"))
meta <- as.data.table(pData(stromal), keep.rownames = "barcode")
plot_data <- merge(coords, meta[, .(barcode, frozen_role)], by = "barcode", sort = FALSE)
stem_tag <- ifelse(tag == "ChinaHCC", "ChinaHCC", ifelse(tag == "USAHCC", "USAHCC", "CCAcombined"))
fwrite(plot_data, file.path(OUTT, sprintf("fig4a_%s_stromal_DDRTree_coordinates.csv", stem_tag)))
saveRDS(stromal, file.path(OUTP, sprintf("fig4a_%s_stromal_DDRTree_unrooted_cds.rds", stem_tag)))

role_labels <- if (tag == "ChinaHCC") {
  c(qHSC_root = "qHSC", activated_HSC_intermediate = "Activated HSC",
    portal_fibroblast_root = "Portal fibroblast", fibrotic_MFB_terminal = "Fibrotic MFB")
} else {
  c(qHSC_root = "qHSC", portal_fibroblast_like_root = "Portal fibroblast-like",
    fibrotic_MFB_terminal = "Fibrotic MFB")
}
role_cols <- c(qHSC_root = "#80C7D8", activated_HSC_intermediate = "#F2A07B",
               portal_fibroblast_root = "#43A77B", portal_fibroblast_like_root = "#43A77B",
               fibrotic_MFB_terminal = "#B65C85")

# plot_cell_trajectory draws the refitted principal tree and cell cloud.  The
# explicit palette/legend keeps all four stromal roles visible at publication
# scale and avoids implying a root or direction.
cohort_label <- ifelse(tag == "ChinaHCC", "China-HCC", ifelse(tag == "USAHCC", "USA-HCC", "CCA-combined"))
p <- plot_cell_trajectory(stromal, color_by = "frozen_role", show_tree = TRUE) +
  scale_color_manual(values = role_cols, breaks = names(role_labels),
                     labels = unname(role_labels), drop = FALSE) +
  labs(title = sprintf("A  %s stromal DDRTree (unrooted)", cohort_label),
       subtitle = sprintf("%d stromal cells; %s", expected_n, paste(unname(role_labels), collapse = ", ")),
       x = "DDRTree component 1", y = "DDRTree component 2", color = NULL,
       caption = "Principal tree is unrooted; topology does not assign lineage direction.") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 9, colour = "#555555"),
        plot.caption = element_text(size = 8, colour = "#555555", hjust = 0),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 9), legend.key.width = unit(12, "pt"))

stem_tag <- ifelse(tag == "ChinaHCC", "ChinaHCC", ifelse(tag == "USAHCC", "USAHCC", "CCAcombined"))
stem <- sprintf("fig4a_%s_stromal_DDRTree_unrooted", stem_tag)
ggsave(file.path(OUTP, paste0(stem, ".pdf")), p, width = 7, height = 7, device = cairo_pdf, bg = "white")
ggsave(file.path(OUTP, paste0(stem, ".png")), p, width = 7, height = 7, dpi = 300, bg = "white")
cat(sprintf("[done] %s unrooted stromal DDRTree: %d cells\n", tag, expected_n))
