#!/usr/bin/env Rscript
# FigureYa-style optimization of panel 10 (epithelial vs mesenchymal phenotype).
# The gene sets and expression summaries are unchanged; this only improves
# annotation, ordering, and publication layout.
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

tag <- Sys.getenv("EMT_HEATMAP_COHORT", "ChinaHCC")
stopifnot(tag %in% c("ChinaHCC", "USAHCC", "CCAcombined"))
cohort_name <- if (tag == "ChinaHCC") "China-HCC" else if (tag == "USAHCC") "USA-HCC" else "CCA-combined"
stem_tag <- if (tag == "ChinaHCC") "ChinaHCC" else if (tag == "USAHCC") "USAHCC" else "CCAcombined"
ROOT <- "emt_invivo_bodyfig_20260901"
OUTP <- "emt_invivo_bodyfig_standalone_20260901/panels"
OUTT <- "emt_invivo_bodyfig_standalone_20260901/tables"
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTT, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(ROOT, sprintf("S2_%s_joint_with_normal_roles_frozen.rds", stem_tag)))
if (tag == "USAHCC") obj <- subset(obj, subset = frozen_role != "normal_cholangiocyte_root")
DefaultAssay(obj) <- "RNA"
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
lib <- Matrix::colSums(cnt); lib[lib <= 0] <- 1
groups <- if (tag == "ChinaHCC") {
  c("normal_cholangiocyte_root", "normal_hepatocyte_root", "malignant_epithelial_EMT_root",
    "qHSC_root", "activated_HSC_intermediate", "portal_fibroblast_root", "fibrotic_MFB_terminal")
} else if (tag == "USAHCC") {
  c("normal_hepatocyte_root", "malignant_epithelial_EMT_root",
    "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal")
} else {
  c("normal_cholangiocyte_root", "malignant_epithelial_EMT_root", "qHSC_root",
    "portal_fibroblast_like_root", "fibrotic_MFB_terminal")
}
labels <- c(normal_cholangiocyte_root = "Normal cholangiocyte",
            normal_hepatocyte_root = "Normal hepatocyte",
            malignant_epithelial_EMT_root = "Malignant epithelial",
            qHSC_root = "qHSC", activated_HSC_intermediate = "Activated HSC",
            portal_fibroblast_root = "Portal fibroblast",
            portal_fibroblast_like_root = "Portal fibroblast-like",
            fibrotic_MFB_terminal = "Fibrotic MFB")
role_cols <- c(normal_cholangiocyte_root = "#4C78A8", normal_hepatocyte_root = "#F2CF5B",
               malignant_epithelial_EMT_root = "#E08214", qHSC_root = "#80C7D8",
               activated_HSC_intermediate = "#F28E2B", portal_fibroblast_root = "#43A77B",
               portal_fibroblast_like_root = "#43A77B", fibrotic_MFB_terminal = "#B65C85")
sets <- list(
  Epithelial = c("EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1","OCLN","TJP1","DSP","JUP","CLDN3","CLDN4","CLDN7","TACSTD2"),
  Mesenchymal = c("VIM","FN1","CDH2","ITGA5","ITGB1","LGALS1","SPARC","COL1A1","COL1A2","COL3A1","DCN","LUM","ACTA2","TAGLN","FAP","POSTN","CTHRC1","COL11A1")
)
allg <- intersect(unique(unlist(sets)), rownames(cnt))
logc <- log1p(t(t(cnt[allg, , drop = FALSE]) / lib) * 1e4)
role <- factor(as.character(obj$frozen_role), levels = groups)
gene_mean <- sapply(groups, function(g) rowMeans(logc[, role == g, drop = FALSE]))
colnames(gene_mean) <- unname(labels[groups])
rownames(gene_mean) <- allg
z <- t(scale(t(gene_mean))); z[!is.finite(z)] <- 0
program <- rep(NA_character_, length(allg)); names(program) <- allg
program[intersect(sets$Epithelial, allg)] <- "Epithelial"
program[intersect(sets$Mesenchymal, allg)] <- "Mesenchymal"
program <- factor(program[rownames(z)], levels = c("Epithelial", "Mesenchymal"))
stopifnot(!anyNA(program), identical(colnames(z), unname(labels[groups])))

fwrite(as.data.table(z, keep.rownames = "gene"),
       file.path(OUTT, sprintf("fig2_%s_phenotype_optimized_zscore.csv", stem_tag)))
role_display <- unname(labels[groups])
role_display_cols <- setNames(unname(role_cols[groups]), role_display)
# Keep the full role names in the legend, while using compact vertical labels
# under the matrix so long names cannot collide between columns.
column_labels <- c(
  normal_cholangiocyte_root = "Cholangiocyte",
  normal_hepatocyte_root = "Hepatocyte",
  malignant_epithelial_EMT_root = "Malignant epithelial",
  qHSC_root = "qHSC",
  activated_HSC_intermediate = "Activated HSC",
  portal_fibroblast_root = "Portal fibroblast",
  portal_fibroblast_like_root = "Portal fibroblast-like",
  fibrotic_MFB_terminal = "Fibrotic MFB"
)
ha <- HeatmapAnnotation(
  Role = factor(role_display, levels = role_display),
  col = list(Role = role_display_cols),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  annotation_height = unit(6, "mm"), show_annotation_name = FALSE,
  annotation_legend_param = list(Role = list(
    at = role_display, labels = role_display,
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)))
)
row_col <- c(Epithelial = "#4C78A8", Mesenchymal = "#E64B35")
left <- rowAnnotation(Program = program, col = list(Program = row_col),
                      annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
                      show_annotation_name = TRUE, width = unit(5, "mm"),
                      annotation_legend_param = list(Program = list(
                        title_gp = gpar(fontsize = 8, fontface = "bold"),
                        labels_gp = gpar(fontsize = 8))))
cell_mm <- 6
body_width_mm <- ncol(z) * cell_mm
body_height_mm <- nrow(z) * cell_mm + 1
ht <- Heatmap(z, name = "Row z-score",
              col = colorRamp2(c(-2.5, 0, 2.5), c("#3C5488", "#F7F7F7", "#E64B35")),
              na_col = "#FFFFFF", cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "#FFFFFF", lwd = 0.5),
              row_split = program, row_gap = unit(1, "mm"), row_title_rot = 0,
              row_title_gp = gpar(fontsize = 9, fontface = "bold"),
              left_annotation = left, top_annotation = ha, row_names_side = "left", column_names_side = "bottom",
              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7.5, fontface = "bold"),
              column_labels = c(normal_cholangiocyte_root = "Cholangiocyte",
                                normal_hepatocyte_root = "Hepatocyte",
                                malignant_epithelial_EMT_root = "Malignant epithelial",
                                qHSC_root = "qHSC",
                                activated_HSC_intermediate = "Activated HSC",
                                portal_fibroblast_root = "Portal fibroblast",
                                portal_fibroblast_like_root = "Portal fibroblast-like",
                                fibrotic_MFB_terminal = "Fibrotic MFB")[groups],
              column_names_rot = 90,
              column_names_max_height = unit(40, "mm"),
              width = unit(body_width_mm, "mm"), height = unit(body_height_mm, "mm"),
              column_title = sprintf("%s: epithelial versus mesenchymal phenotype", cohort_name),
              column_title_gp = gpar(fontsize = 11, fontface = "bold"),
              heatmap_legend_param = list(title = "Row z-score", at = c(-2, 0, 2),
                                          title_gp = gpar(fontsize = 8, fontface = "bold"),
                                          labels_gp = gpar(fontsize = 8)))
device_width_mm <- max(160, body_width_mm + 105)
device_height_mm <- body_height_mm + 52
for (ext in c("pdf", "png")) {
  out <- file.path(OUTP, sprintf("fig2_%s_phenotype.%s", stem_tag, ext))
  if (ext == "pdf") pdf(out, width = device_width_mm / 25.4, height = device_height_mm / 25.4, useDingbats = FALSE)
  else png(out, width = round(device_width_mm / 25.4 * 300), height = round(device_height_mm / 25.4 * 300), res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grid.text(sprintf("Predefined cell roles in %s; descriptive RNA expression (not protein localization).", cohort_name),
            x = unit(3, "mm"), y = unit(2, "mm"), just = c("left", "bottom"),
            gp = gpar(fontsize = 7, col = "#666666"))
  dev.off()
}
cat(sprintf("[done] optimized panel 10: %s (%d cells, %d roles)\n", cohort_name, ncol(obj), length(groups)))
