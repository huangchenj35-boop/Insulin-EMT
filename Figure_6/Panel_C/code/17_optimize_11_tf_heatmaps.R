#!/usr/bin/env Rscript
# FigureYa-style optimization of panel 11 (core epithelial vs mesenchymal TF RNA).
# Uses the current frozen objects and the existing core-TF gene sets; no values
# or cell assignments are changed by this script.
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dorothea)
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
  Epithelial = c("GRHL2", "OVOL1", "OVOL2", "ELF3", "HNF4A", "HNF1A", "HNF1B", "SOX9", "KLF5", "GATA6"),
  Mesenchymal = c("SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2", "PRRX1", "PRRX2", "FOXC2", "TCF4", "FOSL1", "RUNX2")
)
allg <- intersect(unique(unlist(sets)), rownames(cnt))
logc <- log1p(t(t(cnt[allg, , drop = FALSE]) / lib) * 1e4)
role <- factor(as.character(obj$frozen_role), levels = groups)
gene_mean <- sapply(groups, function(g) rowMeans(logc[, role == g, drop = FALSE]))
colnames(gene_mean) <- unname(labels[groups])
rownames(gene_mean) <- allg
z <- t(scale(t(gene_mean))); z[!is.finite(z)] <- 0
program <- rep(NA_character_, length(allg)); names(program) <- allg
program[intersect(sets$Epithelial, allg)] <- "Epithelial core TF"
program[intersect(sets$Mesenchymal, allg)] <- "Mesenchymal core TF"
program <- factor(program[rownames(z)], levels = c("Epithelial core TF", "Mesenchymal core TF"))
stopifnot(!anyNA(program), identical(colnames(z), unname(labels[groups])))

# Report the small normal-cholangiocyte group transparently.  In USA-HCC this
# group contains only four frozen cells, so a zero SOX9 signal is a data
# availability/measurement result rather than evidence that cholangiocytes
# lack this canonical marker.
chol_n <- sum(role == "normal_cholangiocyte_root")
if (tag == "USAHCC") {
  qc_note <- "Normal cholangiocytes excluded from USA-HCC display (available n=4)."
} else {
  chol_sox9_detected <- if ("SOX9" %in% rownames(cnt)) {
    sum(cnt["SOX9", role == "normal_cholangiocyte_root"] > 0)
  } else 0L
  qc_note <- sprintf("Normal cholangiocyte: n=%d; SOX9 detected in %d/%d cells.",
                     chol_n, chol_sox9_detected, chol_n)
}
fwrite(as.data.table(z, keep.rownames = "gene"),
       file.path(OUTT, sprintf("fig2_%s_TF_expression_optimized_zscore.csv", stem_tag)))

# Infer activity for the same core TFs using signed DoRothEA A-C target sets.
# This is a target-expression proxy, not a measurement of TF protein activity.
net <- as.data.table(dorothea_hs)[confidence %in% c("A", "B", "C"), .(tf, target, mor)]
target_genes <- intersect(unique(net$target), rownames(cnt))
target_logc <- log1p(t(t(cnt[target_genes, , drop = FALSE]) / lib) * 1e4)
target_mean <- sapply(groups, function(g) rowMeans(target_logc[, role == g, drop = FALSE]))
rownames(target_mean) <- target_genes
target_z <- t(scale(t(target_mean))); target_z[!is.finite(target_z)] <- 0
activity <- matrix(NA_real_, nrow = length(allg), ncol = length(groups),
                   dimnames = list(allg, unname(labels[groups])))
for (tf_name in allg) {
  e <- net[tf == tf_name & target %in% rownames(target_z)]
  if (nrow(e) == 0) next
  idx <- match(e$target, rownames(target_z))
  activity[tf_name, ] <- colMeans(target_z[idx, , drop = FALSE] * e$mor)
}
activity <- activity[rownames(z), , drop = FALSE]
fwrite(as.data.table(activity, keep.rownames = "TF"),
       file.path(OUTT, sprintf("fig2_%s_TF_activity_optimized.csv", stem_tag)))

role_display <- unname(labels[groups])
role_display_cols <- setNames(unname(role_cols[groups]), role_display)
column_labels <- c(
  normal_cholangiocyte_root = "Cholangiocyte",
  normal_hepatocyte_root = "Hepatocyte",
  malignant_epithelial_EMT_root = "Malignant epithelial",
  qHSC_root = "qHSC", activated_HSC_intermediate = "Activated HSC",
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
    title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)))
)
row_col <- c("Epithelial core TF" = "#4C78A8", "Mesenchymal core TF" = "#E64B35")
left <- rowAnnotation(Program = program, col = list(Program = row_col),
                      annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
                      show_annotation_name = TRUE, width = unit(5, "mm"),
                      annotation_legend_param = list(Program = list(
                        title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8))))
cell_mm <- 7
body_width_mm <- ncol(z) * cell_mm
body_height_mm <- nrow(z) * cell_mm + 1
ht_expr <- Heatmap(z, name = "Row z-score",
              col = colorRamp2(c(-2.5, 0, 2.5), c("#3C5488", "#F7F7F7", "#E64B35")),
              na_col = "#FFFFFF", cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = program, row_gap = unit(1, "mm"), row_title_rot = 0,
              row_title_gp = gpar(fontsize = 9, fontface = "bold"),
              left_annotation = left, top_annotation = ha, row_names_side = "left", column_names_side = "bottom",
              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8, fontface = "bold"),
              column_labels = unname(column_labels[groups]), column_names_rot = 90,
              column_names_max_height = unit(38, "mm"), width = unit(body_width_mm, "mm"),
              height = unit(body_height_mm, "mm"),
              column_title = "TF RNA expression",
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              heatmap_legend_param = list(title = "Expression\nrow z-score", at = c(-2, 0, 2),
                                          title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)))
act_col <- colorRamp2(c(-1.5, 0, 1.5), c("#3C5488", "#F7F7F7", "#E64B35"))
ht_act <- Heatmap(activity, name = "Activity",
                  col = act_col, na_col = "#EEEEEE", cluster_rows = FALSE, cluster_columns = FALSE,
                  row_split = program, row_gap = unit(1, "mm"), row_title_rot = 0,
                  row_title_gp = gpar(fontsize = 9, fontface = "bold"), show_row_names = TRUE,
                  row_names_gp = gpar(fontsize = 8),
                  top_annotation = ha, row_names_side = "left",
                  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
                  column_labels = unname(column_labels[groups]), column_names_rot = 90,
                  column_names_max_height = unit(38, "mm"), width = unit(ncol(z) * cell_mm, "mm"),
                  height = unit(body_height_mm, "mm"),
                  column_title = "DoRothEA A-C activity",
                  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                  heatmap_legend_param = list(title = "Inferred\nactivity", at = c(-1.5, 0, 1.5),
                                              title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)))
# Render expression and activity as two independent figures.
device_width_mm <- max(205, ncol(z) * cell_mm + 160)
device_height_mm <- body_height_mm + 68
plots <- list(
  TF_expression = list(ht = ht_expr, title = sprintf("%s: core TF RNA expression", cohort_name),
                       note = sprintf("Predefined cell roles in %s; TF RNA expression is not TF protein activity or nuclear localization. %s", cohort_name, qc_note)),
  TF_activity = list(ht = ht_act, title = sprintf("%s: DoRothEA A-C TF activity proxy", cohort_name),
                     note = sprintf("Predefined cell roles in %s; signed DoRothEA A-C target-expression proxy. Grey = no supported A-C targets; TF protein activity and nuclear localization are not measured. %s", cohort_name, qc_note))
)
for (plot_name in names(plots)) for (ext in c("pdf", "png")) {
  out <- file.path(OUTP, sprintf("fig2_%s_%s.%s", stem_tag, plot_name, ext))
  if (ext == "pdf") pdf(out, width = device_width_mm / 25.4, height = device_height_mm / 25.4, useDingbats = FALSE)
  else png(out, width = round(device_width_mm / 25.4 * 300), height = round(device_height_mm / 25.4 * 300), res = 300)
  draw(plots[[plot_name]]$ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       column_title = plots[[plot_name]]$title, column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  note_text <- paste(strwrap(plots[[plot_name]]$note, width = 115), collapse = "\n")
  grid.text(note_text, x = unit(3, "mm"), y = unit(2, "mm"),
            just = c("left", "bottom"), gp = gpar(fontsize = 6.5, col = "#666666", lineheight = 1.05))
  dev.off()
}
cat(sprintf("[done] optimized panel 11 expression + activity figures: %s (%d cells, %d roles, %d TFs)\n", cohort_name, ncol(obj), length(groups), nrow(z)))
