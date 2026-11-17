#!/usr/bin/env Rscript
# Current-cell China-HCC heatmaps (all seven frozen roles, no legacy subset).
Sys.setenv(LANGUAGE="en"); options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(dorothea);library(ComplexHeatmap);library(circlize);library(grid)})
tag <- Sys.getenv("EMT_HEATMAP_COHORT", "ChinaHCC")
stopifnot(tag %in% c("ChinaHCC", "USAHCC", "CCAcombined"))
epithelial_only <- identical(Sys.getenv("EMT_METABOLISM_SCOPE", "all"), "epithelial")
file_suffix <- if (epithelial_only) "_epithelial_only" else ""
cohort_name <- if (tag == "ChinaHCC") "China-HCC" else if (tag == "USAHCC") "USA-HCC" else "CCA-combined"
stem_tag <- if (tag == "ChinaHCC") "ChinaHCC" else if (tag == "USAHCC") "USAHCC" else "CCAcombined"
ROOT <- "emt_invivo_bodyfig_20260901"; OUTP <- file.path("emt_invivo_bodyfig_standalone_20260901/panels", if (epithelial_only) "epithelial_only" else ""); OUTT <- "emt_invivo_bodyfig_standalone_20260901/tables"
dir.create(OUTP,recursive=TRUE,showWarnings=FALSE);dir.create(OUTT,recursive=TRUE,showWarnings=FALSE)
obj <- readRDS(file.path(ROOT, sprintf("S2_%s_joint_with_normal_roles_frozen.rds", stem_tag))); DefaultAssay(obj)<-"RNA"
if (tag == "USAHCC") obj <- subset(obj, subset = frozen_role != "normal_cholangiocyte_root")
cnt <- GetAssayData(obj,assay="RNA",layer="counts"); lib <- Matrix::colSums(cnt); lib[lib<=0]<-1
groups <- if (tag == "ChinaHCC") c("normal_cholangiocyte_root","normal_hepatocyte_root","malignant_epithelial_EMT_root","qHSC_root","activated_HSC_intermediate","portal_fibroblast_root","fibrotic_MFB_terminal") else if (tag == "USAHCC") c("normal_hepatocyte_root","malignant_epithelial_EMT_root","qHSC_root","portal_fibroblast_like_root","fibrotic_MFB_terminal") else c("normal_cholangiocyte_root","malignant_epithelial_EMT_root","qHSC_root","portal_fibroblast_like_root","fibrotic_MFB_terminal")
if (epithelial_only) groups <- if (tag == "ChinaHCC") c("normal_cholangiocyte_root", "normal_hepatocyte_root", "malignant_epithelial_EMT_root") else if (tag == "USAHCC") c("normal_hepatocyte_root", "malignant_epithelial_EMT_root") else c("normal_cholangiocyte_root", "malignant_epithelial_EMT_root")
labels <- c(normal_cholangiocyte_root="Normal cholangiocyte",normal_hepatocyte_root="Normal hepatocyte",malignant_epithelial_EMT_root="Malignant epithelial",qHSC_root="qHSC",activated_HSC_intermediate="Activated HSC",portal_fibroblast_root="Portal fibroblast",portal_fibroblast_like_root="Portal fibroblast-like",fibrotic_MFB_terminal="Fibrotic MFB")
role_cols <- c(normal_cholangiocyte_root="#4C78A8", normal_hepatocyte_root="#F2CF5B", malignant_epithelial_EMT_root="#E08214", qHSC_root="#80C7D8", activated_HSC_intermediate="#F28E2B", portal_fibroblast_root="#43A77B", portal_fibroblast_like_root="#43A77B", fibrotic_MFB_terminal="#B65C85")
role_display <- unname(labels[groups])
role_display_cols <- setNames(unname(role_cols[groups]), role_display)
ha <- HeatmapAnnotation(Role = factor(role_display, levels = role_display),
                        col = list(Role = role_display_cols), show_annotation_name = FALSE,
                        annotation_height = unit(6, "mm"),
                        annotation_legend_param = list(Role = list(at = role_display, labels = role_display,
                                                                   title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                                   labels_gp = gpar(fontsize = 8))))
sets <- list(Epithelial=c("EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1","OCLN","TJP1","DSP","JUP","CLDN3","CLDN4","CLDN7","TACSTD2"),Mesenchymal=c("VIM","FN1","CDH2","ITGA5","ITGB1","LGALS1","SPARC","COL1A1","COL1A2","COL3A1","DCN","LUM","ACTA2","TAGLN","FAP","POSTN","CTHRC1","COL11A1"),Epithelial_TF=c("GRHL2","OVOL1","OVOL2","ELF3","HNF4A","HNF1A","HNF1B","SOX9","KLF5","GATA6"),Mesenchymal_TF=c("SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","TWIST2","PRRX1","PRRX2","FOXC2","TCF4","FOSL1","RUNX2"),TGF=c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7","SERPINE1","CTGF","COL1A1"),INSULIN_AKT=c("INSR","IRS1","IRS2","PIK3CA","PIK3CB","AKT1","AKT2","MTOR","FOXO1","GSK3B","RPS6KB1"),HGF_AKT=c("HGF","MET","GAB1","PIK3CA","AKT1","MAPK1","MAPK3","SOS1","GRB2"),EGF_RAS=c("EGF","EGFR","GRB2","SOS1","KRAS","NRAS","HRAS","RAF1","BRAF","MAP2K1","MAPK1","ELK1"),Anabolic=c("FASN","ACACA","ACLY","SCD","HMGCR","HMGCS1","SQLE","ASNS","PHGDH","PSAT1","SHMT2","CAD","DHODH","RRM1","RRM2","EIF4E","RPS6KB1","MYC"),Catabolic=c("FOXO1","FOXO3","CPT1A","ACADM","ACOX1","PPARA","PNPLA2","LPL","LIPE","ULK1","ATG5","ATG7","BECN1","CTSB","CTSD","SQSTM1"))
allg <- unique(unlist(sets)); allg <- intersect(allg,rownames(cnt)); logc <- log1p(t(t(cnt[allg,,drop=FALSE])/lib)*1e4)
role <- as.character(obj$frozen_role); role <- factor(role,levels=groups)
gene_mean <- rbindlist(lapply(allg,function(g) data.table(gene=g,group=groups,mean_expression=sapply(groups,function(r) mean(logc[g, !is.na(role) & role == r])))))
fwrite(gene_mean,file.path(OUTT,sprintf("fig2_%s_current_gene_mean_by_role%s.csv", stem_tag, file_suffix)))
build <- function(prs, stem, title) {
  genes <- unique(unlist(sets[prs])); genes <- intersect(genes, allg)
  mat <- matrix(NA_real_, length(genes), length(groups),
                dimnames = list(genes, role_display))
  split <- setNames(character(length(genes)), genes)
  for (pr in prs) split[intersect(sets[[pr]], genes)] <- pr
  for (j in seq_along(groups)) {
    mat[, j] <- gene_mean[gene %in% genes & group == groups[j]]$mean_expression
  }
  z <- t(scale(t(mat))); z[!is.finite(z)] <- 0
  row_split_levels <- prs
  row_split <- split
  # Keep the original gene-expression rows and append one role-level balance score
  # to the metabolism plot so the relative anabolic/catabolic bias is visible in
  # the same figure. The score is computed by script 19 from the same frozen cells.
  if (stem == "metabolism") {
    score_file <- file.path(OUTT, sprintf("fig2_%s_metabolism_module_scores_zscore%s.csv", stem_tag, file_suffix))
    if (file.exists(score_file)) {
      score_dt <- fread(score_file)
      score_row <- as.numeric(score_dt[module == "Anabolic_minus_Catabolic", ..role_display])
      names(score_row) <- role_display
      z <- rbind(z, `Score: Anabolic - Catabolic` = score_row)
      row_split <- c(row_split, `Score: Anabolic - Catabolic` = "Module score")
      row_split_levels <- c(row_split_levels, "Module score")
    }
    # FOXO1 itself is shown as direct RNA expression (the FOXO1 gene is also
    # retained in the Catabolic gene block above). This is not protein activity.
    if ("FOXO1" %in% rownames(z)) {
      z <- rbind(z, `Expression: FOXO1 (direct RNA)` = z["FOXO1", ])
      row_split <- c(row_split, `Expression: FOXO1 (direct RNA)` = "Direct TF RNA")
      row_split_levels <- c(row_split_levels, "Direct TF RNA")
    }
    viper_file <- file.path(OUTT, sprintf("fig2_%s_FOXO1_VIPER_activity_zscore.csv", stem_tag))
    if (file.exists(viper_file)) {
      viper_dt <- fread(viper_file)
      viper_row <- as.numeric(viper_dt[metric == "FOXO1_VIPER_activity", ..role_display])
      names(viper_row) <- role_display
      z <- rbind(z, `Activity: FOXO1 VIPER` = viper_row)
      row_split <- c(row_split, `Activity: FOXO1 VIPER` = "TF activity")
      row_split_levels <- c(row_split_levels, "TF activity")
    }
  }
  fwrite(as.data.table(z, keep.rownames = "gene"),
         file.path(OUTT, paste0("fig2_", stem_tag, "_", stem, file_suffix, "_zscore.csv")))
  ht <- Heatmap(z, name = "Row z-score",
                col = colorRamp2(c(-2.5, 0, 2.5), c("#3C5488", "#F7F7F7", "#E64B35")),
                na_col = "#FFFFFF", cluster_rows = FALSE, cluster_columns = FALSE,
                rect_gp = gpar(col = "#FFFFFF", lwd = 0.5),
                row_split = factor(row_split, levels = row_split_levels), row_gap = unit(1, "mm"),
                row_title_rot = 0, row_title_gp = gpar(fontsize = 9, fontface = "bold"),
                top_annotation = ha, row_names_side = "left", column_names_side = "bottom",
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 7.5, fontface = "bold"),
                column_names_rot = 90, column_names_max_height = unit(40, "mm"),
                width = unit(length(groups) * 6, "mm"),
                height = unit(length(genes) * 6 + 1, "mm"),
                column_title = title, column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                heatmap_legend_param = list(title = "Row z-score", at = c(-2, 0, 2),
                                            title_gp = gpar(fontsize = 8, fontface = "bold"),
                                            labels_gp = gpar(fontsize = 8)))
  h <- (length(genes) * 6 + 1 + 52) / 25.4
  w <- if (stem == "metabolism") max(180, length(groups) * 6 + 125) / 25.4 else max(160, length(groups) * 6 + 105) / 25.4
  note <- if (stem == "metabolism") {
    sprintf("Predefined cell roles in %s; gene and FOXO1 rows show direct RNA expression. The score row is relative anabolic-minus-catabolic; FOXO1 VIPER is inferred from signed DoRothEA A-C targets. Neither is a direct measurement of FOXO1 protein activity or localization.", cohort_name)
  } else {
    sprintf("Predefined cell roles in %s; descriptive RNA expression (not protein localization).", cohort_name)
  }
  for (ext in c("pdf", "png")) {
    out <- file.path(OUTP, paste0("fig2_", stem_tag, "_", stem, file_suffix, ".", ext))
    if (ext == "pdf") pdf(out, width = w, height = h, useDingbats = FALSE)
    else png(out, width = round(w * 300), height = round(h * 300), res = 300)
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
    grid.text(paste(strwrap(note, width = 105), collapse = "\n"), x = unit(3, "mm"), y = unit(2, "mm"),
              just = c("left", "bottom"), gp = gpar(fontsize = 6.5, col = "#666666", lineheight = 1.05))
    dev.off()
  }
}
build(c("Epithelial","Mesenchymal"),"phenotype",sprintf("%s: epithelial versus mesenchymal phenotype", cohort_name))
build(c("Epithelial_TF","Mesenchymal_TF"),"TF_expression",sprintf("%s: epithelial versus mesenchymal core TFs", cohort_name))
build(c("TGF","INSULIN_AKT","HGF_AKT","EGF_RAS"),"signaling",sprintf("%s: EMT-relevant signaling programs", cohort_name))
build(c("Anabolic","Catabolic"),"metabolism",sprintf("%s: anabolic versus catabolic metabolism", cohort_name))
cat(sprintf("[done] %s current-cell heatmaps (%d cells)\n", cohort_name, ncol(obj)))
