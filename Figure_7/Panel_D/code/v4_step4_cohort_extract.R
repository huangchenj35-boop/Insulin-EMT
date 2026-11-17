# v4 CNV matching 扩展: 提取所有可评估共存患者的四组 counts (患者内)
suppressMessages({library(Seurat); library(data.table)})

ov <- fread("results_v4/tables/S2_malignant_MFB_overlap.csv")
evalpts <- ov[malig_epi>=30 & MFB>=10]$patient
cat("可评估共存患者(malig>=30 & MFB>=10):", paste(evalpts,collapse=", "), "\n")

fc <- fread("results_v4/tables/S0_final_cell_labels.csv")
slim <- fread("results_v4/tables/v4_allcells_meta_slim.csv"); slim[,barcode:=gsub('"','',barcode)]
caf <- readRDS("results_v3/S3f_CAF_clean_annotated.rds")
cm <- caf@meta.data; cm$barcode <- rownames(cm); cm$dataset <- as.character(cm$dataset)
cm$pid <- ifelse(grepl("151530",cm$dataset), as.character(cm$S_ID), as.character(cm$orig.ident))

cat("载入 19G counts ...\n"); flush.console()
seu <- readRDS("results_v3/scRNA_harmony_allcells_curate_v1.rds")
counts <- GetAssayData(seu, assay="RNA", layer="counts")
go <- readRDS("results_v3/cnv/GSE149614_HBV_infercnv/infercnv_obj.rds")@gene_order

ann_all <- list()
for (pt in evalpts) {
  me <- fc[patient==pt & grepl("^malignant",final_cell_label)]$barcode
  mfb <- cm$barcode[cm$pid==pt & cm$annotation=="MFB_contractile"]
  oc  <- cm$barcode[cm$pid==pt & cm$annotation!="MFB_contractile"]
  im  <- slim[patient_id==pt & celltype %in% c("T cell","B cell")]$barcode
  for (grp in c("malignant_epi","MFB","other_CAF","immune_ref")) {
    cells <- switch(grp, malignant_epi=me, MFB=mfb, other_CAF=oc, immune_ref=im)
    cells <- intersect(cells, colnames(counts))
    if (length(cells)) ann_all[[paste(pt,grp)]] <- data.frame(cell=cells, patient=pt, group=grp)
  }
}
ann <- do.call(rbind, ann_all)
mat <- as.matrix(counts[, unique(ann$cell)])
dir.create("results_v4/cnv_matching/cohort", recursive=TRUE, showWarnings=FALSE)
saveRDS(list(counts=mat, ann=ann, gene_order=go), "results_v4/cnv_matching/cohort/input.rds")
cat("\n各患者各组细胞数:\n"); print(table(ann$patient, ann$group))
cat("已存 cohort/input.rds\n")
