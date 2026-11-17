# 为 FigureYa257armCNV + FigureYa307CNVHeatmap 准备数据
# 底层患者内算 arm CNV, 按组(China-HCC/USA-HCC/CCA-combined)每患者等权汇总
suppressMessages({library(data.table)})
D <- readRDS("results_v4/cnv_matching/cohort/input.rds")
mat <- D$counts; ann <- as.data.table(D$ann); go <- D$gene_order
slim <- fread("results_v4/tables/v4_allcells_meta_slim.csv")
patient_group <- unique(slim[, .(patient_id, dataset_group)])
pg <- setNames(patient_group$dataset_group, patient_group$patient_id)
centro <- c(chr1=123.4,chr2=93.9,chr3=90.9,chr4=50,chr5=48.8,chr6=59.8,chr7=60.1,chr8=45.2,chr9=43,
 chr10=39.8,chr11=53.4,chr12=35.5,chr13=17.7,chr14=17.2,chr15=19,chr16=36.8,chr17=25.1,chr18=18.5,
 chr19=26.2,chr20=28.1,chr21=12,chr22=15,chrX=61,chrY=10.4)*1e6
genes <- intersect(rownames(mat), rownames(go)); g <- go[genes,]
g$arm <- paste0(g$chr, ifelse(g$start < centro[as.character(g$chr)],"p","q"))
arm_levels <- unique(g$arm[order(match(g$chr,paste0("chr",c(1:22,"X","Y"))),g$start)])

# 每患者: 每细胞 arm CNV (相对该患者免疫); 保存细胞级(供 FigureYa307) + pseudobulk(供 FigureYa257)
S4 <- fread("results_v4/tables/S4_cohort_cnv_matching.csv")
evalpts <- S4[call!="not_evaluable"]$patient    # 可评估: HCC06T,S042,S007
cellarm_list <- list(); pb_list <- list()
for (pt in evalpts) {
  ap <- ann[patient==pt]; m <- mat[genes, ap$cell]
  libs <- colSums(m); expr <- log2(sweep(m,2,libs,"/")*1e4+1)
  ref <- ap$cell[ap$group=="immune_ref"]
  rel <- expr - rowMeans(expr[,ref,drop=FALSE]); rel <- pmax(pmin(rel,3),-3)
  armmat <- rowsum(rel, g$arm); armmat <- armmat/as.numeric(table(g$arm)[rownames(armmat)])
  armmat <- armmat[arm_levels[arm_levels%in%rownames(armmat)],]
  grp <- setNames(ap$group, ap$cell)[colnames(armmat)]
  dsg <- pg[[pt]]
  # 细胞级 (FigureYa307) 只留 malignant+MFB 代表患者
  cellarm_list[[pt]] <- list(armmat=armmat, group=grp, dataset=dsg, patient=pt)
  # pseudobulk per celltype (每患者)
  for (ct in c("malignant_epi","MFB","other_CAF","immune_ref")) {
    cs <- names(grp)[grp==ct]; if(!length(cs)) next
    pb_list[[paste(pt,ct)]] <- data.table(patient=pt, dataset=dsg, celltype=ct,
      arm=rownames(armmat), cnv=rowMeans(armmat[,cs,drop=FALSE]))
  }
}
pb <- rbindlist(pb_list)
# 按组每患者等权平均 -> arm × celltype per group
grparm <- pb[, .(cnv=mean(cnv)), by=.(dataset,celltype,arm)]
grparm[, arm := factor(arm, levels=arm_levels)]
fwrite(grparm, "results_v4/tables/S4_armcnv_by_group_celltype.csv")
saveRDS(cellarm_list, "results_v4/cnv_matching/cohort/cellarm_list.rds")
cat("已写: S4_armcnv_by_group_celltype.csv (arm x celltype per group), cellarm_list.rds\n")
cat("\n=== 各组 malignant vs MFB 在关键arm的CNV (示例) ===\n")
print(dcast(grparm[celltype%in%c("malignant_epi","MFB") & arm%in%c("chr1q","chr8q","chr18q")],
      dataset+arm~celltype, value.var="cnv"))
