# v4 cohort CNV pattern matching (患者内) + 冻结规则判定
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

res <- list()
for (pt in unique(ann$patient)) {
  ap <- ann[patient==pt]; cells <- ap$cell; grp <- setNames(ap$group, ap$cell)
  m <- mat[genes, cells]; libs <- colSums(m); expr <- log2(sweep(m,2,libs,"/")*1e4 + 1)
  ref <- ap$cell[ap$group=="immune_ref"]
  rel <- expr - rowMeans(expr[,ref,drop=FALSE]); rel <- pmax(pmin(rel,3),-3)
  armmat <- rowsum(rel, g$arm); armmat <- armmat/as.numeric(table(g$arm)[rownames(armmat)])
  malig <- ap$cell[ap$group=="malignant_epi"]; mfb <- ap$cell[ap$group=="MFB"]
  consensus <- rowMeans(armmat[,malig,drop=FALSE])
  csum <- rowSums(armmat[,malig,drop=FALSE]); nM <- length(malig)
  malig_loo <- median(sapply(malig, function(c) suppressWarnings(cor(armmat[,c],(csum-armmat[,c])/(nM-1)))),na.rm=TRUE)
  key <- names(which(abs(consensus) > 0.15))
  cor_by <- function(cs) median(apply(armmat[,cs,drop=FALSE],2,function(v) suppressWarnings(cor(v,consensus))),na.rm=TRUE)
  bur_by <- function(cs) median(colMeans(abs(armmat[,cs,drop=FALSE])))
  mfb_pb <- rowMeans(armmat[,mfb,drop=FALSE])
  concord <- if(length(key)>0) mean(sign(mfb_pb[key])==sign(consensus[key])) else NA
  evaluable <- (length(malig)>=30 & length(mfb)>=10 & malig_loo>=0.6 & length(key)>=3)
  mfb_corr <- cor_by(mfb)
  retained <- if(!evaluable) "not_evaluable" else if(mfb_corr>=0.6 & !is.na(concord) & concord>=0.5) "possible_retention" else "not_retained"
  res[[pt]] <- data.table(patient=pt, dataset=pg[[pt]],
    malig_n=length(malig), MFB_n=length(mfb), immune_n=sum(ap$group=="immune_ref"),
    malig_LOO_corr=round(malig_loo,3), MFB_corr=round(mfb_corr,3),
    CAF_corr=round(cor_by(ap$cell[ap$group=="other_CAF"]),3), immune_corr=round(cor_by(ref),3),
    MFB_burden=round(bur_by(mfb),4), malig_burden=round(bur_by(malig),4),
    n_key_events=length(key), MFB_arm_concordance=round(concord,2), call=retained)
}
R <- rbindlist(res)[order(dataset,patient)]
fwrite(R, "results_v4/tables/S4_cohort_cnv_matching.csv")
cat("=== 每患者 CNV pattern matching (患者内) ===\n"); print(R)
cat("\n=== 按分组汇总 ===\n")
print(R[, .(n_patients=.N, n_not_retained=sum(call=="not_retained"),
            n_possible=sum(call=="possible_retention"), n_not_eval=sum(call=="not_evaluable"),
            median_MFB_corr=round(median(MFB_corr),3)), by=dataset])
cat("\n已写: S4_cohort_cnv_matching.csv\n")
