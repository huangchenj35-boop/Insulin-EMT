# v4 第三步: 间质细胞【重新高分辨率聚类】+ 逐cluster判定, 专门找 HSC 小群
suppressMessages({library(Seurat); library(harmony); library(ggplot2); library(dplyr)})
set.seed(1)

s2 <- readRDS("results_v3/S2c_CAF_subset.rds")   # 未clean全部CAF
cat("间质细胞数:", ncol(s2), " 数据集:", paste(unique(s2$dataset),collapse=","), "\n")

DefaultAssay(s2) <- "RNA"
s2 <- NormalizeData(s2, verbose=FALSE)
s2 <- FindVariableFeatures(s2, nfeatures=2000, verbose=FALSE)
s2 <- ScaleData(s2, verbose=FALSE)
s2 <- RunPCA(s2, npcs=30, verbose=FALSE)
s2 <- RunHarmony(s2, group.by.vars="dataset", verbose=FALSE)
s2 <- RunUMAP(s2, reduction="harmony", dims=1:30, verbose=FALSE)
s2 <- FindNeighbors(s2, reduction="harmony", dims=1:30, verbose=FALSE)
s2 <- FindClusters(s2, resolution=2.0, verbose=FALSE)
cat("高分辨率(res=2.0) cluster 数:", length(unique(s2$seurat_clusters)), "\n")

programs <- list(
  qHSC=c("LRAT","RELN","NGFR","GFAP","VIPR1","ANGPTL6","RBP1","CYGB","DES","HGF","COLEC11","PTH1R"),
  aHSC=c("ACTA2","COL1A1","COL3A1","PDGFRB","TAGLN","LOX","TIMP1","COL1A2"),
  portal=c("COL15A1","PI16","DPT","CFD","FBLN1","MFAP5","SLIT2","ELN","THY1","CD34"),
  pericyte=c("RGS5","CSPG4","MCAM","NOTCH3","KCNJ8","ABCC9","MYH11"),
  endo=c("PECAM1","VWF","CDH5","CLDN5","FLT1","KDR","ENG","EGFL7"),
  MFB=c("ACTA2","TAGLN","MYL9","CNN1","POSTN","FN1","TPM2"))
progs <- lapply(programs, function(g) intersect(g, rownames(s2)))
s2 <- AddModuleScore(s2, features=progs, name=paste0("MS_",names(programs)), seed=1)
sc <- paste0("MS_",names(programs),seq_along(programs)); names(sc) <- names(programs)

m <- s2@meta.data
agg <- m %>% group_by(seurat_clusters) %>% summarise(n=n(), across(all_of(sc), ~round(mean(.),2)))
colnames(agg) <- c("cluster","n",names(programs))
hscg <- intersect(c("LRAT","RELN","DES","CYGB","COLEC11"), rownames(s2))
ex <- as.data.frame(t(as.matrix(AverageExpression(s2, features=hscg, group.by="seurat_clusters", layer="data")$RNA)))
agg <- cbind(agg, round(ex[as.character(agg$cluster),],2))
write.csv(agg, "results_v4/tables/S3_recluster_module_by_cluster.csv", row.names=FALSE)
cat("\n=== 每 cluster: module score + HSC单基因 ===\n"); print(as.data.frame(agg))
cat("\n>>> qHSC score 最高的3个cluster:\n")
print(agg[order(-agg$qHSC)[1:3],])

saveRDS(s2, "results_v4/S3_stroma_reclustered.rds")
p1 <- DimPlot(s2, group.by="seurat_clusters", label=TRUE) + ggtitle("Stroma reclustered res2.0")
ggsave("results_v4/figures/S3_stroma_recluster_UMAP.pdf", p1, width=8, height=7)
ggsave("results_v4/figures/S3_stroma_qHSC_score_UMAP.pdf", FeaturePlot(s2, features="MS_qHSC1"), width=7, height=6)
cat("\n已写: S3_recluster_module_by_cluster.csv, S3_stroma_reclustered.rds, UMAP\n")
