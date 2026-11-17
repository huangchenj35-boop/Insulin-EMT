#!/usr/bin/env Rscript
# Export existing frozen values for Supplement S1; no inference is rerun.
args <- commandArgs(trailingOnly=FALSE)
self <- sub("^--file=", "", args[grepl("^--file=", args)][1])
root <- normalizePath(file.path(dirname(self), "../.."))
out <- file.path(root, "EMT_Resutls/figures/Supplement_S1/inputs")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read <- function(p) read.csv(file.path(root,p), check.names=FALSE, stringsAsFactors=FALSE)
write <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na="")

# Slots are read as attributes so no Seurat installation is needed for this
# read-only extraction of already-computed embeddings and module scores.
s <- readRDS(file.path(root,"results_v4/S3_stroma_reclustered.rds"))
m <- attr(s,"meta.data"); e <- attr(attr(s,"reductions")$umap,"cell.embeddings")
stopifnot(nrow(m)==3018L,identical(rownames(e),rownames(m)))
d <- data.frame(barcode=rownames(m),UMAP_1=e[,1],UMAP_2=e[,2],
                cluster=as.character(m$seurat_clusters),qHSC_score=m$MS_qHSC1)
rolefiles <- c("results_v4/tables/S6_HBV_stromal_frozen_roles_by_cell.csv",
 "results_v4/tables/S6_USAHCC_stromal_frozen_roles_by_cell.csv",
 "results_v4/tables/S6_CCAcombined_stromal_frozen_roles_by_cell.csv")
roles <- do.call(rbind,lapply(seq_along(rolefiles),function(i){
 z<-read(rolefiles[i]); data.frame(barcode=z$barcode,frozen_role=z$frozen_role,
 cohort=c("China-HCC","USA-HCC","CCA-combined")[i])
}))
stopifnot(!anyDuplicated(roles$barcode))
d$frozen_role <- roles$frozen_role[match(d$barcode,roles$barcode)]
d$cohort <- roles$cohort[match(d$barcode,roles$barcode)]
write(d,"stromal_umap_frozen.csv")

write(read("results_v4/tables/S4_cohort_cnv_matching.csv"),"cnv_matching_all_samples.csv")
write(read("results_v4/tables/S4_armcnv_by_group_celltype.csv"),"arm_cnv_by_group.csv")

# Retain exactly the display cell IDs selected in the original final panels.
cnv <- readRDS(file.path(root,"results_v4/cnv_matching/cohort/cellarm_list.rds"))
audits <- c(HCC06T="results_v4/tables/S4_review_cellCNV_ChinaHCC_HCC06T_v4_audit.csv",
 S042="results_v4/tables/S4_review_cellCNV_USAHCC_S042_v4_audit.csv")
allcounts <- list()
for (id in names(audits)) {
 a <- read(audits[[id]]); x<-cnv[[id]]
 stopifnot(identical(colnames(x$armmat),names(x$group)),all(a$cell %in% colnames(x$armmat)))
 z<-x$armmat[,a$cell,drop=FALSE]
 write(data.frame(arm=rownames(z),z,check.names=FALSE),paste0("cellarm_",id,".csv"))
 write(a,paste0("display_cells_",id,".csv"))
 full<-table(x$group); shown<-table(x$group[a$cell])
 allcounts[[id]]<-data.frame(sample=id,group=names(full),full_n=as.integer(full),
   displayed_n=as.integer(shown[names(full)]))
}
write(do.call(rbind,allcounts),"cnv_display_counts.csv")
cat("S1 exported frozen embeddings/scores and original CNV display selections.\n")
print(table(d$frozen_role,useNA="ifany"));print(do.call(rbind,allcounts))
