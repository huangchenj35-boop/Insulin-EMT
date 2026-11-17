#!/usr/bin/env Rscript
# Global UMAP before/after annotation. FigureYa93UMAP-inspired DimPlot core.
# Existing integrated UMAP coordinates are reused; no analysis is recomputed.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel); library(cowplot); library(grid); library(ggrastr)})
options(stringsAsFactors=FALSE); set.seed(20260910)
root <- "/home/sodajay/EMT"
obj_file <- file.path(root,"results_v3/scRNA_harmony_allcells_curate_v1.rds")
meta_file <- file.path(root,"results_v4/tables/v4_allcells_meta_slim.csv")
coord_file <- file.path(root,"results_v4/tables/S1_allcells_UMAP_coordinates_v4.csv.gz")
cache_file <- file.path(root,"EMT_Resutls/data/global_umap_plot_data.csv.gz")
out_dir <- file.path(root,"EMT_Resutls/figures/global_umap")
dir.create(dirname(cache_file),recursive=TRUE,showWarnings=FALSE); dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

if (!file.exists(cache_file)) {
  message("Loading integrated object once to cache cluster labels ...")
  seu <- readRDS(obj_file); stopifnot("umap" %in% Seurat::Reductions(seu))
  emb <- Seurat::Embeddings(seu,"umap"); md <- seu@meta.data; bc <- rownames(md)
  stopifnot("seurat_clusters" %in% colnames(md))
  d <- data.table(barcode=bc, UMAP_1=as.numeric(emb[bc,1]), UMAP_2=as.numeric(emb[bc,2]), cluster=as.character(md$seurat_clusters))
  ann <- fread(meta_file,select=c("barcode","celltype")); d <- merge(d,ann,by="barcode",all.x=TRUE,sort=FALSE)
  stopifnot(nrow(d)==95022L,!anyNA(d$celltype)); fwrite(d,cache_file)
  rm(seu,emb,md,d,ann); gc(verbose=FALSE)
}
d <- fread(cache_file); stopifnot(nrow(d)==95022L,uniqueN(d$barcode)==95022L)
d[, cluster := as.character(cluster)]

cell_levels <- c("Hepatocyte","Epithelial (cholangiocyte)","Fibroblast/CAF","Endothelial","Myeloid","T cell","B cell","Plasma cell")
cell_cols <- c("Hepatocyte"="#E89B86","Epithelial (cholangiocyte)"="#7895B2","Fibroblast/CAF"="#9CCDBD","Endothelial"="#F2B39A","Myeloid"="#72B9A8","T cell"="#83C4D8","B cell"="#A89078","Plasma cell"="#A3ADCA")
cl_levels <- unique(d$cluster)
cl_levels <- cl_levels[order(as.numeric(cl_levels))]
cl_cols <- setNames(grDevices::hcl(h=seq(15,375,length.out=length(cl_levels)+1)[seq_along(cl_levels)],c=55,l=72),cl_levels)
xlim <- range(d$UMAP_1); ylim <- range(d$UMAP_2); xlim <- xlim+c(-diff(xlim)*.025,diff(xlim)*.025); ylim <- ylim+c(-diff(ylim)*.025,diff(ylim)*.025)
th <- theme_void(base_family="Arial",base_size=22)+theme(plot.background=element_rect(fill="white",colour=NA),panel.background=element_rect(fill="white",colour=NA),plot.title=element_text(family="Arial",size=24,face="bold",hjust=.5,colour="#222222",margin=margin(b=6)),plot.subtitle=element_text(family="Arial",size=20,hjust=.5,colour="#555555",margin=margin(b=10)),legend.title=element_text(family="Arial",size=20,face="bold"),legend.text=element_text(family="Arial",size=20),legend.key.height=unit(.65,"cm"),legend.key.width=unit(.65,"cm"),plot.margin=margin(14,14,14,14))

make_before <- function(raster=FALSE, show_legend=TRUE){
  lab <- d[,.(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),by=cluster]
  pts <- if (raster) ggrastr::geom_point_rast(size=.16,alpha=.78,shape=16,raster.dpi=300) else geom_point(size=.16,alpha=.78,shape=16,stroke=0)
  ggplot(d,aes(UMAP_1,UMAP_2,colour=cluster))+pts+geom_label_repel(data=lab,aes(label=cluster),family="Arial",size=6.2,fontface="bold",colour="#333333",fill=scales::alpha("white",.84),label.size=0,box.padding=.3,point.padding=.1,min.segment.length=Inf,seed=20260910,show.legend=FALSE)+scale_colour_manual(values=cl_cols,breaks=cl_levels,name="Cluster")+coord_fixed(xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+labs(title="Global UMAP before annotation",subtitle="Unannotated Seurat clusters (95,022 cells)")+th+guides(colour=guide_legend(ncol=2,override.aes=list(size=4,alpha=1)))+if(!show_legend) theme(legend.position="none") else NULL
}
make_after <- function(raster=FALSE, show_legend=TRUE){
  lab <- d[,.(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),by=celltype]; lab[,celltype:=factor(celltype,levels=cell_levels)]
  pts <- if (raster) ggrastr::geom_point_rast(size=.16,alpha=.82,shape=16,raster.dpi=300) else geom_point(size=.16,alpha=.82,shape=16,stroke=0)
  ggplot(d,aes(UMAP_1,UMAP_2,colour=factor(celltype,levels=cell_levels)))+pts+geom_label_repel(data=lab,aes(label=celltype),family="Arial",size=6.2,fontface="bold",colour="#333333",fill=scales::alpha("white",.84),label.size=0,box.padding=.35,point.padding=.12,min.segment.length=Inf,seed=20260910,show.legend=FALSE)+scale_colour_manual(values=cell_cols,breaks=cell_levels,name=NULL,drop=FALSE)+coord_fixed(xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+labs(title="Global UMAP after annotation",subtitle="Eight audited major cell types (95,022 cells)")+th+guides(colour=guide_legend(ncol=1,override.aes=list(size=4,alpha=1)))+if(!show_legend) theme(legend.position="none") else NULL
}

# The paired figure deliberately removes legends from the map panels. This makes
# the two coordinate panels exactly the same physical size (8 x 8 inches each).
p_before <- make_before(raster=FALSE,show_legend=FALSE); p_after <- make_after(raster=FALSE,show_legend=FALSE)
p_before_r <- make_before(raster=TRUE,show_legend=FALSE); p_after_r <- make_after(raster=TRUE,show_legend=FALSE)
combo <- cowplot::plot_grid(p_before,p_after,ncol=2,rel_widths=c(1,1),align="hv",axis="tb")
combo_r <- cowplot::plot_grid(p_before_r,p_after_r,ncol=2,rel_widths=c(1,1),align="hv",axis="tb")

# Square map-only panels and equal-size paired panels.
ggsave(file.path(out_dir,"global_umap_before_annotation.pdf"),p_before,width=8,height=8,device=cairo_pdf,bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"global_umap_before_annotation.png"),p_before,width=8,height=8,dpi=300,bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"global_umap_after_annotation.pdf"),p_after,width=8,height=8,device=cairo_pdf,bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"global_umap_after_annotation.png"),p_after,width=8,height=8,dpi=300,bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"global_umap_before_after.pdf"),combo,width=16,height=8,device=cairo_pdf,bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"global_umap_before_after.png"),combo,width=16,height=8,dpi=300,bg="white",limitsize=FALSE)

# Lightweight PDFs: point clouds are rasterized at 300 dpi, while labels,
# titles, and axes remain vector text. These are much smaller than full-vector
# PDFs for 95,022 cells.
ggsave(file.path(out_dir,"global_umap_before_annotation_rasterized_points.pdf"),p_before_r,width=8,height=8,device=cairo_pdf,bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"global_umap_after_annotation_rasterized_points.pdf"),p_after_r,width=8,height=8,device=cairo_pdf,bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"global_umap_before_after_rasterized_points.pdf"),combo_r,width=16,height=8,device=cairo_pdf,bg="white",limitsize=FALSE)

# Individual legend files remain available for panel assembly.
legend_before <- cowplot::get_legend(make_before(raster=FALSE,show_legend=TRUE))
legend_after <- cowplot::get_legend(make_after(raster=FALSE,show_legend=TRUE))
ggsave(file.path(out_dir,"global_umap_before_annotation_legend.pdf"),legend_before,width=5,height=7,device=cairo_pdf,bg="white")
ggsave(file.path(out_dir,"global_umap_after_annotation_legend.pdf"),legend_after,width=4,height=7,device=cairo_pdf,bg="white")
message("Done: ",out_dir)
