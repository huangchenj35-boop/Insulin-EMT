#!/usr/bin/env Rscript
# FigureYa311PAM50Heatmap adaptation: fixed paired patient columns.
# Existing patient proportions only. No clustering or row standardization.
suppressPackageStartupMessages({library(ComplexHeatmap);library(circlize);library(grid);library(jsonlite)})
a=commandArgs(trailingOnly=FALSE)
script=sub('^--file=','',a[grepl('^--file=',a)])
work=dirname(dirname(normalizePath(script)));root=dirname(work)
pal=fromJSON(file.path(root,'配色参考主刊_20261010/palette.json'))
d=read.csv(file.path(work,'inputs/patient_percentages.csv'))
stat=read.csv(file.path(work,'inputs/paired_statistics.csv'))
patients=sprintf('HCC%02d',3:10)
metrics=c('EM_high','M_high');pops=c('Non-tumor','Malignant')
W=595.2756;H=215;cell=36;left=156
colorfun=colorRamp2(c(0,100),c(pal$continuous$mid,pal$continuous$high),space='RGB')
quartz(type='pdf',file=file.path(work,'audit/S1A_heatmap_layer.pdf'),width=W/72,height=H/72,family='Arial',bg='white')
grid.newpage()
txt=function(s,x,y,size=7,bold=FALSE,col='#262626',just='left') {
  grid.text(s,unit(x,'bigpts'),unit(H-y,'bigpts'),just=just,
    gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
}
rect=function(x,y,w,h,col) {
  grid.rect(unit(x+w/2,'bigpts'),unit(H-y-h/2,'bigpts'),width=unit(w,'bigpts'),height=unit(h,'bigpts'),gp=gpar(fill=col,col=NA))
}
txt('A',14,18,12,TRUE)
txt('Epithelial and mesenchymal states in paired HCC tissues',32,18,10,TRUE)
txt('HCC-GSE149614 | 8 paired patients',W/2,35,8,TRUE,just='centre')
for(j in seq_along(patients))txt(patients[j],left+(j-.5)*cell,48,7,just='centre')
txt('Paired P\n(Holm)',484,50,7,just='centre')
txt('Cells (%)',551,50,7,just='centre')
annotations=data.frame()
for(i in seq_along(metrics)) {
  top=c(58,138)[i]
  mat=matrix(NA_real_,nrow=2,ncol=8,dimnames=list(pops,patients))
  for(k in 1:2)for(j in 1:8) {
    q=d[d$metric==metrics[i]&d$patient==patients[j]&d$population==pops[k],]
    stopifnot(nrow(q)==1)
    mat[k,j]=q$percentage
    annotations=rbind(annotations,data.frame(metric=metrics[i],patient=patients[j],population=pops[k],
      percentage=q$percentage,label=sprintf('%.1f%%',q$percentage),
      x=left+(j-.5)*cell,y=top+(k-.5)*cell,width=cell,height=cell))
  }
  stopifnot(all(is.finite(mat)),all(mat>=0),all(mat<=100))
  ht=Heatmap(mat,name='Cells (%)',col=colorfun,
    cluster_rows=FALSE,cluster_columns=FALSE,show_row_names=FALSE,show_column_names=FALSE,
    show_heatmap_legend=FALSE,use_raster=FALSE,
    width=unit(8*cell,'bigpts'),height=unit(2*cell,'bigpts'),
    rect_gp=gpar(col='white',lwd=.7),
    cell_fun=function(j,k,x,y,width,height,fill) {
      grid.text(sprintf('%.1f%%',mat[k,j]),x,y,gp=gpar(fontfamily='Arial',fontsize=7,col='#262626'))
    })
  pushViewport(viewport(x=unit(left+4*cell,'bigpts'),y=unit(H-top-cell,'bigpts'),
    width=unit(8*cell,'bigpts'),height=unit(2*cell,'bigpts')))
  draw(ht,newpage=FALSE,padding=unit(rep(0,4),'bigpts'))
  popViewport()
  for(k in 1:2) {
    y=top+(k-.5)*cell
    txt(pops[k],144,y-5,7,just='right')
    txt(c('n = 1,575','n = 11,153')[k],144,y+7,7,col=pal$neutral$muted,just='right')
    rect(150,top+(k-1)*cell+.35,3,cell-.7,
         if(k==1)pal$cell_types$hepatocyte else pal$cell_types$malignant)
  }
  if(i==1) {
    txt('E/M high',14,top+27,8,TRUE)
    txt('Both scores high',14,top+41,7,col=pal$neutral$muted)
  } else {
    txt('M-high total',14,top+22,8,TRUE)
    txt('Mesenchymal high',14,top+36,7,col=pal$neutral$muted)
    txt('+ E/M high',14,top+47,7,col=pal$neutral$muted)
  }
  q=stat[stat$metric==metrics[i],]
  stopifnot(nrow(q)==1,q$n_pairs==8)
  txt(sprintf('%.3f',q$p_holm),484,top+cell,8,just='centre')
}
# Shared, unscaled percentage legend: white is 0%, coral is 100%.
lx=547;ly=73;lh=116;lw=8
grid.raster(as.raster(matrix(colorfun(seq(100,0,length.out=512)),ncol=1)),
  x=unit(lx+lw/2,'bigpts'),y=unit(H-ly-lh/2,'bigpts'),
  width=unit(lw,'bigpts'),height=unit(lh,'bigpts'),interpolate=TRUE)
for(v in seq(0,100,25))txt(as.character(v),562,ly+lh*(1-v/100),7,col=pal$neutral$muted)
dev.off()
write.csv(annotations,file.path(work,'audit/heatmap_annotations.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(work,'audit/R_sessionInfo.txt'))
cat('32 square 36-point tiles; common 0-100% scale; paired patient order retained\n')
