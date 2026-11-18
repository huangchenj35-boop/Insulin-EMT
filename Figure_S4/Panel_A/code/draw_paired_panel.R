#!/usr/bin/env Rscript
# Adapted from FigureYa195PanPaire (paired segments and points).
# Saved patient proportions only; no expression scoring is run.
suppressPackageStartupMessages({library(grid); library(jsonlite)})
arg <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', arg[grepl('^--file=', arg)])
work <- dirname(dirname(normalizePath(script)))
root <- dirname(work)
pal <- fromJSON(file.path(root, '配色参考主刊_20261010/palette.json'))
d <- read.csv(file.path(work, 'inputs/patient_percentages.csv'))
metrics <- c('EM_high','M_high')
stats <- do.call(rbind, lapply(metrics, function(k) {
  x <- d[d$metric==k & d$population=='Non-tumor', c('patient','percentage')]
  y <- d[d$metric==k & d$population=='Malignant', c('patient','percentage')]
  w <- merge(x,y,by='patient',suffixes=c('_non_tumor','_malignant'))
  stopifnot(nrow(w)==8, !anyDuplicated(w$patient))
  tt <- wilcox.test(w$percentage_malignant, w$percentage_non_tumor,
                    paired=TRUE, alternative='two.sided', exact=TRUE)
  data.frame(metric=k, n_pairs=8, statistic=unname(tt$statistic),
    p_raw=tt$p.value,
    median_patient_difference_pp=median(w$percentage_malignant-w$percentage_non_tumor),
    malignant_higher=sum(w$percentage_malignant>w$percentage_non_tumor))
}))
stats$p_holm <- p.adjust(stats$p_raw, method='holm')
write.csv(stats, file.path(work,'audit/paired_statistics.csv'), row.names=FALSE)

W <- 595.2756; H <- 215
quartz(type='pdf', file=file.path(work,'audit/S1A_paired_layer.pdf'),
       width=W/72, height=H/72, family='Arial', bg='white')
grid.newpage()
txt <- function(label,x,y,size=7,bold=FALSE,col='#262626',just='left',rot=0) {
  grid.text(label,x=unit(x,'bigpts'),y=unit(H-y,'bigpts'),just=just,rot=rot,
            gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
}
seg <- function(x1,y1,x2,y2,col='#D0D0D0',lwd=.65) {
  grid.segments(unit(x1,'bigpts'),unit(H-y1,'bigpts'),unit(x2,'bigpts'),unit(H-y2,'bigpts'),
                gp=gpar(col=col,lwd=lwd))
}
dot <- function(x,y,col,r=2.8) {
  grid.circle(unit(x,'bigpts'),unit(H-y,'bigpts'),r=unit(r,'bigpts'),
              gp=gpar(fill=col,col='white',lwd=.55))
}
txt('A',14,18,12,TRUE)
txt('Epithelial and mesenchymal states in paired HCC tissues',32,18,10,TRUE)
txt('HCC-GSE149614 | 8 paired patients',W/2,35,8,TRUE,just='centre')
titles <- c('E/M high','M-high total')
definitions <- c('E high and M high','Mesenchymal high + E/M high')
point_audit <- data.frame()
for (i in seq_along(metrics)) {
  mid <- c(164,447)[i]
  xleft <- mid-63; xright <- mid+51
  axisx <- mid-91
  y <- function(p) 174-p*1.04
  txt(titles[i],mid,52,8,TRUE,just='centre')
  txt(definitions[i],mid,64,7,col=pal$neutral$muted,just='centre')
  # A common 0-100% scale across both metrics.
  for (tick in seq(0,100,25)) {
    seg(axisx,y(tick),mid+89,y(tick),col=if(tick==0)'#666666' else '#E6E6E6',lwd=.45)
    txt(as.character(tick),axisx-7,y(tick),7,col=pal$neutral$muted,just='right')
  }
  x <- d[d$metric==metrics[i] & d$population=='Non-tumor',]
  z <- d[d$metric==metrics[i] & d$population=='Malignant',]
  w <- merge(x,z,by='patient',suffixes=c('_non_tumor','_malignant'))
  w <- w[order(w$patient),]
  # Shared deterministic horizontal offsets keep close endpoints visible
  # while preserving the same patient's pairing and exact percentage.
  offsets <- seq(-3.5,3.5,length.out=nrow(w))
  for(j in seq_len(nrow(w))) {
    seg(xleft+offsets[j],y(w$percentage_non_tumor[j]),
        xright+offsets[j],y(w$percentage_malignant[j]),col='#BBBBBB',lwd=.7)
  }
  for(j in seq_len(nrow(w))) {
    dot(xleft+offsets[j],y(w$percentage_non_tumor[j]),pal$cell_types$hepatocyte)
    dot(xright+offsets[j],y(w$percentage_malignant[j]),pal$cell_types$malignant)
    for(pop in c('Non-tumor','Malignant')) {
      pct <- if(pop=='Non-tumor')w$percentage_non_tumor[j] else w$percentage_malignant[j]
      px <- if(pop=='Non-tumor')xleft+offsets[j] else xright+offsets[j]
      point_audit <- rbind(point_audit,data.frame(metric=metrics[i],patient=w$patient[j],population=pop,percentage=pct,x=px,y=y(pct)))
    }
  }
  txt('Non-tumor',xleft,185,7,just='centre')
  txt('Malignant',xright,185,7,just='centre')
  txt('1,575 cells',xleft,195,7,col=pal$neutral$muted,just='centre')
  txt('11,153 cells',xright,195,7,col=pal$neutral$muted,just='centre')
  txt(sprintf('Paired P (Holm) = %.3f', stats$p_holm[i]),mid,207,7,just='centre')
}
txt('Cells in state (%)',22,122,7,just='centre',rot=90)
dev.off()
write.csv(point_audit,file.path(work,'audit/plotted_points.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(work,'audit/R_sessionInfo.txt'))
print(stats)
