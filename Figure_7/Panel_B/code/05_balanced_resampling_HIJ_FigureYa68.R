#!/usr/bin/env Rscript
# FigureYa68friendsV2 raincloud layout, adapted to uncapped balanced PAGA resamples.
# https://github.com/ying-ge/FigureYa/blob/main/FigureYa68friendsV2/FigureYa68friendsV2.Rmd
# First run 05a_rebuild_uncapped_balanced_PAGA.py. Half-violin + jitter + median/95% interval.
# R grid keeps the same physical Arial sizes in all exports.
suppressPackageStartupMessages({library(data.table); library(grid)})
argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
project <- normalizePath(file.path(dirname(script), "../.."))
out <- file.path(project, "EMT_Resutls/figures/Balanced_resampling_HIJ")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
keys <- c("ChinaHCC", "USAHCC", "CCAcombined")
cohorts <- c("China-HCC", "USA-HCC", "CCA-combined")
inputs <- file.path(out, "uncapped_data",
  sprintf("S4_%s_balanced_resampling_uncapped_by_iteration.csv", keys))
role_map <- c(normal_hepatocyte_root = "hepatocyte", normal_cholangiocyte_root = "cholangiocyte",
  malignant_epithelial_EMT_root = "malignant", qHSC_root = "qHSC",
  activated_HSC_intermediate = "activated_HSC", portal_fibroblast_root = "portal",
  portal_fibroblast_like_root = "portal")
roles <- c("hepatocyte", "cholangiocyte", "malignant", "qHSC", "activated_HSC", "portal")
dat <- rbindlist(lapply(seq_along(inputs), function(i) {
  z <- fread(inputs[i]); stopifnot(all(z$cohort == cohorts[i]))
  z[, `:=`(source_file = inputs[i], source_md5 = unname(tools::md5sum(inputs[i])))]
  z
}))
dat[, role := unname(role_map[sub("--MFB$", "", edge)])]
dat[, paga_connectivity := paga_observed_expected_ratio]
stopifnot(nrow(dat) == 5600L, !anyNA(dat), all(dat$n_neighbors %in% c(10L, 20L)),
  all(dat$paga_connectivity >= 0), any(dat$paga_connectivity > 1),
  !anyDuplicated(dat[, .(cohort, n_neighbors, replicate, role)]))
counts <- dat[, .(n = .N, repeats = uniqueN(replicate)), by = .(cohort, n_neighbors, role)]
stopifnot(all(counts$n == 200L), all(counts$repeats == 200L))
for (i in 1:3) {
  z <- dat[cohort == cohorts[i]]
  stopifnot(uniqueN(z$role) == c(6,4,4)[i], identical(unique(z$n_per_role), c(29L,59L,16L)[i]),
    identical(sort(unique(z$replicate)), 1:200))
}

summary <- dat[, .(n_resamples = .N, n_per_role = unique(n_per_role),
  median = median(paga_connectivity), q025 = as.numeric(quantile(paga_connectivity, .025)),
  q975 = as.numeric(quantile(paga_connectivity, .975)),
  minimum = min(paga_connectivity), maximum = max(paga_connectivity),
  n_at_zero = sum(paga_connectivity == 0), n_above_one = sum(paga_connectivity > 1)),
  by = .(cohort, n_neighbors, role, edge)]
# Verify that re-capping the recomputed values reproduces upstream summary statistics.
reference <- rbindlist(lapply(keys, function(key) fread(file.path(project,
  "emt_invivo_bodyfig_20260901/tables", sprintf("S4_%s_balanced_resampling_edge_summary.csv", key)))))
capped_summary <- dat[, .(median=median(pmin(paga_connectivity,1)),
  q025=as.numeric(quantile(pmin(paga_connectivity,1),.025)),
  q975=as.numeric(quantile(pmin(paga_connectivity,1),.975))),by=.(cohort,n_neighbors,edge)]
comparison <- merge(capped_summary, reference, by = c("cohort", "n_neighbors", "edge"), suffixes = c("", "_saved"))
stopifnot(nrow(comparison) == 28L, all(abs(comparison$median-comparison$median_saved)<1e-12),
  all(abs(comparison$q025-comparison$q025_saved)<1e-12), all(abs(comparison$q975-comparison$q975_saved)<1e-12))
paired <- dcast(dat, cohort + n_neighbors + replicate + n_per_role ~ role, value.var = "paga_connectivity")
paired[, epithelial_max := pmax(hepatocyte, cholangiocyte, malignant, na.rm = TRUE)]
paired[, `:=`(portal_gt_epithelial = portal > epithelial_max + 1e-12,
  portal_gt_qHSC = portal > qHSC + 1e-12, portal_ties_qHSC = abs(portal-qHSC)<1e-12)]
stability <- paired[, .(n_resamples = .N, n_per_role = unique(n_per_role),
  portal_gt_epithelial = sum(portal_gt_epithelial), portal_gt_qHSC = sum(portal_gt_qHSC),
  portal_ties_qHSC = sum(portal_ties_qHSC)), by = .(cohort, n_neighbors)]
stopifnot(all(stability[n_neighbors==20, portal_gt_epithelial] == 200L),
  stability[cohort=="China-HCC" & n_neighbors==10, portal_gt_epithelial] == 185L)
fwrite(dat, file.path(out, "balanced_resampling_uncapped_values.csv"))
fwrite(summary, file.path(out, "resampling_summary_verified.csv"))
fwrite(stability, file.path(out, "paired_stability_by_cohort_and_k.csv"))
fwrite(paired, file.path(out, "paired_resample_comparisons.csv"))

palette <- c(hepatocyte = "#E89B86", cholangiocyte = "#7895B2", malignant = "#F2B39A",
  qHSC = "#83C4D8", activated_HSC = "#72B9A8", portal = "#9CCDBD")
labels <- c(hepatocyte="Non-malignant\nhepatocyte", cholangiocyte="Non-malignant\ncholangiocyte",
  malignant="Malignant\nepithelial", qHSC="qHSC", activated_HSC="Activated HSC",
  portal="Portal fibroblast\n/ fibroblast-like")
font <- "Arial"; ink <- "#222222"; grey <- "#555555"
xmax <- ceiling(max(dat$paga_connectivity)*2)/2
fwrite(data.table(role=roles, colour=unname(palette[roles]), font=font,
  title_pt=14, label_pt=11, tick_pt=10, panel_letter_pt=16), file.path(out,"style_palette_and_font_sizes.csv"))

draw_chart <- function(indices, k, width, height, standalone=FALSE, supplementary=FALSE) {
  grid.newpage(); pushViewport(viewport(xscale=c(0,width), yscale=c(0,height)))
  txt <- function(label,x,y,size=10,colour=ink,face="plain",just="left") {
    grid.text(label,x,y,default.units="native",just=just,
      gp=gpar(fontfamily=font,fontsize=size,col=colour,fontface=face,lineheight=.96))
  }
  line <- function(x0,y0,x1,y1,colour="#D8DCE0",lwd=.5,lty="solid") {
    grid.segments(x0,y0,x1,y1,default.units="native",gp=gpar(col=colour,lwd=lwd,lty=lty))
  }
  shown <- if(standalone) roles[roles %in% dat[cohort==cohorts[indices],role]] else roles
  yvals <- setNames(seq(height-1.05,1.42,length.out=length(shown)),shown)
  label_width <- if(standalone) 1.57 else 1.69
  gap <- if(standalone) 0 else .16
  panel_width <- (width-label_width-.09-gap*(length(indices)-1L))/length(indices)
  axis_y <- 1.12
  txt("Cell population",.10,height-.72,size=10,colour=grey)
  for(role_name in shown) {
    lab <- labels[role_name]
    if(standalone && role_name=="portal") lab <- if(indices==1L) "Portal fibroblast" else "Portal\nfibroblast-like"
    txt(lab,.10,yvals[role_name],size=11)
  }
  for(j in seq_along(indices)) {
    i <- indices[j]; start <- label_width+(j-1L)*(panel_width+gap)
    left <- start+.06; right <- start+panel_width-.09
    xp <- function(value) left+value/xmax*(right-left)
    title_x <- if(standalone) width/2 else start+panel_width/2+.025
    letter <- if(supplementary) LETTERS[i] else LETTERS[i+7L]
    txt(letter,if(standalone) .10 else start-.025,height-.15,size=16,face="bold")
    txt(cohorts[i],title_x,height-.16,size=14,face="bold",just="centre")
    ss <- stability[cohort==cohorts[i] & n_neighbors==k]
    txt(sprintf("%d cells per population",ss$n_per_role),title_x,height-.43,size=10,colour=grey,just="centre")
    txt(sprintf("200 resamples; k = %d",k),title_x,height-.65,size=9,colour=grey,just="centre")
    # Ratio 1 is the null expectation, not an upper bound.
    line(xp(1),axis_y,xp(1),max(yvals)+.19,colour="#BFC5CA",lwd=.6,lty="dashed")
    line(left,axis_y,right,axis_y,colour=ink,lwd=.65)
    for(value in seq(0,floor(xmax),by=1)) {
      line(xp(value),axis_y,xp(value),axis_y-.04,colour=ink,lwd=.65)
      txt(sprintf("%g",value),xp(value),axis_y-.145,size=10,just="centre")
    }
    for(role_name in shown) {
      yy <- yvals[role_name]
      z <- dat[cohort==cohorts[i] & n_neighbors==k & role==role_name]
      if(!nrow(z)) {
        txt("-",xp(xmax/2),yy,size=10,colour="#A0A5AA",just="centre"); next
      }
      values <- z[order(replicate),paga_connectivity]
      col <- palette[role_name]
      # Half-density is bounded to the observed range and omitted if variance is zero.
      # There is no horizontal jitter: the x coordinates are the actual frozen values.
      if(diff(range(values))>1e-12) {
        den <- density(values,n=512,from=min(values),to=max(values),bw="nrd0",adjust=.8)
        grid.polygon(x=c(xp(den$x),rev(xp(den$x))),
          y=c(yy+.025+den$y/max(den$y)*.155,rep(yy+.025,length(den$x))),
          default.units="native",gp=gpar(fill=adjustcolor(col,alpha.f=.50),col=col,lwd=.6))
      }
      set.seed(20260910L+1000L*i+100L*k+match(role_name,roles))
      grid.points(xp(values),yy-runif(length(values),.035,.13),default.units="native",pch=16,
        size=unit(.65,"mm"),gp=gpar(col=adjustcolor(col,alpha.f=.40)))
      stat <- summary[cohort==cohorts[i] & n_neighbors==k & role==role_name]
      line(xp(stat$q025),yy,xp(stat$q975),yy,colour=ink,lwd=.95)
      line(xp(stat$q025),yy-.021,xp(stat$q025),yy+.021,colour=ink,lwd=.65)
      line(xp(stat$q975),yy-.021,xp(stat$q975),yy+.021,colour=ink,lwd=.65)
      grid.points(xp(stat$median),yy,default.units="native",pch=21,size=unit(1.6,"mm"),
        gp=gpar(fill=ink,col="white",lwd=.35))
    }
    stability_x <- if(standalone) width/2 else start+panel_width/2
    txt("Portal > strongest epithelial",stability_x,.46,size=8.5,colour=grey,just="centre")
    txt(sprintf("%d / %d resamples",ss$portal_gt_epithelial,ss$n_resamples),
      stability_x,.265,size=11,face="bold",just="centre")
  }
  txt("Uncapped PAGA observed/expected ratio",width/2,.735,size=if(standalone) 10 else 11,just="centre")
  if(standalone) {
    txt("Black: median and 95% resampling interval",width/2,.105,size=8.5,colour=grey,just="centre")
  } else {
    txt("Dots: individual resamples; black: median and 95% resampling interval; dash: not included.",
      width/2,.105,size=8.5,colour=grey,just="centre")
  }
  popViewport()
}

save_chart <- function(stem,indices,k,width,height,standalone=FALSE,supplementary=FALSE) {
  quartz(type="pdf",file=file.path(out,paste0(stem,".pdf")),width=width,height=height,family=font,bg="white")
  draw_chart(indices,k,width,height,standalone,supplementary); invisible(dev.off())
  ragg::agg_png(file.path(out,paste0(stem,".png")),width=width,height=height,units="in",res=400,background="white")
  draw_chart(indices,k,width,height,standalone,supplementary); invisible(dev.off())
}
save_chart("Fig1_HIJ_balanced_resampling_FigureYa68_Arial",1:3,20,7.5,4.65)
for(i in 1:3) save_chart(paste0("Fig1",LETTERS[i+7L],"_balanced_resampling_",keys[i],"_Arial"),
  i,20,4.15,if(i==1L) 4.65 else 3.85,standalone=TRUE)
save_chart("Supplement_balanced_resampling_k10_FigureYa68_Arial",1:3,10,7.5,4.65,supplementary=TRUE)
audit <- dat[,.(n_rows=.N,n_resamples=uniqueN(replicate),n_per_role=unique(n_per_role),
  n_candidates=uniqueN(role),source_file=source_file[1],source_md5=source_md5[1]),by=.(cohort,n_neighbors)]
audit[,`:=`(analysis_recomputed=TRUE,values_changed=TRUE,roles_dropped=FALSE,
  metric="uncapped PAGA observed/expected ratio",recapped_summary_matches_upstream=TRUE,
  display_x_min=0,display_x_max=xmax)]
fwrite(audit,file.path(out,"plot_audit.csv"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
message("Completed H-J and k10 supplement: ",out)
