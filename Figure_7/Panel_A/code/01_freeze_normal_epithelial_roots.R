Sys.setenv(LANGUAGE="en"); options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table)})
set.seed(20260901)
OUTT <- "emt_invivo_bodyfig_20260901/tables"
dir.create(OUTT, recursive=TRUE, showWarnings=FALSE)
TARGET_N <- 150L
meta <- fread("results_v4/tables/v4_allcells_meta_slim.csv")
meta[, barcode := gsub('"','',barcode,fixed=TRUE)]
norm <- meta[is_malignant==FALSE & celltype %in% c("Hepatocyte","Epithelial (cholangiocyte)")]
plan <- list(
  list(cc="China-HCC",  ct="Hepatocyte",                   role="normal_hepatocyte_root"),
  list(cc="China-HCC",  ct="Epithelial (cholangiocyte)",  role="normal_cholangiocyte_root"),
  list(cc="USA-HCC",    ct="Hepatocyte",                   role="normal_hepatocyte_root"),
  list(cc="USA-HCC",    ct="Epithelial (cholangiocyte)",  role="normal_cholangiocyte_root"),
  list(cc="CCA-combined",ct="Epithelial (cholangiocyte)",  role="normal_cholangiocyte_root")
)
summ <- list()
for (p in plan) {
  pool <- norm[dataset_group==p$cc & celltype==p$ct]
  n <- nrow(pool)
  # Keep every available normal epithelial cell.  HCC has only a handful of
  # non-malignant cholangiocytes in the USA cohort (and a few dozen in China),
  # but these cells are explicitly part of the requested analysis and must not
  # be discarded merely for being below the 30-cell threshold.
  if (n == 0L) { cat("[skip]", p$cc, p$ct, "no cells available\n"); next }
  keep <- if (n > TARGET_N) pool[sample(.N, TARGET_N)] else pool
  d <- data.table(barcode=keep$barcode, cohort=p$cc, celltype=p$ct, frozen_role=p$role)
  stopifnot(nrow(d) <= TARGET_N, !anyDuplicated(d$barcode))
  tag <- if (p$cc=="China-HCC") "ChinaHCC" else if (p$cc=="USA-HCC") "USAHCC" else "CCAcombined"
  suffix <- if (grepl("chol", p$role)) "cholangiocyte" else "hepatocyte"
  fwrite(d, file.path(OUTT, sprintf("S1_%s_normal_%s_root_frozen_by_cell.csv", tag, suffix)))
  summ[[length(summ)+1]] <- data.table(cohort=p$cc, celltype=p$ct, role=p$role, pool_n=n, frozen_n=nrow(d))
}
fwrite(rbindlist(summ), file.path(OUTT,"S1_normal_root_freeze_summary.csv"))
print(rbindlist(summ))
