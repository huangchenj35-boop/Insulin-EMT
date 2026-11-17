#!/usr/bin/env Rscript

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

ROOT <- "emt_invivo_bodyfig_20260901"
OUTT <- file.path(ROOT, "tables")
dir.create(OUTT, recursive = TRUE, showWarnings = FALSE)

cohorts <- data.table(
  tag = c("ChinaHCC", "USAHCC", "CCAcombined"),
  cohort = c("China-HCC", "USA-HCC", "CCA-combined")
)
cohort_filter <- Sys.getenv("EMT_COHORT", "")
if (nzchar(cohort_filter)) {
  cohorts <- cohorts[tag == cohort_filter]
  if (nrow(cohorts) == 0L) stop("Unknown EMT_COHORT: ", cohort_filter)
}

for (i in seq_len(nrow(cohorts))) {
  tag <- cohorts$tag[i]
  cohort_name <- cohorts$cohort[i]
  object_file <- file.path(ROOT, sprintf("S2_%s_joint_with_normal_roles_frozen.rds", tag))
  input_file <- file.path(OUTT, sprintf("S2_%s_joint_with_normal_roles_input.csv", tag))
  output_file <- file.path(OUTT, sprintf("S3_%s_PAGA_Harmony_input.csv", tag))

  cells <- fread(input_file)
  cells[, barcode := gsub('"', '', barcode, fixed = TRUE)]
  if (tag == "USAHCC") {
    cells <- cells[frozen_role != "normal_cholangiocyte_root"]
  }
  stopifnot(
    all(c("barcode", "cohort", "frozen_role") %in% names(cells)),
    !anyDuplicated(cells$barcode),
    all(cells$cohort == cohort_name)
  )

  seu <- readRDS(object_file)
  stopifnot("harmony" %in% Reductions(seu))
  harmony <- Embeddings(seu, "harmony")
  stopifnot(ncol(harmony) == 30L, all(cells$barcode %in% rownames(harmony)))
  harmony <- harmony[cells$barcode, , drop = FALSE]
  colnames(harmony) <- sprintf("Harmony_%02d", seq_len(ncol(harmony)))

  out <- cbind(cells[, .(barcode, cohort, frozen_role)], as.data.table(harmony))
  stopifnot(nrow(out) == nrow(cells), identical(out$barcode, rownames(harmony)))
  fwrite(out, output_file)

  cat(sprintf("[done] %s: %d cells, %d Harmony dimensions\n", tag, nrow(out), ncol(harmony)))
  print(out[, .N, by = frozen_role])
  rm(seu, harmony, out, cells)
  gc(verbose = FALSE)
}
