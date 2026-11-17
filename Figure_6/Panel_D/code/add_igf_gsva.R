#!/usr/bin/env Rscript
# Add only IGF1R GSVA; retain the existing four pathway rows unchanged.
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(GSVA)
})
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)]))
work <- dirname(dirname(script)); root <- dirname(work)
inputs <- file.path(work, "audit", "inputs")
out <- file.path(work, "audit", "IGF")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
archive <- "/Users/sodajay/Desktop/Deutsch_Project/03_EMT/archive/EMT_reproduce_final_20260910"
obj_root <- file.path(archive, "emt_invivo_bodyfig_20260901")
frozen <- file.path(root, "Figure1_静态表达", "拼图_20261009", "inputs")

gmt <- readLines(file.path(inputs, "Reactome_R-HSA-2404192_20261024.gmt"), warn = FALSE)
stopifnot(length(gmt) == 1L)
fields <- strsplit(gmt, "\t", fixed = TRUE)[[1]]
stopifnot(fields[2] == "R-HSA-2404192")
igf <- fields[-c(1, 2)]
stopifnot(length(igf) == 54L, !anyDuplicated(igf), all(c("IGF1", "IGF2", "IGF1R") %in% igf))
fwrite(data.table(pathway = "IGF", source_pathway = fields[2], symbol = igf),
       file.path(out, "IGF_gene_membership.csv"))

cohorts <- c(ChinaHCC = "HCC-GSE149614", USAHCC = "HCC-GSE151530", CCAcombined = "CCA-GSE138709_GSE151530")
expected_md5 <- c(ChinaHCC = "06910d7e500ac6da7526787396966e96", USAHCC = "d621adf1a8b5bc081a0ad938cdbb9329", CCAcombined = "e19b2ca07f12f2985b37cdc256b48d50")
groups <- list(
  ChinaHCC = c("normal_cholangiocyte_root", "normal_hepatocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "activated_HSC_intermediate", "portal_fibroblast_root", "fibrotic_MFB_terminal"),
  USAHCC = c("normal_hepatocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal"),
  CCAcombined = c("normal_cholangiocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal")
)
labels <- c(normal_cholangiocyte_root = "Normal cholangiocyte", normal_hepatocyte_root = "Normal hepatocyte", malignant_epithelial_EMT_root = "Malignant epithelial", qHSC_root = "qHSC", activated_HSC_intermediate = "Activated HSC", portal_fibroblast_root = "Portal fibroblast", portal_fibroblast_like_root = "Portal fibroblast-like", fibrotic_MFB_terminal = "Fibrotic MFB")
original_rows <- c("TGF", "INSULIN-AKT", "HGF-AKT", "EGF-RAS")
order_five <- c("TGF", "INSULIN-AKT", "HGF-AKT", "IGF", "EGF-RAS")
input_info <- coverage_all <- role_all <- preserved <- list()
set.seed(1)
for (key in names(cohorts)) {
  path <- file.path(obj_root, paste0("S2_", key, "_joint_with_normal_roles_frozen.rds"))
  stopifnot(unname(tools::md5sum(path)) == expected_md5[[key]])
  message("[start] ", cohorts[[key]])
  obj <- readRDS(path)
  # Match the original four-pathway calculation exactly: this exclusion was
  # applied to USAHCC only; ChinaHCC's unused reference remains in scoring.
  if (key == "USAHCC") obj <- subset(obj, subset = frozen_role != "normal_cholangiocyte_root")
  DefaultAssay(obj) <- "RNA"
  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
  stopifnot(!anyDuplicated(rownames(cnt)), !anyDuplicated(colnames(cnt)))
  lib <- Matrix::colSums(cnt); lib[lib <= 0] <- 1
  keep <- Matrix::rowSums(cnt) > 0
  expr <- as.matrix(log1p(t(t(cnt[keep, , drop = FALSE]) / lib) * 1e4))
  matched <- intersect(igf, rownames(expr))
  variable <- rownames(expr)[matrixStats::rowVars(expr) > 0]
  stopifnot(length(matched) >= 3L)
  # Only the new IGF set is scored. Existing pathways are not recalculated.
  param <- GSVA::gsvaParam(expr, list(IGF = matched), kcdf = "Gaussian", minSize = 3L, maxSize = Inf)
  scores <- as.matrix(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
  stopifnot(identical(rownames(scores), "IGF"), identical(colnames(scores), colnames(cnt)), all(is.finite(scores)))
  role <- as.character(obj$frozen_role)
  stopifnot(setequal(unique(role), groups[[key]]))
  means <- vapply(groups[[key]], function(g) mean(scores[1, role == g]), numeric(1))
  z <- as.numeric(scale(means)); z[!is.finite(z)] <- 0
  names(z) <- groups[[key]]
  displayed <- if (key == "ChinaHCC") groups[[key]] != "normal_cholangiocyte_root" else rep(TRUE, length(means))
  role_table <- data.table(cohort = cohorts[[key]], role_id = groups[[key]], role_label = unname(labels[groups[[key]]]),
    n_cells = as.integer(table(factor(role, levels = groups[[key]]))), mean_gsva = as.numeric(means),
    within_cohort_role_zscore = z, displayed_in_figure6 = displayed)
  fwrite(role_table, file.path(out, paste0(key, "_IGF_role_scores.csv")))
  fwrite(data.table(cohort = cohorts[[key]], barcode = colnames(cnt), role_id = role, IGF = as.numeric(scores[1, ])),
    file.path(out, paste0(key, "_IGF_per_cell.csv")))
  role_all[[key]] <- role_table
  coverage_all[[key]] <- data.table(cohort = cohorts[[key]], symbol = igf,
    present_in_counts = igf %in% rownames(cnt), expressed = igf %in% rownames(expr), variable = igf %in% variable)
  input_info[[key]] <- data.table(cohort = cohorts[[key]], source = path, md5 = unname(tools::md5sum(path)), cells_scored = ncol(cnt),
    roles_scored = length(groups[[key]]), genes_in_signature = length(igf), genes_used = sum(igf %in% variable))
  for (kind in c("by_role", "zscore")) {
    original_path <- file.path(frozen, paste0("fig2_", key, "_signaling_GSVA_activity_", kind, ".csv"))
    original <- fread(original_path)
    expected_columns <- c("pathway", unname(labels[groups[[key]]]))
    stopifnot(identical(names(original), expected_columns), identical(original$pathway, original_rows))
    before_hash <- unname(tools::md5sum(original_path))
    new <- as.data.table(as.list(c(list(pathway = "IGF"), setNames(as.list(if (kind == "by_role") as.numeric(means) else z), expected_columns[-1]))))
    fwrite(new, file.path(out, paste0("fig6_", key, "_IGF_GSVA_activity_", kind, ".csv")))
    combined <- rbind(original, new)
    combined <- combined[match(order_five, pathway)]
    stopifnot(identical(combined[pathway %in% original_rows], original))
    combined_path <- file.path(out, paste0("fig6_", key, "_signaling_GSVA_activity_5pathways_", kind, ".csv"))
    fwrite(combined, combined_path)
    # Keep the original numeric tokens byte-for-byte in the combined CSV too.
    original_lines <- readLines(original_path, warn = FALSE)
    igf_lines <- readLines(file.path(out, paste0("fig6_", key, "_IGF_GSVA_activity_", kind, ".csv")), warn = FALSE)
    stopifnot(length(original_lines) == 5L, length(igf_lines) == 2L)
    writeLines(c(original_lines[1:4], igf_lines[2], original_lines[5]), combined_path, useBytes = TRUE)
    # Confirm serialization preserves the original values to CSV precision.
    back <- fread(combined_path)[pathway %in% original_rows]
    delta <- max(abs(as.matrix(back[, -1]) - as.matrix(original[, -1])))
    stopifnot(delta < 1e-14, unname(tools::md5sum(original_path)) == before_hash)
    preserved[[paste(key, kind)]] <- data.table(cohort = cohorts[[key]], kind = kind, original_file = original_path,
      original_md5 = before_hash, old_values = nrow(original) * (ncol(original) - 1L), old_rows_recomputed = FALSE,
      maximum_csv_difference = delta, original_file_unchanged = TRUE)
    file.copy(original_path, file.path(inputs, basename(original_path)), overwrite = FALSE)
  }
  message("[done] ", key, " ", ncol(cnt), " cells; ", sum(igf %in% variable), "/54 genes used")
  rm(obj, cnt, expr, scores, param); invisible(gc())
}
fwrite(rbindlist(role_all), file.path(out, "IGF_all_role_scores.csv"))
fwrite(rbindlist(coverage_all), file.path(out, "IGF_gene_coverage.csv"))
fwrite(rbindlist(input_info), file.path(out, "input_objects.csv"))
fwrite(rbindlist(preserved), file.path(out, "original_four_rows_preserved.csv"))
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
writeLines(c(
  "Added pathway: IGF, official Reactome R-HSA-2404192, 54 exact symbols from downloaded GMT.",
  "Original four rows are copied from frozen CSV files, not recomputed or renamed.",
  "Input: same frozen RNA count matrices as original 18_signaling_activity_gsva.R.",
  "USAHCC normal_cholangiocyte_root excluded before scoring; ChinaHCC retains all seven original roles.",
  "Normalization: log1p(count / complete cell library total * 10000); remove genes with zero total count.",
  "GSVA Gaussian kernel, minSize=3, maxSize=Inf, installed-version defaults, SerialParam.",
  "Role means are standardized using sample SD across the original roles within each cohort.",
  "ChinaHCC Normal cholangiocyte is hidden only at plotting, after calculating its IGF row z-score.",
  "Row order for drawing: TGF-beta, existing PI3K-AKT display row, HGF, IGF, EGF."
), file.path(out, "calculation_parameters.txt"))
message("[verified] IGF-only computation complete; all 136 original CSV values retained unchanged.")
