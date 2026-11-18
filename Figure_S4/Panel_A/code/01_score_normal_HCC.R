#!/usr/bin/env Rscript
# Score non-tumor liver epithelial cells of GSE149614 with the Figure 2A method.
#
# Same as 10abc_review_Cook2020_two_axis_EMT_v4.R: Cook 17-gene epithelial and
# Hallmark EMT (200 genes) AddModuleScore-equivalent scores, LogNormalize
# (scale 1e4) on unfiltered raw counts, the saved control genes, and the saved
# pooled malignant-cell medians as cutoffs. Genes absent from a dataset count as
# zero, as in the original union-merged object. The saved malignant scores are
# recomputed first to confirm the reproduction.
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(Matrix); library(data.table)})

argv <- commandArgs(trailingOnly = FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", argv[grepl("^--file=", argv)]))))
up <- normalizePath(file.path(root, "..", "..", "EMT_reproduce_final_20260910"))
out <- file.path(root, "output"); dir.create(out, showWarnings = FALSE)
source(file.path(root, "code", "score_functions.R"))

ref <- load_reference(up)

meta <- fread(file.path(up, "data/GSE149614_HBV_HCC/GSE149614_HCC.metadata.updated.txt.gz"))
cnt <- readRDS(file.path(up, "data/GSE149614_HBV_HCC/count1_sparse.rds"))
stopifnot(!anyDuplicated(colnames(cnt)), !anyDuplicated(rownames(cnt)))

# 1) reproduction check on the saved HCC-GSE149614 malignant cells
mal <- ref$saved[cohort == "China-HCC"]
check <- validate_against_saved(cnt[, mal$barcode], mal)
check[, dataset := "GSE149614"]
fwrite(check, file.path(out, "validation_GSE149614_malignant.csv"))
print(check)
stopifnot(check$max_abs_diff_E < 1e-8, check$max_abs_diff_M < 1e-8)

# 2) non-tumor liver epithelial cells (author annotation "Hepatocyte", site "Normal")
nt <- meta[site == "Normal" & celltype == "Hepatocyte"]
stopifnot(all(nt$Cell %in% colnames(cnt)))
s <- score_cells(lognorm(cnt[, nt$Cell]), ref)
s <- merge(s, nt[, .(barcode = Cell, sample, patient)], by = "barcode")
s[, group := "HCC-GSE149614"][, tissue := "Non-tumor liver"][, source_dataset := "GSE149614"]
s <- classify(s, ref)
fwrite(s, file.path(out, "normal_epithelium_scores_GSE149614.csv"))
print(s[, .N, by = .(patient)][order(patient)])
print(s[, .N, by = state][, pct := round(100 * N / sum(N), 1)][])
