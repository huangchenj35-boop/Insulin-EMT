#!/usr/bin/env Rscript
# Describe GSE149614 non-tumor vs malignant epithelial E/M states: per-patient pairs,
# hepatocyte/cholangiocyte marker split, and Hallmark genes driving the M-score difference.
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(Matrix); library(data.table)})
argv <- commandArgs(trailingOnly = FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", argv[grepl("^--file=", argv)]))))
up <- normalizePath(file.path(root, "..", "..", "EMT_reproduce_final_20260910"))
out <- file.path(root, "output")
source(file.path(root, "code", "score_functions.R"))
ref <- load_reference(up)

nt <- fread(file.path(out, "normal_epithelium_scores_GSE149614.csv"))
mal <- classify(copy(ref$saved[cohort == "China-HCC"]), ref)
mal[, patient := sub("T$", "", patient)]
states <- c("Epithelial high", "E/M high", "Mesenchymal high", "Low both")

pair <- rbind(nt[, .(patient, tissue = "Non-tumor", state)], mal[, .(patient, tissue = "Malignant", state)])
pp <- dcast(pair[, .N, by = .(patient, tissue, state)], patient + tissue ~ state, value.var = "N", fill = 0)
pp[, cells := rowSums(.SD), .SDcols = intersect(states, names(pp))]
pp[, M_high_pct := round(100 * (`E/M high` + `Mesenchymal high`) / cells, 1)]
fwrite(pp, file.path(out, "HCC_GSE149614_state_by_patient.csv"))
print(pp[, .(patient, tissue, cells, M_high_pct)])

# Paired patients (both tissues): M-high % non-tumor vs malignant
w <- dcast(pp, patient ~ tissue, value.var = "M_high_pct")[!is.na(`Non-tumor`) & !is.na(Malignant)]
print(w); cat("non-tumor > malignant in", sum(w$`Non-tumor` > w$Malignant), "of", nrow(w), "patients\n")
print(wilcox.test(w$`Non-tumor`, w$Malignant, paired = TRUE, exact = TRUE))

# Gene contributions to the mean M-score difference (non-tumor minus malignant)
cnt <- readRDS(file.path(up, "data/GSE149614_HBV_HCC/count1_sparse.rds"))
x_nt <- lognorm(cnt[, nt$barcode]); x_m <- lognorm(cnt[, mal$barcode])
f <- intersect(ref$emt_features, rownames(x_nt))
contrib <- data.table(gene = f,
                      nontumor_mean = Matrix::rowMeans(x_nt[f, ]),
                      malignant_mean = Matrix::rowMeans(x_m[f, ]))
contrib[, contribution := (nontumor_mean - malignant_mean) / length(ref$emt_features)]
ctl_diff <- (mean(Matrix::colSums(x_nt[intersect(ref$emt_controls, rownames(x_nt)), ])) -
             mean(Matrix::colSums(x_m[intersect(ref$emt_controls, rownames(x_m)), ]))) / length(ref$emt_controls)
setorder(contrib, -contribution)
fwrite(contrib, file.path(out, "HCC_GSE149614_Hallmark_gene_contributions.csv"))
cat(sprintf("mean M: non-tumor %.4f, malignant %.4f; feature part %.4f, control part %.4f\n",
            mean(nt$hallmark_EMT_score), mean(mal$hallmark_EMT_score), sum(contrib$contribution), ctl_diff))
print(head(contrib, 12))

# Marker split of non-tumor epithelium
mk <- function(x, g) Matrix::colMeans(x[intersect(g, rownames(x)), , drop = FALSE])
nt[, chol_marker := mk(x_nt, c("KRT19", "KRT7", "EPCAM", "SOX9"))[barcode]]
nt[, hep_marker := mk(x_nt, c("ALB", "APOA1", "TF", "HP", "TTR"))[barcode]]
nt[, subtype := fifelse(chol_marker > hep_marker, "Cholangiocyte-like", "Hepatocyte-like")]
print(nt[, .N, by = .(subtype, state)][order(subtype, state)])
fwrite(nt, file.path(out, "normal_epithelium_scores_GSE149614.csv"))
