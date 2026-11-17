#!/usr/bin/env Rscript
# Six canonical Reactome member sets, fixed before scoring. No figure output.
# All output is isolated from the earlier Hallmark comparison and V16.
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")
packages <- c("Seurat", "Matrix", "data.table", "GSVA", "BiocParallel", "matrixStats", "jsonlite")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages; no installation attempted: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table); library(GSVA)})
arg <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", arg[grepl("^--file=", arg)]), mustWork = TRUE)
stage <- dirname(dirname(script)); root <- dirname(stage)
work <- file.path(stage, "reactome_six")
tables <- file.path(work, "tables"); audit <- file.path(work, "audit")
for (path in c(tables, audit)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
previous <- file.path(root, "V16_Insulin与PI3K分开_20261026")
archive <- "/Users/sodajay/Desktop/Deutsch_Project/03_EMT/archive/EMT_reproduce_final_20260910/emt_invivo_bodyfig_20260901"
started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
md5 <- function(path) unname(tools::md5sum(path))
reuse_exports <- "--reuse-exported-scores" %in% commandArgs(trailingOnly = TRUE)
recovery <- if (reuse_exports) jsonlite::fromJSON(file.path(audit, "export_recovery_checkpoint.json")) else NULL
gmt_path <- file.path(stage, "inputs", "reactome_six", "official_six_pathways.gmt")
control_path <- file.path(previous, "inputs", "KEGG_hsa04151_PI3K_AKT.gmt")
old_igf_zip <- file.path(root, "V6_最终图_20261024", "audit", "inputs", "ReactomePathways_current_20261024.gmt.zip")
old_igf_directory <- file.path(root, "V6_最终图_20261024", "audit", "IGF")
definitions <- data.table(
  label = c("TGF-β", "PI3K-AKT", "Insulin", "HGF-MET", "EGFR", "IGF1R"),
  reactome_id = c("R-HSA-170834", "R-HSA-1257604", "R-HSA-74752", "R-HSA-6806834", "R-HSA-177929", "R-HSA-2404192"),
  V16_label = c("TGF-β", "PI3K-AKT", "Insulin", "HGF", "EGF", "IGF"))
row_order <- definitions$label; control_name <- "KEGG_PI3K_AKT_control"
stopifnot(file.exists(gmt_path), file.exists(control_path), file.exists(old_igf_zip), packageVersion("GSVA") == "2.6.6")
raw_lines <- readLines(gmt_path, warn = FALSE); raw_lines <- raw_lines[nzchar(trimws(raw_lines))]
stopifnot(length(raw_lines) == 6L)
parsed <- strsplit(raw_lines, "\t", fixed = TRUE)
extract_id <- function(fields) {
  stopifnot(length(fields) >= 5L)
  ids <- unique(unlist(regmatches(fields[1:2], gregexpr("R-HSA-[0-9]+", fields[1:2]))))
  stopifnot(length(ids) == 1L); ids
}
ids <- vapply(parsed, extract_id, character(1))
stopifnot(!anyDuplicated(ids), setequal(ids, definitions$reactome_id))
parsed <- parsed[match(definitions$reactome_id, ids)]
gene_sets <- setNames(lapply(parsed, function(fields) {
  genes <- fields[-c(1, 2)]
  stopifnot(all(nzchar(genes)), !anyDuplicated(genes), all(genes == trimws(genes)))
  genes
}), row_order)
definitions[, original_gmt_name := vapply(parsed, `[`, character(1), 1L)]
definitions[, original_gmt_description := vapply(parsed, `[`, character(1), 2L)]
definitions[, complete_member_count := lengths(gene_sets)]
stopifnot(identical(definitions$complete_member_count, c(94L, 284L, 82L, 79L, 53L, 54L)))
control_fields <- strsplit(readLines(control_path, warn = FALSE), "\t", fixed = TRUE)
stopifnot(length(control_fields) == 1L, grepl("hsa04151", control_fields[[1]][2], fixed = TRUE))
control_genes <- control_fields[[1]][-c(1, 2)]
stopifnot(length(control_genes) == 362L, !anyDuplicated(control_genes))
gene_sets[[control_name]] <- control_genes

# Compare the actual archived V16 IGF membership, not only its stable pathway ID.
old_igf_lines <- strsplit(readLines(unz(old_igf_zip, "ReactomePathways.gmt"), warn = FALSE), "\t", fixed = TRUE)
for (i in seq_len(nrow(definitions))) {
  source_match <- old_igf_lines[vapply(old_igf_lines, function(x) length(x) > 2L && x[2] == definitions$reactome_id[i], logical(1))]
  stopifnot(length(source_match) == 1L, identical(parsed[[i]], source_match[[1]]))
}
old_igf_match <- old_igf_lines[vapply(old_igf_lines, function(x) length(x) > 2L && x[2] == "R-HSA-2404192", logical(1))]
stopifnot(length(old_igf_match) == 1L)
old_igf_genes <- old_igf_match[[1]][-c(1, 2)]
igf_membership_same <- setequal(old_igf_genes, gene_sets[["IGF1R"]])
igf_union <- sort(union(old_igf_genes, gene_sets[["IGF1R"]]))
fwrite(data.table(gene_symbol = igf_union, in_V16_IGF = igf_union %in% old_igf_genes,
  in_current_IGF1R = igf_union %in% gene_sets[["IGF1R"]]), file.path(tables, "V16_IGF_vs_current_IGF1R_membership.csv"))
membership <- rbindlist(lapply(names(gene_sets), function(name) data.table(pathway = name,
  stable_id = if (name == control_name) "hsa04151" else definitions$reactome_id[match(name, definitions$label)],
  scoring_identifier = gene_sets[[name]], has_gene_symbol = !grepl("^hsa:", gene_sets[[name]]),
  identifier_policy = "Official identifier unchanged; no guessed aliases or overlap removal",
  use = if (name == control_name) "frozen_KEGG_reproduction_control" else "preselected_canonical_Reactome_comparison")))
fwrite(definitions, file.path(tables, "pathway_definitions.csv"))
fwrite(membership, file.path(tables, "complete_gene_set_membership.csv"))
cohorts <- c(ChinaHCC = "HCC-GSE149614", USAHCC = "HCC-GSE151530", CCAcombined = "CCA-GSE138709_GSE151530")
old_igf_percell_paths <- file.path(old_igf_directory, paste0(names(cohorts), "_IGF_per_cell.csv"))
protected_paths <- unique(c(gmt_path, control_path, old_igf_zip, old_igf_percell_paths,
  list.files(file.path(previous, "tables"), pattern = "\\.csv$", full.names = TRUE),
  list.files(file.path(stage, "tables"), pattern = "\\.csv$", full.names = TRUE)))
protected_before <- setNames(vapply(protected_paths, md5, character(1)), protected_paths)
stopifnot(!anyNA(protected_before))
if (reuse_exports) stopifnot(md5(gmt_path) == recovery$gmt_md5, md5(control_path) == recovery$control_gmt_md5)
expected_md5 <- c(ChinaHCC = "06910d7e500ac6da7526787396966e96", USAHCC = "d621adf1a8b5bc081a0ad938cdbb9329", CCAcombined = "e19b2ca07f12f2985b37cdc256b48d50")
expected_scored_cells <- c(ChinaHCC = 705L, USAHCC = 454L, CCAcombined = 369L)
groups <- list(
  ChinaHCC = c("normal_cholangiocyte_root", "normal_hepatocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "activated_HSC_intermediate", "portal_fibroblast_root", "fibrotic_MFB_terminal"),
  USAHCC = c("normal_hepatocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal"),
  CCAcombined = c("normal_cholangiocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal"))
labels <- c(normal_cholangiocyte_root = "Normal cholangiocyte", normal_hepatocyte_root = "Normal hepatocyte", malignant_epithelial_EMT_root = "Malignant epithelial", qHSC_root = "qHSC", activated_HSC_intermediate = "Activated HSC", portal_fibroblast_root = "Portal fibroblast", portal_fibroblast_like_root = "Portal fibroblast-like", fibrotic_MFB_terminal = "Fibrotic MFB")
safe_spearman <- function(x, y) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)))
  if (length(x) < 3L) return(list(rho = NA_real_, status = "fewer_than_3_roles"))
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(list(rho = NA_real_, status = "constant_role_scores"))
  list(rho = unname(stats::cor(x, y, method = "spearman")), status = "ok")
}
input_records <- coverage_records <- control_records <- igf_records <- role_records <- comparison_records <- correlation_records <- list()

for (key in names(cohorts)) {
  message("[start] ", key)
  object_path <- file.path(archive, paste0("S2_", key, "_joint_with_normal_roles_frozen.rds"))
  before_md5 <- md5(object_path); stopifnot(identical(before_md5, unname(expected_md5[[key]])))
  obj <- readRDS(object_path); original_cells <- ncol(obj)
  if (key == "USAHCC") obj <- obj[, as.character(obj$frozen_role) != "normal_cholangiocyte_root"]
  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts"); role <- as.character(obj$frozen_role)
  stopifnot(ncol(cnt) == expected_scored_cells[[key]], identical(colnames(cnt), colnames(obj)),
    !anyDuplicated(rownames(cnt)), !anyDuplicated(colnames(cnt)), !anyNA(role), setequal(unique(role), groups[[key]]))
  libraries <- Matrix::colSums(cnt); libraries[libraries <= 0] <- 1
  keep <- Matrix::rowSums(cnt) > 0
  expr <- as.matrix(log1p(t(t(cnt[keep, , drop = FALSE]) / libraries) * 1e4))
  variable <- rownames(expr)[matrixStats::rowVars(expr) > 0]
  used_sets <- lapply(gene_sets, intersect, y = rownames(expr))
  stopifnot(all(vapply(used_sets, length, integer(1)) >= 3L))
  for (name in names(gene_sets)) {
    genes <- gene_sets[[name]]
    coverage_records[[paste(key, name)]] <- data.table(cohort_key = key, cohort = cohorts[[key]], pathway = name,
      stable_id = if (name == control_name) "hsa04151" else definitions$reactome_id[match(name, definitions$label)],
      scoring_identifier = genes, has_gene_symbol = !grepl("^hsa:", genes),
      present_in_counts = genes %in% rownames(cnt), detected_in_any_cell = genes %in% rownames(cnt)[keep],
      variable = genes %in% variable, supplied_to_GSVA = genes %in% used_sets[[name]])
  }
  set.seed(1)
  param <- GSVA::gsvaParam(expr, used_sets, kcdf = "Gaussian", minSize = 3L, maxSize = Inf)
  capture.output(show(param), file = file.path(audit, paste0(key, "_GSVA_parameter_object.txt")))
  if (reuse_exports) {
    cache_path <- file.path(tables, paste0(key, "_six_pathway_per_cell_scores.csv"))
    stopifnot(md5(cache_path) == recovery$per_cell_md5[[key]], before_md5 == recovery$frozen_object_md5[[key]])
    cached <- fread(cache_path)
    stopifnot(identical(cached$barcode, colnames(cnt)), identical(cached$role_id, role))
    scores <- t(as.matrix(cached[, names(gene_sets), with = FALSE]))
    dimnames(scores) <- list(names(gene_sets), cached$barcode)
    message("[recover] ", key, ": using checksummed completed per-cell scores; GSVA is not rerun")
  } else {
    scores <- as.matrix(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
  }
  stopifnot(identical(rownames(scores), names(gene_sets)), identical(colnames(scores), colnames(cnt)), all(is.finite(scores)))
  old_percell <- fread(file.path(previous, "tables", paste0(key, "_two_pathway_per_cell_scores.csv")))
  stopifnot(nrow(old_percell) == ncol(scores), !anyDuplicated(old_percell$barcode), setequal(old_percell$barcode, colnames(scores)))
  old_percell <- old_percell[match(colnames(scores), old_percell$barcode)]
  stopifnot(identical(old_percell$barcode, colnames(scores)), identical(old_percell$role_id, role))
  control_error <- max(abs(scores[control_name, ] - old_percell[["PI3K-AKT"]]))
  stopifnot(is.finite(control_error), control_error <= 1e-10)
  old_igf <- fread(file.path(old_igf_directory, paste0(key, "_IGF_per_cell.csv")))
  stopifnot(nrow(old_igf) == ncol(scores), !anyDuplicated(old_igf$barcode), setequal(old_igf$barcode, colnames(scores)))
  old_igf <- old_igf[match(colnames(scores), old_igf$barcode)]
  stopifnot(identical(old_igf$barcode, colnames(scores)), identical(old_igf$role_id, role))
  igf_error <- max(abs(scores["IGF1R", ] - old_igf$IGF))
  same_effective_igf_members <- setequal(intersect(old_igf_genes, rownames(expr)), used_sets[["IGF1R"]])
  if (same_effective_igf_members) stopifnot(igf_error <= 1e-10)
  displayed <- !(key == "ChinaHCC" & groups[[key]] == "normal_cholangiocyte_root")
  per_cell <- data.table(cohort_key = key, cohort = cohorts[[key]], barcode = colnames(cnt), role_id = role,
    role_label = unname(labels[role]), displayed_in_figure6 = role %in% groups[[key]][displayed])
  for (name in names(gene_sets)) per_cell[[name]] <- as.numeric(scores[name, ])
  per_cell[["V16_KEGG_PI3K_AKT"]] <- old_percell[["PI3K-AKT"]]
  per_cell[["KEGG_control_difference"]] <- as.numeric(scores[control_name, ] - old_percell[["PI3K-AKT"]])
  per_cell[["V16_IGF"]] <- old_igf$IGF
  per_cell[["IGF1R_minus_V16_IGF"]] <- as.numeric(scores["IGF1R", ] - old_igf$IGF)
  fwrite(per_cell, file.path(tables, paste0(key, "_six_pathway_per_cell_scores.csv")))
  means_all <- vapply(groups[[key]], function(g) rowMeans(scores[, role == g, drop = FALSE]), numeric(length(gene_sets)))
  dimnames(means_all) <- list(names(gene_sets), unname(labels[groups[[key]]]))
  z_all <- t(scale(t(means_all))); z_all[!is.finite(z_all)] <- 0
  means <- means_all[row_order, , drop = FALSE]; zscores <- z_all[row_order, , drop = FALSE]
  fwrite(as.data.table(means, keep.rownames = "pathway"), file.path(tables, paste0(key, "_six_pathway_role_means.csv")))
  fwrite(as.data.table(zscores, keep.rownames = "pathway"), file.path(tables, paste0(key, "_six_pathway_role_zscores.csv")))
  old_tables <- list()
  for (kind in c("by_role", "zscore")) {
    old <- fread(file.path(previous, "tables", paste0("fig6_", key, "_signaling_GSVA_activity_6pathways_", kind, ".csv")))
    stopifnot(identical(names(old), c("pathway", unname(labels[groups[[key]]]))),
      identical(old$pathway, c("TGF-β", "PI3K-AKT", "Insulin", "HGF", "EGF", "IGF")))
    old_mat <- as.matrix(old[, -1]); rownames(old_mat) <- old$pathway; old_tables[[kind]] <- old_mat
    result <- if (kind == "by_role") means else zscores
    control_role_error <- max(abs((if (kind == "by_role") means_all else z_all)[control_name, ] - old_mat["PI3K-AKT", ]))
    stopifnot(control_role_error <= 1e-10)
    igf_role_error <- max(abs(result["IGF1R", ] - old_mat["IGF", ]))
    if (same_effective_igf_members) stopifnot(igf_role_error <= 1e-10)
    fwrite(as.data.table(result, keep.rownames = "pathway"),
      file.path(tables, paste0("fig6_", key, "_signaling_GSVA_activity_6pathways_", kind, ".csv")))
    fwrite(as.data.table(result[, displayed, drop = FALSE], keep.rownames = "pathway"),
      file.path(tables, paste0("fig6_", key, "_signaling_GSVA_activity_6pathways_", kind, "_displayed.csv")))
    control_records[[paste(key, kind)]] <- data.table(cohort_key = key, comparison = kind,
      values_compared = length(groups[[key]]), maximum_absolute_difference = control_role_error, tolerance = 1e-10, pass = TRUE)
    igf_records[[paste(key, kind)]] <- data.table(cohort_key = key, comparison = kind,
      same_complete_membership = igf_membership_same, same_effective_membership = same_effective_igf_members,
      values_compared = length(groups[[key]]), maximum_absolute_difference = igf_role_error,
      matches_within_1e_10 = igf_role_error <= 1e-10,
      explanation = if (same_effective_igf_members) "Same effective members and pipeline; only floating-point precision differences allowed" else "Official member set differs; inspect saved membership comparison and gene coverage")
  }
  control_records[[paste(key, "per_cell")]] <- data.table(cohort_key = key, comparison = "per_cell",
    values_compared = ncol(scores), maximum_absolute_difference = control_error, tolerance = 1e-10, pass = TRUE)
  igf_records[[paste(key, "per_cell")]] <- data.table(cohort_key = key, comparison = "per_cell",
    same_complete_membership = igf_membership_same, same_effective_membership = same_effective_igf_members,
    values_compared = ncol(scores), maximum_absolute_difference = igf_error, matches_within_1e_10 = igf_error <= 1e-10,
    explanation = if (same_effective_igf_members) "Same effective members and pipeline; only floating-point precision differences allowed" else "Official member set differs; inspect saved membership comparison and gene coverage")
  role_records[[key]] <- rbindlist(lapply(row_order, function(name) data.table(cohort_key = key, cohort = cohorts[[key]],
    pathway = name, stable_id = definitions$reactome_id[match(name, definitions$label)], role_id = groups[[key]],
    role_label = unname(labels[groups[[key]]]), n_cells = as.integer(table(factor(role, levels = groups[[key]]))),
    mean_gsva = as.numeric(means[name, ]), within_cohort_role_zscore = as.numeric(zscores[name, ]),
    displayed_in_figure6 = displayed, recomputed = TRUE)))
  for (i in seq_len(nrow(definitions))) {
    name <- definitions$label[i]; old_name <- definitions$V16_label[i]
    old_mean <- old_tables$by_role[old_name, ]; new_mean <- means[name, ]
    old_z <- old_tables$zscore[old_name, ]; new_z <- zscores[name, ]
    old_display_rank <- new_display_rank <- rep(NA_real_, length(displayed))
    old_display_rank[displayed] <- rank(-old_mean[displayed], ties.method = "average")
    new_display_rank[displayed] <- rank(-new_mean[displayed], ties.method = "average")
    comparison_records[[paste(key, name)]] <- data.table(cohort_key = key, cohort = cohorts[[key]],
      previous_pathway = old_name, new_pathway = name, reactome_id = definitions$reactome_id[i],
      role_id = groups[[key]], role_label = unname(labels[groups[[key]]]), displayed_in_figure6 = displayed,
      previous_mean_GSVA = as.numeric(old_mean), new_mean_GSVA = as.numeric(new_mean),
      previous_Z = as.numeric(old_z), new_Z = as.numeric(new_z), Z_difference = as.numeric(new_z - old_z),
      previous_Z_sign = sign(old_z), new_Z_sign = sign(new_z), opposite_nonzero_Z_sign = old_z * new_z < 0,
      any_Z_sign_change = sign(old_z) != sign(new_z), previous_rank_all_original_roles = rank(-old_mean, ties.method = "average"),
      new_rank_all_original_roles = rank(-new_mean, ties.method = "average"),
      previous_rank_displayed_roles = old_display_rank, new_rank_displayed_roles = new_display_rank,
      original_role_count = length(displayed), displayed_role_count = sum(displayed))
    for (scope in c("all_original_roles", "displayed_roles")) {
      selected <- if (scope == "all_original_roles") rep(TRUE, length(displayed)) else displayed
      correlation <- safe_spearman(old_mean[selected], new_mean[selected])
      correlation_records[[paste(key, name, scope)]] <- data.table(cohort_key = key, cohort = cohorts[[key]],
        previous_pathway = old_name, new_pathway = name, reactome_id = definitions$reactome_id[i], scope = scope,
        n_roles = sum(selected), spearman_rho_role_means = correlation$rho, status = correlation$status,
        p_values_computed = FALSE, interpretation = "Descriptive comparison of role-level score ordering only")
    }
  }
  after_md5 <- md5(object_path); stopifnot(identical(before_md5, after_md5))
  input_records[[key]] <- data.table(cohort_key = key, source = object_path, expected_md5 = expected_md5[[key]],
    input_md5_before = before_md5, input_md5_after = after_md5, input_unchanged = TRUE,
    original_cells = original_cells, cells_scored = ncol(cnt), excluded_cells = original_cells - ncol(cnt),
    roles_scored = length(groups[[key]]), roles_displayed = sum(displayed), genes_in_counts = nrow(cnt),
    nonzero_genes = nrow(expr), variable_genes = length(variable))
  message("[done] ", key, ": ", ncol(cnt), " cells; KEGG per-cell error ", format(control_error),
    "; IGF per-cell error ", format(igf_error), "; IGF effective members identical: ", same_effective_igf_members)
  rm(obj, cnt, expr, scores, param, per_cell); invisible(gc())
}
coverage <- rbindlist(coverage_records)
coverage_summary <- coverage[, .(complete_member_count = .N, identifiers_with_symbol = sum(has_gene_symbol),
  present_in_counts = sum(present_in_counts), detected_in_any_cell = sum(detected_in_any_cell),
  variable = sum(variable), supplied_to_GSVA = sum(supplied_to_GSVA)), by = .(cohort_key, cohort, pathway, stable_id)]
comparisons <- rbindlist(comparison_records); correlations <- rbindlist(correlation_records)
fwrite(coverage, file.path(tables, "gene_coverage.csv"))
fwrite(coverage_summary, file.path(tables, "gene_coverage_summary.csv"))
fwrite(coverage[supplied_to_GSVA == FALSE], file.path(tables, "unavailable_members_by_cohort.csv"))
fwrite(rbindlist(role_records), file.path(tables, "six_pathways_all_role_scores.csv"))
fwrite(comparisons, file.path(tables, "V16_vs_Reactome_role_comparison.csv"))
fwrite(correlations, file.path(tables, "V16_vs_Reactome_role_spearman_descriptive.csv"), na = "NA")
fwrite(comparisons[role_id == "malignant_epithelial_EMT_root"], file.path(tables, "V16_vs_Reactome_malignant_cell_summary.csv"))
fwrite(comparisons[, .(roles = .N, opposite_nonzero_Z_sign_count = sum(opposite_nonzero_Z_sign),
  any_Z_sign_change_count = sum(any_Z_sign_change)), by = .(cohort_key, previous_pathway, new_pathway, displayed_in_figure6)],
  file.path(tables, "V16_vs_Reactome_Z_sign_change_summary.csv"))
fwrite(rbindlist(control_records), file.path(audit, "KEGG_reproduction_checks.csv"))
fwrite(rbindlist(igf_records), file.path(audit, "IGF_reproduction_checks.csv"))
fwrite(rbindlist(input_records), file.path(audit, "input_objects_verified.csv"))
protected_after <- setNames(vapply(protected_paths, md5, character(1)), protected_paths)
stopifnot(identical(protected_before, protected_after))
fwrite(data.table(path = protected_paths, md5_before = unname(protected_before), md5_after = unname(protected_after), unchanged = TRUE),
  file.path(audit, "protected_input_hashes.csv"))
capture.output(sessionInfo(), file = file.path(audit, "R_sessionInfo.txt"))
report <- list(status = "completed", purpose = "Six preselected canonical Reactome pathway memberships scored uniformly; no figure output",
  started = started, completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), script = script, script_md5 = md5(script),
  recovered_completed_exports = reuse_exports,
  recovery_reason = if (reuse_exports) "All three GSVA calculations completed and controls passed; resumed existing checksummed scores after a data.table export-filter syntax correction, with no parameter or membership changes" else NULL,
  preselection = "The six stable IDs and labels were fixed by pathway definitions before these scores were inspected; no result-based alternative set selection",
  definitions = definitions, gmt_path = normalizePath(gmt_path), gmt_md5 = protected_before[[gmt_path]],
  membership_policy = "Complete official symbols exactly as supplied; only expression availability filtering; no alias guessing, arbitrary member deletion, or shared-member removal",
  KEGG_control = list(standard_name = "KEGG hsa04151", complete_member_count = 362L,
    original_companion_set = "KEGG Insulin hsa04910", new_companion_sets = row_order,
    maximum_per_cell_difference = max(rbindlist(control_records)[comparison == "per_cell"]$maximum_absolute_difference),
    maximum_role_mean_or_Z_difference = max(rbindlist(control_records)[comparison != "per_cell"]$maximum_absolute_difference), tolerance = 1e-10),
  IGF_comparison = list(stable_id = "R-HSA-2404192", old_label = "IGF", new_label = "IGF1R",
    previous_archive = normalizePath(old_igf_zip), previous_archive_md5 = protected_before[[old_igf_zip]],
    old_member_count = length(old_igf_genes), new_member_count = length(gene_sets[["IGF1R"]]), complete_membership_identical = igf_membership_same,
    added_members = setdiff(gene_sets[["IGF1R"]], old_igf_genes), removed_members = setdiff(old_igf_genes, gene_sets[["IGF1R"]]),
    checks = rbindlist(igf_records), note = "Identical effective membership must reproduce at <=1e-10; membership changes are explicitly recorded rather than attributed to biology"),
  normalization = "log1p(count / complete original per-cell library total * 10000); zero-total genes removed; complete expression background retained",
  gsva = list(kcdf = "Gaussian", minSize = 3L, maxSize = "Inf", seed = 1L, BPPARAM = "BiocParallel::SerialParam",
    other_parameters = "Same GSVA 2.6.6 defaults as V16; parameter object captured per cohort", package_version = as.character(packageVersion("GSVA"))),
  cell_selection = "Same three frozen objects and frozen_role; USAHCC normal cholangiocytes excluded before scoring",
  aggregation = "Per-cell GSVA mean within each original role; sample-SD row Z across all original roles within each cohort",
  display = "ChinaHCC Normal cholangiocyte hidden after row standardization only; full role retained in full tables",
  comparison = "Old V16 versus current per-pathway role means, within-row Z, ranks and descriptive role-level Spearman; no P values and no cell-level pseudo-replication inference",
  row_order = row_order, coverage_summary = coverage_summary, input_objects = rbindlist(input_records),
  interpretation = "RNA enrichment of official pathway members; receptor-specific networks and the downstream PI3K network have different biological levels; no phosphorylation or ligand-specific attribution",
  checks = list(all_frozen_MD5_match = TRUE, V16_Hallmark_and_original_inputs_unchanged = TRUE,
    all_six_official_memberships_retained = TRUE, all_six_recomputed = TRUE, KEGG_control_pass = TRUE,
    same_effective_IGF_members_reproduce_if_applicable = TRUE, all_scores_finite = TRUE, no_figures_written = TRUE))
jsonlite::write_json(report, file.path(audit, "calculation_notes.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16, na = "null")
message("[verified] Six Reactome pathways recomputed; KEGG and IGF controls checked; V16, Hallmark and original files unchanged.")
