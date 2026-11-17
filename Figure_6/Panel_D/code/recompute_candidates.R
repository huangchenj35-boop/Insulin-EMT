#!/usr/bin/env Rscript
# Score all prespecified candidates; selection is read from a frozen input file.
# No result-based selection, member tuning, P values, or figure generation.
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")
packages <- c("Seurat", "Matrix", "data.table", "GSVA", "BiocParallel", "matrixStats", "jsonlite", "digest")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages; no installation attempted: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table); library(GSVA)})
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)]), mustWork = TRUE)
work <- dirname(dirname(script)); root <- dirname(work)
inputs <- file.path(work, "inputs"); tables <- file.path(work, "tables"); audit <- file.path(work, "audit")
for (path in c(tables, audit)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
previous <- file.path(root, "V16_Insulin与PI3K分开_20261026")
archive <- "/Users/sodajay/Desktop/Deutsch_Project/03_EMT/archive/EMT_reproduce_final_20260910/emt_invivo_bodyfig_20260901"
md5 <- function(path) unname(tools::md5sum(path))
sha256 <- function(path) digest::digest(file = path, algo = "sha256")
started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
spec_path <- file.path(inputs, "candidate_specifications.json")
manifest_path <- file.path(inputs, "official_source_manifest.json")
selection_path <- file.path(audit, "selection_before_scores.json")
stopifnot(file.exists(spec_path), file.exists(manifest_path), file.exists(selection_path), packageVersion("GSVA") == "2.6.6")
definitions <- as.data.table(jsonlite::fromJSON(spec_path))
selected <- as.data.table(jsonlite::fromJSON(selection_path)$selected)
stopifnot(all(c("standard_name", "topic", "source", "expected_members") %in% names(definitions)),
  all(c("label", "standard_name", "V16_label") %in% names(selected)), nrow(selected) == 6L,
  !anyDuplicated(definitions$standard_name), !anyDuplicated(selected$label), !anyDuplicated(selected$standard_name),
  all(selected$standard_name %in% definitions$standard_name))
definitions[, expected_members := as.integer(expected_members)]
selected_topics <- definitions$topic[match(selected$standard_name, definitions$standard_name)]
stopifnot(!anyDuplicated(selected_topics), setequal(selected_topics, unique(definitions$topic)))
definitions[, V16_label := selected$V16_label[match(topic, selected_topics)]]
definitions[, selected_label := selected$label[match(standard_name, selected$standard_name)]]
definitions[, selected_for_figure := standard_name %in% selected$standard_name]
stopifnot(!anyNA(definitions$V16_label), all(definitions$V16_label %in% c("TGF-β", "PI3K-AKT", "Insulin", "HGF", "EGF", "IGF")))
row_order <- selected$label; candidate_order <- definitions$standard_name
control_name <- "KEGG_PI3K_AKT_control"; stopifnot(!control_name %in% candidate_order)
manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
gene_sets <- list(); source_files <- c(spec_path, manifest_path, selection_path)
for (i in seq_len(nrow(definitions))) {
  standard <- definitions$standard_name[i]
  gmt_path <- file.path(inputs, paste0(standard, ".gmt")); json_path <- file.path(inputs, paste0(standard, ".json"))
  stopifnot(file.exists(gmt_path), file.exists(json_path), !is.null(manifest[[standard]]),
    sha256(gmt_path) == manifest[[standard]]$gmt$sha256, sha256(json_path) == manifest[[standard]]$json$sha256)
  lines <- readLines(gmt_path, warn = FALSE); lines <- lines[nzchar(trimws(lines))]
  stopifnot(length(lines) == 1L)
  fields <- strsplit(lines, "\t", fixed = TRUE)[[1]]; genes <- fields[-c(1, 2)]
  meta <- jsonlite::fromJSON(json_path)[[standard]]
  stopifnot(fields[1] == standard, length(genes) == definitions$expected_members[i], !anyDuplicated(genes),
    all(nzchar(genes)), all(genes == trimws(genes)), identical(genes, meta$geneSymbols),
    identical(genes, unlist(manifest[[standard]]$members, use.names = FALSE)),
    meta$systematicName == manifest[[standard]]$systematic_name)
  if (!is.null(manifest[[standard]]$collection)) stopifnot(meta$collection == manifest[[standard]]$collection)
  gene_sets[[standard]] <- genes
  definitions[i, `:=`(systematic_name = meta$systematicName, collection = meta$collection,
    official_url = meta$msigdbURL, complete_member_count = length(genes), GMT_SHA256 = sha256(gmt_path), JSON_SHA256 = sha256(json_path))]
  source_files <- c(source_files, gmt_path, json_path)
}
control_path <- file.path(previous, "inputs", "KEGG_hsa04151_PI3K_AKT.gmt")
control_fields <- strsplit(readLines(control_path, warn = FALSE), "\t", fixed = TRUE)
stopifnot(length(control_fields) == 1L, grepl("hsa04151", control_fields[[1]][2], fixed = TRUE))
control_genes <- control_fields[[1]][-c(1, 2)]
stopifnot(length(control_genes) == 362L, !anyDuplicated(control_genes))
gene_sets[[control_name]] <- control_genes
protected_paths <- unique(c(source_files, control_path, list.files(file.path(previous, "tables"), pattern = "\\.csv$", full.names = TRUE)))
protected_before <- setNames(vapply(protected_paths, md5, character(1)), protected_paths)
stopifnot(!anyNA(protected_before)); selection_sha256 <- sha256(selection_path)
membership <- rbindlist(lapply(names(gene_sets), function(name) data.table(standard_name = name,
  scoring_identifier = gene_sets[[name]], has_gene_symbol = !grepl("^hsa:", gene_sets[[name]]),
  member_order = seq_along(gene_sets[[name]]),
  purpose = if (name == control_name) "frozen_KEGG_reproduction_control" else "prespecified_candidate")))
selected_definitions <- definitions[match(selected$standard_name, standard_name)]
selected_definitions[, label := selected$label]
fwrite(definitions, file.path(tables, "candidate_definitions.csv"))
fwrite(selected_definitions, file.path(tables, "selected_pathway_definitions.csv"))
fwrite(membership, file.path(tables, "complete_gene_set_membership.csv"))
cohorts <- c(ChinaHCC = "HCC-GSE149614", USAHCC = "HCC-GSE151530", CCAcombined = "CCA-GSE138709_GSE151530")
expected_md5 <- c(ChinaHCC = "06910d7e500ac6da7526787396966e96", USAHCC = "d621adf1a8b5bc081a0ad938cdbb9329", CCAcombined = "e19b2ca07f12f2985b37cdc256b48d50")
expected_cells <- c(ChinaHCC = 705L, USAHCC = 454L, CCAcombined = 369L)
groups <- list(
  ChinaHCC = c("normal_cholangiocyte_root", "normal_hepatocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "activated_HSC_intermediate", "portal_fibroblast_root", "fibrotic_MFB_terminal"),
  USAHCC = c("normal_hepatocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal"),
  CCAcombined = c("normal_cholangiocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal"))
labels <- c(normal_cholangiocyte_root = "Normal cholangiocyte", normal_hepatocyte_root = "Normal hepatocyte", malignant_epithelial_EMT_root = "Malignant epithelial", qHSC_root = "qHSC", activated_HSC_intermediate = "Activated HSC", portal_fibroblast_root = "Portal fibroblast", portal_fibroblast_like_root = "Portal fibroblast-like", fibrotic_MFB_terminal = "Fibrotic MFB")
safe_spearman <- function(x, y) {
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  unname(stats::cor(x, y, method = "spearman"))
}
coverage_records <- input_records <- control_records <- role_records <- comparison_records <- correlation_records <- equal_stroma_records <- pairwise_stroma_records <- list()
for (key in names(cohorts)) {
  message("[start] ", key)
  object_path <- file.path(archive, paste0("S2_", key, "_joint_with_normal_roles_frozen.rds"))
  before_md5 <- md5(object_path); stopifnot(identical(before_md5, unname(expected_md5[[key]])))
  obj <- readRDS(object_path); original_cells <- ncol(obj)
  if (key == "USAHCC") obj <- obj[, as.character(obj$frozen_role) != "normal_cholangiocyte_root"]
  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts"); role <- as.character(obj$frozen_role)
  stopifnot(ncol(cnt) == expected_cells[[key]], identical(colnames(cnt), colnames(obj)),
    !anyDuplicated(rownames(cnt)), !anyDuplicated(colnames(cnt)), !anyNA(role), setequal(unique(role), groups[[key]]))
  libraries <- Matrix::colSums(cnt); libraries[libraries <= 0] <- 1
  keep <- Matrix::rowSums(cnt) > 0
  expr <- as.matrix(log1p(t(t(cnt[keep, , drop = FALSE]) / libraries) * 1e4))
  variable <- rownames(expr)[matrixStats::rowVars(expr) > 0]
  used_sets <- lapply(gene_sets, intersect, y = rownames(expr))
  stopifnot(all(lengths(used_sets) >= 3L))
  for (name in names(gene_sets)) {
    genes <- gene_sets[[name]]
    coverage_records[[paste(key, name)]] <- data.table(cohort_key = key, cohort = cohorts[[key]], standard_name = name,
      scoring_identifier = genes, has_gene_symbol = !grepl("^hsa:", genes), present_in_counts = genes %in% rownames(cnt),
      detected_in_any_cell = genes %in% rownames(cnt)[keep], variable = genes %in% variable, supplied_to_GSVA = genes %in% used_sets[[name]])
  }
  set.seed(1)
  param <- GSVA::gsvaParam(expr, used_sets, kcdf = "Gaussian", minSize = 3L, maxSize = Inf)
  capture.output(show(param), file = file.path(audit, paste0(key, "_GSVA_parameter_object.txt")))
  scores <- as.matrix(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
  stopifnot(identical(rownames(scores), names(gene_sets)), identical(colnames(scores), colnames(cnt)), all(is.finite(scores)))
  old_percell <- fread(file.path(previous, "tables", paste0(key, "_two_pathway_per_cell_scores.csv")))
  stopifnot(nrow(old_percell) == ncol(scores), !anyDuplicated(old_percell$barcode), setequal(old_percell$barcode, colnames(scores)))
  old_percell <- old_percell[match(colnames(scores), old_percell$barcode)]
  stopifnot(identical(old_percell$barcode, colnames(scores)), identical(old_percell$role_id, role))
  control_error <- max(abs(scores[control_name, ] - old_percell[["PI3K-AKT"]])); stopifnot(control_error <= 1e-10)
  displayed <- !(key == "ChinaHCC" & groups[[key]] == "normal_cholangiocyte_root")
  per_cell <- data.table(cohort_key = key, cohort = cohorts[[key]], barcode = colnames(cnt), role_id = role,
    role_label = unname(labels[role]), displayed_in_figure6 = role %in% groups[[key]][displayed])
  for (name in names(gene_sets)) per_cell[[name]] <- as.numeric(scores[name, ])
  per_cell[["V16_KEGG_PI3K_AKT"]] <- old_percell[["PI3K-AKT"]]
  per_cell[["KEGG_control_difference"]] <- as.numeric(scores[control_name, ] - old_percell[["PI3K-AKT"]])
  fwrite(per_cell, file.path(tables, paste0(key, "_all_candidate_per_cell_scores.csv")))
  means_all <- vapply(groups[[key]], function(g) rowMeans(scores[, role == g, drop = FALSE]), numeric(length(gene_sets)))
  dimnames(means_all) <- list(names(gene_sets), unname(labels[groups[[key]]]))
  z_all <- t(scale(t(means_all))); z_all[!is.finite(z_all)] <- 0
  means <- means_all[candidate_order, , drop = FALSE]; zscores <- z_all[candidate_order, , drop = FALSE]
  fwrite(as.data.table(means, keep.rownames = "standard_name"), file.path(tables, paste0(key, "_all_candidate_role_means.csv")))
  fwrite(as.data.table(zscores, keep.rownames = "standard_name"), file.path(tables, paste0(key, "_all_candidate_role_zscores.csv")))
  old_tables <- list()
  for (kind in c("by_role", "zscore")) {
    old <- fread(file.path(previous, "tables", paste0("fig6_", key, "_signaling_GSVA_activity_6pathways_", kind, ".csv")))
    stopifnot(identical(names(old), c("pathway", unname(labels[groups[[key]]]))),
      identical(old$pathway, c("TGF-β", "PI3K-AKT", "Insulin", "HGF", "EGF", "IGF")))
    old_mat <- as.matrix(old[, -1]); rownames(old_mat) <- old$pathway; old_tables[[kind]] <- old_mat
    result_all <- if (kind == "by_role") means_all else z_all
    error <- max(abs(result_all[control_name, ] - old_mat["PI3K-AKT", ])); stopifnot(error <= 1e-10)
    result <- result_all[selected$standard_name, , drop = FALSE]; rownames(result) <- selected$label
    fwrite(as.data.table(result, keep.rownames = "pathway"),
      file.path(tables, paste0("fig6_", key, "_signaling_GSVA_activity_6pathways_", kind, ".csv")))
    fwrite(as.data.table(result[, displayed, drop = FALSE], keep.rownames = "pathway"),
      file.path(tables, paste0("fig6_", key, "_signaling_GSVA_activity_6pathways_", kind, "_displayed.csv")))
    control_records[[paste(key, kind)]] <- data.table(cohort_key = key, comparison = kind,
      values_compared = length(groups[[key]]), maximum_absolute_difference = error, tolerance = 1e-10, pass = TRUE)
  }
  control_records[[paste(key, "per_cell")]] <- data.table(cohort_key = key, comparison = "per_cell",
    values_compared = ncol(scores), maximum_absolute_difference = control_error, tolerance = 1e-10, pass = TRUE)
  epithelial_indices <- which(grepl("^normal_", groups[[key]]) | groups[[key]] == "malignant_epithelial_EMT_root")
  stromal_indices <- setdiff(seq_along(groups[[key]]), epithelial_indices)
  stopifnot(length(stromal_indices) == if (key == "ChinaHCC") 4L else 3L)
  for (i in seq_len(nrow(definitions))) {
    name <- definitions$standard_name[i]; old_name <- definitions$V16_label[i]
    display_name <- if (definitions$selected_for_figure[i]) definitions$selected_label[i] else name
    identity <- data.table(standard_name = name, topic = definitions$topic[i], source = definitions$source[i],
      selected_for_figure = definitions$selected_for_figure[i], selected_label = definitions$selected_label[i])
    old_mean <- old_tables$by_role[old_name, ]; new_mean <- means[name, ]; old_z <- old_tables$zscore[old_name, ]; new_z <- zscores[name, ]
    role_records[[paste(key, name)]] <- cbind(identity, data.table(cohort_key = key, cohort = cohorts[[key]],
      pathway = display_name, role_id = groups[[key]], role_label = unname(labels[groups[[key]]]),
      n_cells = as.integer(table(factor(role, levels = groups[[key]]))), mean_gsva = as.numeric(new_mean),
      within_cohort_role_zscore = as.numeric(new_z), displayed_in_figure6 = displayed))
    old_display_rank <- new_display_rank <- rep(NA_real_, length(displayed))
    old_display_rank[displayed] <- rank(-old_mean[displayed], ties.method = "average")
    new_display_rank[displayed] <- rank(-new_mean[displayed], ties.method = "average")
    comparison_records[[paste(key, name)]] <- cbind(identity, data.table(cohort_key = key, cohort = cohorts[[key]],
      previous_pathway = old_name, new_pathway = display_name, role_id = groups[[key]], role_label = unname(labels[groups[[key]]]),
      displayed_in_figure6 = displayed, previous_mean_GSVA = as.numeric(old_mean), new_mean_GSVA = as.numeric(new_mean),
      previous_Z = as.numeric(old_z), new_Z = as.numeric(new_z), Z_difference = as.numeric(new_z - old_z),
      opposite_nonzero_Z_sign = old_z * new_z < 0, any_Z_sign_change = sign(old_z) != sign(new_z),
      previous_rank_all_original_roles = rank(-old_mean, ties.method = "average"), new_rank_all_original_roles = rank(-new_mean, ties.method = "average"),
      previous_rank_displayed_roles = old_display_rank, new_rank_displayed_roles = new_display_rank))
    for (scope in c("all_original_roles", "displayed_roles")) {
      use <- if (scope == "all_original_roles") rep(TRUE, length(displayed)) else displayed
      correlation_records[[paste(key, name, scope)]] <- cbind(identity, data.table(cohort_key = key, previous_pathway = old_name,
        new_pathway = display_name, scope = scope, n_roles = sum(use),
        spearman_rho_role_means = safe_spearman(old_mean[use], new_mean[use]), p_values_computed = FALSE))
    }
    for (e in epithelial_indices) {
      shared <- cbind(identity, data.table(cohort_key = key, cohort = cohorts[[key]], previous_pathway = old_name,
        new_pathway = display_name, epithelial_role_id = groups[[key]][e], epithelial_role_label = names(new_mean)[e],
        epithelial_displayed_in_figure6 = displayed[e], epithelial_type = if (groups[[key]][e] == "malignant_epithelial_EMT_root") "malignant" else "normal"))
      equal_stroma_records[[paste(key, name, e)]] <- cbind(shared, data.table(stromal_role_count = length(stromal_indices),
        stromal_roles = paste(groups[[key]][stromal_indices], collapse = ";"),
        previous_epithelial_mean = old_mean[e], new_epithelial_mean = new_mean[e],
        previous_equal_weight_stromal_role_mean = mean(old_mean[stromal_indices]), new_equal_weight_stromal_role_mean = mean(new_mean[stromal_indices]),
        previous_epithelial_minus_equal_weight_stroma = old_mean[e] - mean(old_mean[stromal_indices]),
        new_epithelial_minus_equal_weight_stroma = new_mean[e] - mean(new_mean[stromal_indices]),
        previous_epithelial_Z = old_z[e], new_epithelial_Z = new_z[e],
        previous_equal_weight_stromal_Z = mean(old_z[stromal_indices]), new_equal_weight_stromal_Z = mean(new_z[stromal_indices]),
        previous_Z_difference = old_z[e] - mean(old_z[stromal_indices]), new_Z_difference = new_z[e] - mean(new_z[stromal_indices]),
        weighting = "One equal weight per stromal role; not weighted by cell count", p_values_computed = FALSE))
      for (s in stromal_indices) pairwise_stroma_records[[paste(key, name, e, s)]] <- cbind(shared, data.table(
        stromal_role_id = groups[[key]][s], stromal_role_label = names(new_mean)[s],
        previous_epithelial_mean = old_mean[e], new_epithelial_mean = new_mean[e], previous_stromal_mean = old_mean[s], new_stromal_mean = new_mean[s],
        previous_epithelial_minus_stromal_role = old_mean[e] - old_mean[s], new_epithelial_minus_stromal_role = new_mean[e] - new_mean[s],
        previous_epithelial_Z = old_z[e], new_epithelial_Z = new_z[e], previous_stromal_Z = old_z[s], new_stromal_Z = new_z[s],
        previous_Z_difference = old_z[e] - old_z[s], new_Z_difference = new_z[e] - new_z[s], p_values_computed = FALSE))
    }
  }
  after_md5 <- md5(object_path); stopifnot(identical(before_md5, after_md5))
  input_records[[key]] <- data.table(cohort_key = key, source = object_path, expected_md5 = expected_md5[[key]],
    input_md5_before = before_md5, input_md5_after = after_md5, input_unchanged = TRUE, original_cells = original_cells,
    cells_scored = ncol(cnt), excluded_cells = original_cells - ncol(cnt), roles_scored = length(groups[[key]]),
    roles_displayed = sum(displayed), genes_in_counts = nrow(cnt), nonzero_genes = nrow(expr), variable_genes = length(variable))
  message("[done] ", key, ": ", ncol(cnt), " cells, ", nrow(definitions), " candidates; KEGG maximum per-cell difference ", format(control_error))
  rm(obj, cnt, expr, scores, param, per_cell); invisible(gc())
}
coverage <- rbindlist(coverage_records)
coverage_summary <- coverage[, .(complete_member_count = .N, identifiers_with_symbol = sum(has_gene_symbol),
  present_in_counts = sum(present_in_counts), detected_in_any_cell = sum(detected_in_any_cell),
  variable = sum(variable), supplied_to_GSVA = sum(supplied_to_GSVA)), by = .(cohort_key, cohort, standard_name)]
all_roles <- rbindlist(role_records); comparisons <- rbindlist(comparison_records); correlations <- rbindlist(correlation_records)
equal_stroma <- rbindlist(equal_stroma_records); pairwise_stroma <- rbindlist(pairwise_stroma_records)
fwrite(coverage, file.path(tables, "gene_coverage.csv"))
fwrite(coverage_summary, file.path(tables, "gene_coverage_summary.csv"))
fwrite(coverage[supplied_to_GSVA == FALSE], file.path(tables, "unavailable_members_by_cohort.csv"))
fwrite(all_roles, file.path(tables, "all_candidates_all_role_scores.csv"))
fwrite(all_roles[selected_for_figure == TRUE], file.path(tables, "six_pathways_all_role_scores.csv"))
fwrite(comparisons, file.path(tables, "V16_vs_all_candidates_role_comparison.csv"))
fwrite(comparisons[selected_for_figure == TRUE], file.path(tables, "V16_vs_selected_role_comparison.csv"))
fwrite(comparisons[role_id == "malignant_epithelial_EMT_root"], file.path(tables, "V16_vs_all_candidates_malignant_cell_summary.csv"))
fwrite(comparisons[role_id == "malignant_epithelial_EMT_root" & selected_for_figure == TRUE], file.path(tables, "V16_vs_selected_malignant_cell_summary.csv"))
fwrite(correlations, file.path(tables, "V16_vs_all_candidates_role_spearman_descriptive.csv"), na = "NA")
fwrite(correlations[selected_for_figure == TRUE], file.path(tables, "V16_vs_selected_role_spearman_descriptive.csv"), na = "NA")
fwrite(equal_stroma, file.path(tables, "all_candidates_epithelial_vs_equal_weight_stromal_role_means.csv"))
fwrite(pairwise_stroma, file.path(tables, "all_candidates_epithelial_vs_each_stromal_role.csv"))
fwrite(equal_stroma[selected_for_figure == TRUE], file.path(tables, "epithelial_vs_equal_weight_stromal_role_means.csv"))
fwrite(pairwise_stroma[selected_for_figure == TRUE], file.path(tables, "epithelial_vs_each_stromal_role.csv"))
fwrite(rbindlist(control_records), file.path(audit, "KEGG_reproduction_checks.csv"))
fwrite(rbindlist(input_records), file.path(audit, "input_objects_verified.csv"))
protected_after <- setNames(vapply(protected_paths, md5, character(1)), protected_paths)
stopifnot(identical(protected_before, protected_after), sha256(selection_path) == selection_sha256)
fwrite(data.table(path = protected_paths, md5_before = unname(protected_before), md5_after = unname(protected_after), unchanged = TRUE),
  file.path(audit, "protected_input_hashes.csv"))
capture.output(sessionInfo(), file = file.path(audit, "R_sessionInfo.txt"))
report <- list(status = "completed", purpose = "Score all prespecified candidates and export the six choices fixed before scoring; no figure output",
  started = started, completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), script = script, script_md5 = md5(script),
  candidate_specifications = normalizePath(spec_path), candidate_specifications_md5 = protected_before[[spec_path]],
  preselection_record = normalizePath(selection_path), selection_SHA256_before = selection_sha256, selection_SHA256_after = sha256(selection_path),
  candidate_count = nrow(definitions), definitions = definitions, selected = selected_definitions,
  membership_policy = "Exact official complete GMT and JSON members, verified against manifest; no alias guessing, member tuning, overlap removal, or selection changes based on resulting scores",
  normalization = "log1p(count / complete original per-cell library total * 10000); zero-total genes removed; complete expression background retained",
  gsva = list(kcdf = "Gaussian", minSize = 3L, maxSize = "Inf", seed = 1L, BPPARAM = "BiocParallel::SerialParam",
    other_parameters = "Same installed GSVA 2.6.6 defaults as V16/V19; parameter object saved for every cohort", package_version = as.character(packageVersion("GSVA"))),
  cell_selection = "Same frozen RDS and frozen_role; only USAHCC normal cholangiocytes excluded before scoring",
  aggregation = "Mean per-cell GSVA within each original role; sample-SD row Z across all original roles within each cohort",
  display = "ChinaHCC normal cholangiocyte omitted only after standardization; retained in full tables and contrasts with display flag",
  control = list(standard_name = "KEGG hsa04151", complete_member_count = 362L, tolerance = 1e-10,
    maximum_per_cell_difference = max(rbindlist(control_records)[comparison == "per_cell"]$maximum_absolute_difference),
    maximum_role_mean_or_Z_difference = max(rbindlist(control_records)[comparison != "per_cell"]$maximum_absolute_difference)),
  contrasts = "Every normal epithelial and malignant role minus the equal-weight arithmetic mean of stromal-role means, and minus each individual stromal role; all original role labels retained; no P values",
  previous_comparison = "Each candidate mapped by its prespecified topic to the V16 label recorded in the frozen six-selection input",
  interpretation = "RNA gene-set enrichment; not direct phosphorylation or independent ligand-specific activity",
  row_order = row_order, input_objects = rbindlist(input_records), coverage_summary = coverage_summary,
  checks = list(official_GMT_JSON_and_manifest_verified = TRUE, all_candidates_retained = TRUE,
    six_selection_hash_unchanged = TRUE, complete_membership_unchanged = TRUE, all_frozen_MD5_match = TRUE,
    all_protected_inputs_unchanged = TRUE, KEGG_reproduction_pass = TRUE, all_scores_finite = TRUE,
    no_P_values = TRUE, no_figures_written = TRUE))
jsonlite::write_json(report, file.path(audit, "calculation_notes.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16, na = "null")
message("[verified] All candidate scores retained; fixed selected six exported; KEGG reproduced; selection and protected inputs unchanged")
