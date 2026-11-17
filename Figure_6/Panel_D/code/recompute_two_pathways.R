#!/usr/bin/env Rscript
# Recompute only the two supplied complete KEGG gene sets. Existing TGF, HGF,
# EGF and IGF results are copied from the frozen, verified five-pathway tables.
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")
packages <- c("Seurat", "Matrix", "data.table", "GSVA", "BiocParallel", "matrixStats", "jsonlite")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages; no installation attempted: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table); library(GSVA)})

arg <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", arg[grepl("^--file=", arg)]), mustWork = TRUE)
work <- dirname(dirname(script))
root <- dirname(work)
inputs <- file.path(work, "inputs")
tables <- file.path(work, "tables")
audit <- file.path(work, "audit")
for (path in c(tables, audit)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
archive <- "/Users/sodajay/Desktop/Deutsch_Project/03_EMT/archive/EMT_reproduce_final_20260910/emt_invivo_bodyfig_20260901"
old_tables <- file.path(root, "V6_最终图_20261024", "audit", "IGF")
started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
md5 <- function(path) unname(tools::md5sum(path))

definitions <- list(
  Insulin = list(file = file.path(inputs, "KEGG_hsa04910_Insulin.gmt"), id = "hsa04910"),
  `PI3K-AKT` = list(file = file.path(inputs, "KEGG_hsa04151_PI3K_AKT.gmt"), id = "hsa04151")
)
gene_sets <- list()
definition_info <- list()
membership_by_pathway <- list()
membership_path <- file.path(inputs, "KEGG_pathway_membership.csv")
if (!file.exists(membership_path)) stop("Required official membership table is not ready: ", membership_path)
membership_input_md5 <- md5(membership_path)
official_membership <- fread(membership_path, colClasses = c(ncbi_gene_id = "character"))
stopifnot(all(c("pathway", "kegg_pathway_id", "ncbi_gene_id", "gene_symbol", "scoring_identifier") %in% names(official_membership)))
official_membership[is.na(gene_symbol), gene_symbol := ""]
for (name in names(definitions)) {
  def <- definitions[[name]]
  if (!file.exists(def$file)) stop("Required GMT is not ready: ", def$file)
  lines <- readLines(def$file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  stopifnot(length(lines) == 1L)
  fields <- strsplit(lines, "\t", fixed = TRUE)[[1]]
  stopifnot(length(fields) >= 5L, any(grepl(def$id, fields[1:2], fixed = TRUE)))
  genes <- fields[-c(1, 2)]
  stopifnot(all(nzchar(genes)), !anyDuplicated(genes), all(genes == trimws(genes)))
  member_rows <- official_membership[official_membership$pathway == name]
  stopifnot(nrow(member_rows) == length(genes), !anyDuplicated(member_rows$ncbi_gene_id),
    !anyDuplicated(member_rows$scoring_identifier), all(member_rows$kegg_pathway_id == def$id),
    setequal(member_rows$scoring_identifier, genes))
  member_rows <- member_rows[match(genes, member_rows$scoring_identifier)]
  has_symbol <- nzchar(member_rows$gene_symbol)
  stopifnot(identical(member_rows$scoring_identifier,
    ifelse(has_symbol, member_rows$gene_symbol, paste0("hsa:", member_rows$ncbi_gene_id))))
  gene_sets[[name]] <- genes
  membership_by_pathway[[name]] <- data.table(pathway = name, kegg_id = def$id,
    official_kegg_gene_id = paste0("hsa:", member_rows$ncbi_gene_id),
    ncbi_gene_id = member_rows$ncbi_gene_id, scoring_identifier = genes,
    gene_symbol = member_rows$gene_symbol, has_gene_symbol = has_symbol,
    mapping_status = ifelse(has_symbol, "official_symbol_available", "official_id_without_symbol"))
  definition_info[[name]] <- list(pathway = name, kegg_id = def$id,
    gmt_path = normalizePath(def$file), gmt_name = fields[1], gmt_description = fields[2],
    gmt_md5 = md5(def$file), complete_pathway_member_count = length(genes),
    mapped_gene_symbol_count = sum(has_symbol), unmapped_official_id_count = sum(!has_symbol),
    scoring_identifiers = genes)
}
shared_genes <- intersect(gene_sets$Insulin, gene_sets[["PI3K-AKT"]])
membership <- rbindlist(membership_by_pathway)
membership[, shared_between_Insulin_and_PI3K_AKT := scoring_identifier %in% shared_genes]
fwrite(membership, file.path(tables, "complete_gene_set_membership.csv"))

cohorts <- c(ChinaHCC = "HCC-GSE149614", USAHCC = "HCC-GSE151530", CCAcombined = "CCA-GSE138709_GSE151530")
expected_md5 <- c(ChinaHCC = "06910d7e500ac6da7526787396966e96", USAHCC = "d621adf1a8b5bc081a0ad938cdbb9329", CCAcombined = "e19b2ca07f12f2985b37cdc256b48d50")
expected_scored_cells <- c(ChinaHCC = 705L, USAHCC = 454L, CCAcombined = 369L)
groups <- list(
  ChinaHCC = c("normal_cholangiocyte_root", "normal_hepatocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "activated_HSC_intermediate", "portal_fibroblast_root", "fibrotic_MFB_terminal"),
  USAHCC = c("normal_hepatocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal"),
  CCAcombined = c("normal_cholangiocyte_root", "malignant_epithelial_EMT_root", "qHSC_root", "portal_fibroblast_like_root", "fibrotic_MFB_terminal")
)
labels <- c(normal_cholangiocyte_root = "Normal cholangiocyte", normal_hepatocyte_root = "Normal hepatocyte", malignant_epithelial_EMT_root = "Malignant epithelial", qHSC_root = "qHSC", activated_HSC_intermediate = "Activated HSC", portal_fibroblast_root = "Portal fibroblast", portal_fibroblast_like_root = "Portal fibroblast-like", fibrotic_MFB_terminal = "Fibrotic MFB")
row_order <- c("TGF-β", "PI3K-AKT", "Insulin", "HGF", "EGF", "IGF")
source_rows <- c(`TGF-β` = "TGF", HGF = "HGF-AKT", EGF = "EGF-RAS", IGF = "IGF")
input_records <- coverage_records <- role_records <- corr_records <- preserved_records <- list()

safe_spearman <- function(x, y) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)))
  if (length(x) < 3L) return(list(rho = NA_real_, status = "fewer_than_3_cells"))
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(list(rho = NA_real_, status = "constant_score"))
  list(rho = unname(stats::cor(x, y, method = "spearman")), status = "ok")
}

for (key in names(cohorts)) {
  message("[start] ", key)
  object_path <- file.path(archive, paste0("S2_", key, "_joint_with_normal_roles_frozen.rds"))
  before_md5 <- md5(object_path)
  stopifnot(identical(before_md5, unname(expected_md5[[key]])))
  obj <- readRDS(object_path)
  original_cells <- ncol(obj)
  if (key == "USAHCC") obj <- obj[, as.character(obj$frozen_role) != "normal_cholangiocyte_root"]
  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
  role <- as.character(obj$frozen_role)
  stopifnot(ncol(cnt) == expected_scored_cells[[key]], identical(colnames(cnt), colnames(obj)),
    !anyDuplicated(rownames(cnt)), !anyDuplicated(colnames(cnt)), !anyNA(role), setequal(unique(role), groups[[key]]))
  libraries <- Matrix::colSums(cnt)
  libraries[libraries <= 0] <- 1
  keep <- Matrix::rowSums(cnt) > 0
  # Denominators include all genes in each original cell, as in the frozen run.
  expr <- as.matrix(log1p(t(t(cnt[keep, , drop = FALSE]) / libraries) * 1e4))
  variable <- rownames(expr)[matrixStats::rowVars(expr) > 0]
  used_sets <- lapply(gene_sets, intersect, y = rownames(expr))
  stopifnot(all(vapply(used_sets, length, integer(1)) >= 3L))
  for (name in names(gene_sets)) {
    genes <- gene_sets[[name]]
    member_rows <- membership_by_pathway[[name]]
    coverage_records[[paste(key, name)]] <- data.table(
      cohort_key = key, cohort = cohorts[[key]], pathway = name,
      kegg_id = definitions[[name]]$id,
      official_kegg_gene_id = member_rows$official_kegg_gene_id,
      ncbi_gene_id = member_rows$ncbi_gene_id, scoring_identifier = genes,
      gene_symbol = member_rows$gene_symbol, has_gene_symbol = member_rows$has_gene_symbol,
      mapping_status = member_rows$mapping_status, in_complete_gene_set = TRUE,
      present_in_counts = genes %in% rownames(cnt),
      detected_in_any_cell = genes %in% rownames(cnt)[keep],
      nonzero_normalized_expression = genes %in% rownames(expr),
      variable = genes %in% variable,
      supplied_to_GSVA = genes %in% used_sets[[name]],
      shared_between_Insulin_and_PI3K_AKT = genes %in% shared_genes)
  }

  set.seed(1)
  param <- GSVA::gsvaParam(expr, used_sets, kcdf = "Gaussian", minSize = 3L, maxSize = Inf)
  capture.output(show(param), file = file.path(audit, paste0(key, "_GSVA_parameter_object.txt")))
  scores <- as.matrix(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
  stopifnot(identical(rownames(scores), names(gene_sets)), identical(colnames(scores), colnames(cnt)), all(is.finite(scores)))
  displayed <- !(key == "ChinaHCC" & groups[[key]] == "normal_cholangiocyte_root")
  per_cell <- data.table(cohort_key = key, cohort = cohorts[[key]], barcode = colnames(cnt),
    role_id = role, role_label = unname(labels[role]), displayed_in_figure6 = role %in% groups[[key]][displayed],
    Insulin = as.numeric(scores["Insulin", ]), `PI3K-AKT` = as.numeric(scores["PI3K-AKT", ]))
  fwrite(per_cell, file.path(tables, paste0(key, "_two_pathway_per_cell_scores.csv")))
  means <- vapply(groups[[key]], function(g) rowMeans(scores[, role == g, drop = FALSE]), numeric(length(gene_sets)))
  dimnames(means) <- list(names(gene_sets), unname(labels[groups[[key]]]))
  new_z <- t(scale(t(means)))
  new_z[!is.finite(new_z)] <- 0
  fwrite(as.data.table(means, keep.rownames = "pathway"), file.path(tables, paste0(key, "_two_pathway_role_means.csv")))
  fwrite(as.data.table(new_z, keep.rownames = "pathway"), file.path(tables, paste0(key, "_two_pathway_role_zscores.csv")))

  expected_columns <- c("pathway", unname(labels[groups[[key]]]))
  old <- list()
  old_paths <- old_hashes <- c()
  for (kind in c("by_role", "zscore")) {
    source_path <- file.path(old_tables, paste0("fig6_", key, "_signaling_GSVA_activity_5pathways_", kind, ".csv"))
    old_paths[kind] <- source_path
    old_hashes[kind] <- md5(source_path)
    source_table <- fread(source_path)
    stopifnot(identical(names(source_table), expected_columns),
      identical(source_table$pathway, c("TGF", "INSULIN-AKT", "HGF-AKT", "IGF", "EGF-RAS")))
    old[[kind]] <- as.matrix(source_table[, -1])
    rownames(old[[kind]]) <- source_table$pathway
  }
  # Validate the previous normalization, but preserve the existing four rows.
  old_z_check <- t(scale(t(old$by_role)))
  old_z_check[!is.finite(old_z_check)] <- 0
  stopifnot(max(abs(old_z_check - old$zscore)) < 1e-12)
  combined <- list()
  for (kind in c("by_role", "zscore")) {
    result <- matrix(NA_real_, nrow = length(row_order), ncol = length(groups[[key]]),
      dimnames = list(row_order, expected_columns[-1]))
    result[c("Insulin", "PI3K-AKT"), ] <- if (kind == "by_role") means else new_z
    for (name in names(source_rows)) result[name, ] <- old[[kind]][source_rows[[name]], ]
    stopifnot(all(is.finite(result)))
    combined[[kind]] <- result
    destination <- file.path(tables, paste0("fig6_", key, "_signaling_GSVA_activity_6pathways_", kind, ".csv"))
    fwrite(as.data.table(result, keep.rownames = "pathway"), destination)
    fwrite(as.data.table(result[, displayed, drop = FALSE], keep.rownames = "pathway"),
      file.path(tables, paste0("fig6_", key, "_signaling_GSVA_activity_6pathways_", kind, "_displayed.csv")))
    reread <- fread(destination)
    for (name in names(source_rows)) {
      delta <- max(abs(as.numeric(reread[reread$pathway == name, -1]) - old[[kind]][source_rows[[name]], ]))
      stopifnot(delta < 1e-14, md5(old_paths[[kind]]) == old_hashes[[kind]])
      preserved_records[[paste(key, kind, name)]] <- data.table(cohort_key = key, kind = kind,
        output_pathway = name, source_pathway = source_rows[[name]],
        source_file = old_paths[[kind]], source_md5 = old_hashes[[kind]],
        copied_value_count = length(groups[[key]]), maximum_serialized_numeric_difference = delta,
        source_unchanged = TRUE, pathway_recomputed = FALSE)
    }
  }
  role_records[[key]] <- rbindlist(lapply(row_order, function(name) data.table(
    cohort_key = key, cohort = cohorts[[key]], pathway = name,
    role_id = groups[[key]], role_label = unname(labels[groups[[key]]]),
    n_cells = as.integer(table(factor(role, levels = groups[[key]]))),
    mean_gsva = as.numeric(combined$by_role[name, ]),
    within_cohort_role_zscore = as.numeric(combined$zscore[name, ]),
    displayed_in_figure6 = displayed, recomputed = name %in% names(gene_sets))))

  for (g in c("ALL_CELLS", groups[[key]])) {
    selected <- if (g == "ALL_CELLS") rep(TRUE, length(role)) else role == g
    correlation <- safe_spearman(scores["Insulin", selected], scores["PI3K-AKT", selected])
    corr_records[[paste(key, g)]] <- data.table(cohort_key = key, cohort = cohorts[[key]],
      scope = if (g == "ALL_CELLS") "all_scored_cells" else "within_role", role_id = g,
      role_label = if (g == "ALL_CELLS") "All scored cells" else unname(labels[g]),
      n_cells = sum(selected), spearman_rho = correlation$rho, status = correlation$status,
      p_value_computed = FALSE)
  }
  after_md5 <- md5(object_path)
  stopifnot(identical(after_md5, before_md5))
  input_records[[key]] <- data.table(cohort_key = key, cohort = cohorts[[key]],
    source = object_path, expected_md5 = expected_md5[[key]], input_md5_before = before_md5,
    input_md5_after = after_md5, input_unchanged = TRUE, original_cells = original_cells,
    cells_scored = ncol(cnt), excluded_normal_cholangiocyte_cells = original_cells - ncol(cnt),
    roles_scored = length(groups[[key]]), roles_displayed = sum(displayed),
    genes_in_counts = nrow(cnt), nonzero_genes = nrow(expr), variable_genes = length(variable))
  message("[done] ", key, ": ", ncol(cnt), " cells, Insulin ", length(used_sets$Insulin), "/", length(gene_sets$Insulin),
    "; PI3K-AKT ", length(used_sets[["PI3K-AKT"]]), "/", length(gene_sets[["PI3K-AKT"]]))
  rm(obj, cnt, expr, scores, param, per_cell)
  invisible(gc())
}

coverage <- rbindlist(coverage_records)
coverage_summary <- coverage[, .(complete_pathway_member_count = .N,
  mapped_gene_symbol_count = sum(has_gene_symbol), unmapped_official_id_count = sum(!has_gene_symbol),
  present_in_counts = sum(present_in_counts), detected_in_any_cell = sum(detected_in_any_cell),
  nonzero_normalized_expression = sum(nonzero_normalized_expression), variable = sum(variable),
  supplied_to_GSVA = sum(supplied_to_GSVA), shared_genes_in_complete_set = sum(shared_between_Insulin_and_PI3K_AKT),
  shared_variable_genes = sum(shared_between_Insulin_and_PI3K_AKT & variable)), by = .(cohort_key, cohort, pathway, kegg_id)]
fwrite(coverage, file.path(tables, "gene_coverage.csv"))
fwrite(coverage_summary, file.path(tables, "gene_coverage_summary.csv"))
fwrite(rbindlist(role_records), file.path(tables, "six_pathways_all_role_scores.csv"))
fwrite(rbindlist(corr_records), file.path(tables, "Insulin_PI3K_AKT_spearman_descriptive.csv"), na = "NA")
fwrite(rbindlist(input_records), file.path(audit, "input_objects_verified.csv"))
fwrite(rbindlist(preserved_records), file.path(audit, "four_existing_pathways_preserved.csv"))
for (name in names(definitions)) stopifnot(md5(definitions[[name]]$file) == definition_info[[name]]$gmt_md5)
stopifnot(md5(membership_path) == membership_input_md5)
capture.output(sessionInfo(), file = file.path(audit, "R_sessionInfo.txt"))
report <- list(
  status = "completed", started = started, completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script = script, only_recomputed_pathways = names(gene_sets),
  complete_standard_gene_sets = definition_info,
  membership_mapping = list(source_file = normalizePath(membership_path), source_md5 = membership_input_md5,
    complete_counts_unit = "Distinct official KEGG/NCBI gene IDs within each pathway",
    unresolved_members = membership[has_gene_symbol == FALSE],
    unresolved_member_policy = "Retain the official hsa:geneID as scoring_identifier; leave gene_symbol empty and mark official_id_without_symbol. Do not guess a symbol or manually delete the member."),
  shared_genes = list(count = length(shared_genes), symbols = shared_genes,
    jaccard = length(shared_genes) / length(union(gene_sets$Insulin, gene_sets[["PI3K-AKT"]])),
    retained_in_both_gene_sets = TRUE, pathway_overlap_pruning_performed = FALSE),
  normalization = "log1p(count / complete original per-cell library total * 10000); zero-total genes removed as in the original script",
  gsva = list(kcdf = "Gaussian", minSize = 3L, maxSize = "Inf", seed = 1L,
    BPPARAM = "BiocParallel::SerialParam", other_parameters = "Installed-version gsvaParam defaults; parameter-object snapshots saved per cohort",
    package_version = as.character(packageVersion("GSVA")),
    gene_set_handling = "Full GMT memberships archived unchanged; only genes present in the nonzero expression matrix supplied to GSVA; constant-gene filtering is handled by GSVA"),
  cell_selection = "Same frozen RNA counts and frozen_role as the verified original; normal_cholangiocyte_root excluded only from USAHCC before scoring",
  aggregation = "Mean per-cell GSVA score for each original role; sample-SD row z-score across all original roles within each cohort",
  figure_display = "ChinaHCC Normal cholangiocyte omitted only after scoring and row standardization; all original roles remain in the full tables",
  final_row_order = row_order,
  preserved_pathways = list(names = names(source_rows), source_directory = old_tables,
    means_and_zscores_copied = TRUE, recomputed = FALSE,
    maximum_serialized_numeric_difference = max(rbindlist(preserved_records)$maximum_serialized_numeric_difference)),
  correlation = list(method = "Spearman", scopes = c("all scored cells within each cohort", "each original role within each cohort"),
    descriptive_only = TRUE, p_values_computed = FALSE, cell_scores_not_patient_independent = TRUE),
  interpretation = "RNA-based gene-program enrichment, not direct evidence of receptor phosphorylation or attribution to a specific ligand; shared genes are retained",
  input_objects = rbindlist(input_records), coverage_summary = coverage_summary,
  checks = list(all_frozen_object_MD5_match = TRUE, original_input_files_unchanged = TRUE,
    four_existing_pathways_preserved = TRUE, complete_gene_memberships_unchanged = TRUE,
    new_score_values_finite = TRUE, no_old_figures_or_outputs_written = TRUE))
jsonlite::write_json(report, file.path(audit, "calculation_notes.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16, na = "null")
message("[verified] Two KEGG pathways computed; four existing pathway rows retained; all input hashes unchanged.")
