# Task 2: joint Harmony re-embedding with normal epithelial roots added.
# For each of 3 cohorts: take the existing epithelial+stromal joint frozen
# object (5 roles incl. malignant_epithelial_EMT_root), add the Task-1-frozen
# normal epithelial root(s).  HCC cohorts receive normal hepatocytes and any
# available normal cholangiocytes; CCA receives normal cholangiocytes only.
# Normal cells are subset from the raw expression
# objects, unify the 5 meta columns, merge, JoinLayers, then rerun the
# v4_step6zzf/g/l Harmony recipe. Mirrors that script's structure exactly.
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(data.table)
})
set.seed(20260901)

TABLES_DIR <- "emt_invivo_bodyfig_20260901/tables"
OUT_DIR <- "emt_invivo_bodyfig_20260901"
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

keep_meta <- c("cohort", "batch_id", "frozen_role", "percent.mt", "nCount_RNA")

# --- Lazy, cached loader for the 19G curate_v1 fallback object ---------------
.curate_v1_cache <- new.env()
.count1_cache <- new.env()
get_curate_v1 <- function() {
  if (is.null(.curate_v1_cache$obj)) {
    cat("[load] results_v3/scRNA_harmony_allcells_curate_v1.rds (19G, fallback source) ...\n")
    .curate_v1_cache$obj <- readRDS("results_v3/scRNA_harmony_allcells_curate_v1.rds")
  }
  .curate_v1_cache$obj
}
drop_curate_v1 <- function() {
  if (!is.null(.curate_v1_cache$obj)) {
    .curate_v1_cache$obj <- NULL
    gc(verbose = FALSE)
  }
}
get_count1 <- function() {
  if (is.null(.count1_cache$obj)) {
    cat("[load] data/GSE149614_HBV_HCC/count1_sparse.rds (sparse fallback) ...\n")
    .count1_cache$obj <- readRDS("data/GSE149614_HBV_HCC/count1_sparse.rds")
  }
  .count1_cache$obj
}

# --- Build one normal-root subset with unified 5 meta columns ---------------
# Checks the primary expression source for barcode coverage; if incomplete,
# falls back to curate_v1 (per brief: check sum(barcode %in% colnames(obj)),
# print it, decide which source to use).
build_normal_sub <- function(tag, cohort_label, role_name, normal_tab,
                              primary_obj, primary_name, target_n) {
  match_primary <- sum(normal_tab$barcode %in% colnames(primary_obj))
  cat(sprintf("[%s/%s] primary source match (%s): %d/%d\n",
              tag, role_name, primary_name, match_primary, nrow(normal_tab)))
  if (match_primary == nrow(normal_tab)) {
    src_obj <- primary_obj
    src_name <- primary_name
  } else {
    # The China-HCC primary object is a filtered Seurat object and can miss
    # some frozen normal cells.  Use the compact raw sparse matrix first;
    # loading the historical 19-Gb curate_v1 object is unnecessary and can
    # exhaust memory on this supplemental run.
    if (grepl("GSE149614_HBV_HCC", primary_name, fixed = TRUE)) {
      raw <- get_count1()
      match_raw <- sum(normal_tab$barcode %in% colnames(raw))
      cat(sprintf("[%s/%s] count1_sparse match: %d/%d\n", tag, role_name, match_raw, nrow(normal_tab)))
      stopifnot(match_raw == nrow(normal_tab))
      src_obj <- CreateSeuratObject(counts = raw[, normal_tab$barcode, drop = FALSE],
                                    project = "ChinaHCC_normal")
      # Preserve sample-level batches for Harmony (barcode prefix is the
      # donor/sample identifier in the HBV-HCC matrix).
      src_obj$orig.ident <- sub("_[^_]+$", "", colnames(src_obj))
      src_name <- "data/GSE149614_HBV_HCC/count1_sparse.rds"
    } else {
      cat(sprintf("[%s/%s] primary source incomplete -> fallback to results_v3/scRNA_harmony_allcells_curate_v1.rds\n",
                  tag, role_name))
      curate <- get_curate_v1()
      match_curate <- sum(normal_tab$barcode %in% colnames(curate))
      cat(sprintf("[%s/%s] curate_v1 match: %d/%d\n", tag, role_name, match_curate, nrow(normal_tab)))
      stopifnot(match_curate == nrow(normal_tab))
      src_obj <- curate
      src_name <- "results_v3/scRNA_harmony_allcells_curate_v1.rds"
    }
  }

  sub <- subset(src_obj, cells = normal_tab$barcode)
  sub$frozen_role <- setNames(normal_tab$frozen_role, normal_tab$barcode)[colnames(sub)]
  sub$cohort <- cohort_label
  DefaultAssay(sub) <- "RNA"
  sub[["percent.mt"]] <- PercentageFeatureSet(sub, pattern = "^MT-")

  # batch_id: prefer orig.ident unless it is a single constant value (i.e. not
  # informative for Harmony batch correction), then fall back to a
  # patient/sample column, then a single-batch placeholder.
  md <- sub@meta.data
  if ("orig.ident" %in% colnames(md) && length(unique(as.character(md$orig.ident))) > 1) {
    sub$batch_id <- as.character(md$orig.ident)
    batch_col <- "orig.ident"
  } else if ("S_ID" %in% colnames(md) && !all(is.na(md$S_ID))) {
    sub$batch_id <- as.character(md$S_ID)
    batch_col <- "S_ID"
  } else if ("orig.ident" %in% colnames(md)) {
    sub$batch_id <- as.character(md$orig.ident)
    batch_col <- "orig.ident (constant)"
  } else {
    sub$batch_id <- paste0(tag, "_normal")
    batch_col <- "(none available; single-batch placeholder)"
  }
  cat(sprintf("[%s/%s] batch_id column chosen: %s  (source=%s)\n", tag, role_name, batch_col, src_name))

  sub@meta.data <- sub@meta.data[, keep_meta, drop = FALSE]
  stopifnot(ncol(sub) == target_n)
  sub
}

# --- Per-cohort driver --------------------------------------------------------
run_cohort <- function(tag, cohort_label, existing_path,
                        expr_source_fn, expr_source_name,
                        normal_specs, target_total) {
  cat("\n===== ", tag, " =====\n", sep = "")

  existing <- readRDS(existing_path)
  DefaultAssay(existing) <- "RNA"
  stopifnot(all(keep_meta %in% colnames(existing@meta.data)))
  existing@meta.data <- existing@meta.data[, keep_meta, drop = FALSE]
  cat("[existing role_counts] (n=", ncol(existing), ")\n", sep = "")
  print(table(existing$frozen_role))

  primary_obj <- expr_source_fn()
  normal_subs <- lapply(normal_specs, function(spec) {
    tab <- fread(spec$csv)
    stopifnot(nrow(tab) == spec$target_n, !anyDuplicated(tab$barcode))
    # Some barcodes appear in both the Task-1 normal-root pool and the
    # existing frozen malignant_epithelial_EMT_root pool (differing
    # malignancy-classification pipelines: meta_slim celltype vs. the
    # copykat/inferCNV-based EMT-root freeze). The existing frozen object's
    # label takes precedence; drop these barcodes from the incoming normal
    # pool so the zero-overlap invariant holds, and report the deviation.
    conflict <- intersect(tab$barcode, colnames(existing))
    if (length(conflict) > 0L) {
      cat(sprintf("[%s/%s] DEVIATION: %d barcode(s) overlap existing frozen roles, dropped from normal pool: %s\n",
                  tag, spec$role, length(conflict), paste(conflict, collapse = ", ")))
      tab <- tab[!barcode %in% conflict]
    }
    build_normal_sub(tag, cohort_label, spec$role, tab, primary_obj, expr_source_name, nrow(tab))
  })
  rm(primary_obj)

  normal_colnames <- unlist(lapply(normal_subs, colnames))
  stopifnot(length(intersect(colnames(existing), normal_colnames)) == 0L)

  merged <- merge(existing, normal_subs)
  merged <- JoinLayers(merged)
  actual_total <- ncol(merged)
  if (actual_total != target_total) {
    cat(sprintf("[%s] DEVIATION: total cells %d != brief target %d (due to dropped overlap barcodes)\n",
                tag, actual_total, target_total))
  }
  role_counts <- table(merged$frozen_role)
  cat("[merged role_counts] (n=", ncol(merged), ")\n", sep = "")
  print(role_counts)

  rm(existing, normal_subs)
  gc(verbose = FALSE)

  # --- v4_step6zzf Harmony recipe on the merged object -----------------------
  DefaultAssay(merged) <- "RNA"
  merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)
  merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2500,
                                 verbose = FALSE)
  hvg <- VariableFeatures(merged)
  hvg <- hvg[!grepl("^(MT-|RPS|RPL|HBA|HBB|IG[HKL])", hvg)]
  VariableFeatures(merged) <- hvg
  merged <- ScaleData(merged, features = hvg,
                      vars.to.regress = c("nCount_RNA", "percent.mt"),
                      verbose = FALSE)
  merged <- RunPCA(merged, features = hvg, npcs = 30, verbose = FALSE)
  merged <- RunHarmony(merged, group.by.vars = "batch_id", reduction.use = "pca",
                       dims.use = 1:30, verbose = FALSE)
  stopifnot("harmony" %in% Reductions(merged), ncol(Embeddings(merged, "harmony")) == 30L)
  cat("[harmony dims]", ncol(Embeddings(merged, "harmony")), "\n")

  out_rds <- file.path(OUT_DIR, sprintf("S2_%s_joint_with_normal_roles_frozen.rds", tag))
  saveRDS(merged, out_rds)
  input_tab <- data.table(
    barcode = colnames(merged),
    cohort = as.character(merged$cohort),
    frozen_role = as.character(merged$frozen_role)
  )
  out_csv <- file.path(TABLES_DIR, sprintf("S2_%s_joint_with_normal_roles_input.csv", tag))
  fwrite(input_tab, out_csv)
  cat("[saved]", out_rds, "\n")
  cat("[saved]", out_csv, "\n")

  rm(merged)
  gc(verbose = FALSE)
  role_counts
}

# --- Cohort plans --------------------------------------------------------------
plans <- list(
  list(
    tag = "ChinaHCC", cohort_label = "China-HCC",
    existing_path = "results_v4/S6_ChinaHCC_joint_epithelial_stromal_roles_frozen.rds",
    expr_source_fn = function() readRDS("data/GSE149614_HBV_HCC/seu1.rds"),
    expr_source_name = "data/GSE149614_HBV_HCC/seu1.rds",
    normal_specs = list(
      list(csv = file.path(TABLES_DIR, "S1_ChinaHCC_normal_hepatocyte_root_frozen_by_cell.csv"),
           role = "normal_hepatocyte_root", target_n = 150L),
      list(csv = file.path(TABLES_DIR, "S1_ChinaHCC_normal_cholangiocyte_root_frozen_by_cell.csv"),
           role = "normal_cholangiocyte_root", target_n = 29L)
    ),
    target_total = 705L
  ),
  list(
    tag = "USAHCC", cohort_label = "USA-HCC",
    existing_path = "results_v4/S6_USAHCC_joint_epithelial_stromal_roles_frozen.rds",
    expr_source_fn = function() readRDS("data/GSE151530_nonHBV_HCC/seu2.rds"),
    expr_source_name = "data/GSE151530_nonHBV_HCC/seu2.rds",
    normal_specs = list(
      list(csv = file.path(TABLES_DIR, "S1_USAHCC_normal_hepatocyte_root_frozen_by_cell.csv"),
           role = "normal_hepatocyte_root", target_n = 150L),
      list(csv = file.path(TABLES_DIR, "S1_USAHCC_normal_cholangiocyte_root_frozen_by_cell.csv"),
           role = "normal_cholangiocyte_root", target_n = 4L)
    ),
    target_total = 458L
  ),
  list(
    tag = "CCAcombined", cohort_label = "CCA-combined",
    existing_path = "results_v4/S6_CCAcombined_joint_epithelial_stromal_roles_frozen.rds",
    expr_source_fn = function() readRDS("results_v4/S6_CCA_diploid_reclustered.rds"),
    expr_source_name = "results_v4/S6_CCA_diploid_reclustered.rds",
    normal_specs = list(
      list(csv = file.path(TABLES_DIR, "S1_CCAcombined_normal_cholangiocyte_root_frozen_by_cell.csv"),
           role = "normal_cholangiocyte_root", target_n = 150L)
    ),
    target_total = 369L
  )
)

# Set EMT_COHORT=ChinaHCC (or USAHCC/CCAcombined) to rebuild one cohort at a
# time while preserving the frozen inputs for the other cohorts.
cohort_filter <- Sys.getenv("EMT_COHORT", "")
if (nzchar(cohort_filter)) {
  plans <- Filter(function(p) identical(p$tag, cohort_filter), plans)
  if (length(plans) == 0L) stop("Unknown EMT_COHORT: ", cohort_filter)
}

all_role_counts <- list()
for (p in plans) {
  all_role_counts[[p$tag]] <- run_cohort(
    tag = p$tag, cohort_label = p$cohort_label, existing_path = p$existing_path,
    expr_source_fn = p$expr_source_fn, expr_source_name = p$expr_source_name,
    normal_specs = p$normal_specs, target_total = p$target_total
  )
  gc(verbose = FALSE)
}
drop_curate_v1()

cat("\n===== SUMMARY =====\n")
for (tag in names(all_role_counts)) {
  cat("--", tag, "--\n")
  print(all_role_counts[[tag]])
}
