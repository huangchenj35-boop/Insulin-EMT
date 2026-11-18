# Shared helpers reproducing the Figure 2A E/M scoring (10abc_review_Cook2020_two_axis_EMT_v4.R).

load_reference <- function(up) {
  tab <- file.path(up, "results_v4/tables")
  sig <- fread(file.path(tab, "S5_Cook2020_two_axis_EMT_signature_definition.csv"))
  ctl <- fread(file.path(tab, "S5_Cook2020_two_axis_EMT_AddModuleScore_controls.csv"))
  saved <- fread(file.path(tab, "S5_Cook2020_two_axis_EMT_scores_by_cell.csv"))
  stopifnot(nrow(saved) == 29648L)
  epi <- "Cook_conserved_epithelial_associated"; emt <- "MSigDB_Hallmark_EMT"
  list(
    epi_features = unique(sig[program == epi & available == TRUE, assay_gene]),
    emt_features = unique(sig[program == emt & available == TRUE, assay_gene]),
    epi_controls = unique(ctl[program == epi, control_gene]),
    emt_controls = unique(ctl[program == emt, control_gene]),
    saved = saved,
    E_cut = median(saved$epithelial_score),
    M_cut = median(saved$hallmark_EMT_score)
  )
}

# Seurat LogNormalize: log1p(count / library size * 1e4), library size over all genes.
# CreateSeuratObject renames "_" to "-" in feature names; do the same.
lognorm <- function(cnt) {
  cnt <- as(cnt, "CsparseMatrix")
  rownames(cnt) <- gsub("_", "-", rownames(cnt))
  stopifnot(!anyDuplicated(rownames(cnt)))
  lib <- Matrix::colSums(cnt)
  cnt@x <- log1p(cnt@x / rep.int(lib, diff(cnt@p)) * 1e4)
  cnt
}

# Mean over the full gene list; genes absent from this matrix contribute zero
gene_mean <- function(x, genes) {
  g <- intersect(genes, rownames(x))
  Matrix::colSums(x[g, , drop = FALSE]) / length(genes)
}

score_cells <- function(x, ref) {
  data.table(barcode = colnames(x),
             epithelial_score = unname(gene_mean(x, ref$epi_features) - gene_mean(x, ref$epi_controls)),
             hallmark_EMT_score = unname(gene_mean(x, ref$emt_features) - gene_mean(x, ref$emt_controls)))
}

classify <- function(d, ref) {
  d[, state := fifelse(epithelial_score > ref$E_cut & hallmark_EMT_score <= ref$M_cut, "Epithelial high",
               fifelse(epithelial_score > ref$E_cut & hallmark_EMT_score > ref$M_cut, "E/M high",
               fifelse(hallmark_EMT_score > ref$M_cut, "Mesenchymal high", "Low both")))]
  d[]
}

validate_against_saved <- function(cnt, saved_rows) {
  s <- score_cells(lognorm(cnt), ref = get("ref", envir = parent.frame()))
  m <- merge(s, saved_rows[, .(barcode, E0 = epithelial_score, M0 = hallmark_EMT_score)], by = "barcode")
  stopifnot(nrow(m) == nrow(saved_rows))
  data.table(cells = nrow(m),
             max_abs_diff_E = max(abs(m$epithelial_score - m$E0)),
             max_abs_diff_M = max(abs(m$hallmark_EMT_score - m$M0)),
             cor_E = cor(m$epithelial_score, m$E0), cor_M = cor(m$hallmark_EMT_score, m$M0))
}
