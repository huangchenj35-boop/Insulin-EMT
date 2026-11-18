#!/usr/bin/env Rscript
# Cook 2020-derived two-axis EMT landscape for V4 malignant epithelial cells.
# Y: MSigDB Hallmark EMT score, scored by AddModuleScore as in Cook et al.
# X: Cook et al. conserved downregulated epithelial-associated gene-set score.
# Pooled medians are used only as a descriptive four-quadrant display rule.

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

object_file <- "results_v3/scRNA_harmony_allcells_curate_v1.rds"
malignant_file <- "results_v4/tables/S5_malignant_cells_v4.csv"
geneset_file <- "results_backup_v1/03_EMT/geneSets.rds"
review_dir <- "results_v4/Figure_Review"
table_dir <- "results_v4/tables"

cohort_levels <- c("China-HCC", "USA-HCC", "CCA-combined")
figure_stems <- c(
  `China-HCC` = "Fig3a_Cook2020_two_axis_EMT_ChinaHCC",
  `USA-HCC` = "Fig3b_Cook2020_two_axis_EMT_USAHCC",
  `CCA-combined` = "Fig3c_Cook2020_two_axis_EMT_CCAcombined"
)
display_cap <- 6000L
sampling_seed <- 20260827L

state_levels <- c("Epithelial-high", "E/M-high", "EMT-high", "Low-both")
state_colors <- c(
  `Epithelial-high` = "#5AB4D0",
  `E/M-high` = "#F28E62",
  `EMT-high` = "#55BDAE",
  `Low-both` = "#9E9EAD"
)

# Supplementary Data 2, conserved_downregulated_genes (17 genes), Cook et al.
cook_epithelial <- c(
  "KRT18", "ATP1B1", "PLA2G16", "DHCR24", "HIST1H4C", "FOS",
  "KRT19", "FTL", "ELF3", "KRT8", "CDKN3", "KIAA0101", "GINS2",
  "HCFC1R1", "EPHX1", "ACAT2", "AGR2"
)
hallmark_emt <- unique(
  readRDS(geneset_file)$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
)
stopifnot(length(cook_epithelial) == 17L, length(hallmark_emt) == 200L)

dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

mal <- fread(malignant_file)
stopifnot(all(c("barcode", "cohort", "patient") %in% names(mal)))
stopifnot(!anyDuplicated(mal$barcode))
mal <- mal[cohort %in% cohort_levels]
mal[, cohort := factor(cohort, levels = cohort_levels)]

cat("Loading curated object and scoring Cook-derived EMT programs ...\n")
seu <- readRDS(object_file)
cells <- intersect(mal$barcode, colnames(seu))
stopifnot(length(cells) == nrow(mal))
data_all <- GetAssayData(seu, assay = "RNA", layer = "data")
pool <- rownames(data_all)

map_features <- function(genes) {
  ans <- data.table(source_gene = genes, assay_gene = genes)
  aliases <- c(
    CCN1 = "CYR61",
    CCN2 = "CTGF",
    KIAA0101 = "PCLAF",
    HIST1H4C = "H4C3"
  )
  for (old in names(aliases)) {
    new <- unname(aliases[old])
    ans[source_gene == old & !old %in% pool & new %in% pool,
        assay_gene := new]
  }
  ans[, available := assay_gene %in% pool]
  ans[]
}

epi_map <- map_features(cook_epithelial)
emt_map <- map_features(hallmark_emt)
epi_features <- unique(epi_map[available == TRUE, assay_gene])
emt_features <- unique(emt_map[available == TRUE, assay_gene])
if (length(epi_features) < 15L || length(emt_features) < 190L) {
  stop("Insufficient Cook-derived signature genes in the RNA assay")
}

# Exact one-module equivalent of Seurat::AddModuleScore defaults. Both modules
# use the same pooled malignant-cell expression universe and random seed.
score_add_module <- function(features, cells, data_all, pool,
                             nbin = 24L, ctrl = 100L, seed = 1L) {
  set.seed(seed)
  data_avg <- Matrix::rowMeans(data_all[pool, cells, drop = FALSE])
  data_avg <- data_avg[order(data_avg)]
  data_cut <- ggplot2::cut_number(
    data_avg + rnorm(length(data_avg)) / 1e30,
    n = nbin, labels = FALSE, right = FALSE
  )
  names(data_cut) <- names(data_avg)
  controls <- character()
  for (gene in features) {
    same_bin <- names(data_cut)[data_cut == data_cut[gene]]
    controls <- c(controls, sample(same_bin, size = ctrl, replace = FALSE))
  }
  controls <- unique(controls)
  feature_score <- Matrix::colMeans(data_all[features, cells, drop = FALSE])
  control_score <- Matrix::colMeans(data_all[controls, cells, drop = FALSE])
  list(score = feature_score - control_score, controls = controls)
}

epi <- score_add_module(epi_features, cells, data_all, pool)
emt <- score_add_module(emt_features, cells, data_all, pool)
rm(data_all, seu)
gc(verbose = FALSE)

cell_lookup <- mal[match(cells, barcode)]
stopifnot(identical(cell_lookup$barcode, cells))
d <- cell_lookup[, .(barcode, cohort, patient)]
d[, epithelial_score := unname(epi$score[barcode])]
d[, hallmark_EMT_score := unname(emt$score[barcode])]
stopifnot(all(is.finite(d$epithelial_score)), all(is.finite(d$hallmark_EMT_score)))

epi_cut <- median(d$epithelial_score)
emt_cut <- median(d$hallmark_EMT_score)
d[, quadrant := fifelse(
  epithelial_score > epi_cut & hallmark_EMT_score <= emt_cut,
  "Epithelial-high",
  fifelse(
    epithelial_score > epi_cut & hallmark_EMT_score > emt_cut,
    "E/M-high",
    fifelse(hallmark_EMT_score > emt_cut, "EMT-high", "Low-both")
  )
)]
d[, quadrant := factor(quadrant, levels = state_levels)]

summary <- d[, .(cells = .N), by = .(cohort, quadrant)]
summary[, percentage := 100 * cells / sum(cells), by = cohort]
summary[, epithelial_cutoff := epi_cut]
summary[, hallmark_EMT_cutoff := emt_cut]

fwrite(d, file.path(table_dir, "S5_Cook2020_two_axis_EMT_scores_by_cell.csv"))
fwrite(summary, file.path(table_dir, "S5_Cook2020_two_axis_EMT_quadrant_by_cohort.csv"))
epi_map[, program := "Cook_conserved_epithelial_associated"]
emt_map[, program := "MSigDB_Hallmark_EMT"]
fwrite(rbindlist(list(epi_map, emt_map), use.names = TRUE),
       file.path(table_dir, "S5_Cook2020_two_axis_EMT_signature_definition.csv"))
fwrite(rbindlist(list(
  data.table(program = "Cook_conserved_epithelial_associated",
             control_gene = epi$controls),
  data.table(program = "MSigDB_Hallmark_EMT", control_gene = emt$controls)
)), file.path(table_dir, "S5_Cook2020_two_axis_EMT_AddModuleScore_controls.csv"))

global_x <- as.numeric(quantile(d$epithelial_score, c(0.005, 0.995)))
global_y <- as.numeric(quantile(d$hallmark_EMT_score, c(0.005, 0.995)))

plot_one <- function(cohort_name) {
  z <- copy(d[cohort == cohort_name])
  zs <- copy(summary[cohort == cohort_name])
  set.seed(sampling_seed + match(cohort_name, cohort_levels))
  point_z <- z[sample.int(.N, min(.N, display_cap))]

  dx <- 0.025 * diff(global_x)
  dy <- 0.028 * diff(global_y)
  labels <- data.table(
    x = c(global_x[2] - dx, global_x[2] - dx,
          global_x[1] + dx, global_x[1] + dx),
    y = c(global_y[1] + dy, global_y[2] - dy,
          global_y[2] - dy, global_y[1] + dy),
    quadrant = factor(state_levels, levels = state_levels),
    hjust = c(1, 1, 0, 0),
    vjust = c(0, 1, 1, 0)
  )
  labels <- merge(labels, zs[, .(quadrant, percentage)], by = "quadrant",
                  sort = FALSE)
  labels[, label := sprintf("%s\n%.1f%%", as.character(quadrant), percentage)]

  p_main <- ggplot(point_z, aes(epithelial_score, hallmark_EMT_score)) +
    stat_density_2d(
      data = z, bins = 6, colour = "#424242",
      linewidth = 0.42, alpha = 0.42, show.legend = FALSE
    ) +
    geom_point(aes(colour = quadrant), size = 0.80, alpha = 0.32, stroke = 0) +
    geom_vline(xintercept = epi_cut, linetype = "dashed", linewidth = 0.5,
               colour = "#333333") +
    geom_hline(yintercept = emt_cut, linetype = "dashed", linewidth = 0.5,
               colour = "#333333") +
    geom_text(
      data = labels,
      aes(x, y, label = label, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE, size = 2.9, fontface = "bold", colour = "#333333"
    ) +
    scale_colour_manual(values = state_colors, breaks = state_levels,
                        name = NULL, drop = FALSE) +
    coord_cartesian(xlim = global_x, ylim = global_y, expand = FALSE) +
    labs(x = "Cook conserved epithelial-associated score",
         y = "Hallmark EMT score") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "#333333", linewidth = 0.5),
      axis.title = element_text(size = 10.2, face = "bold"),
      axis.text = element_text(size = 8.5, colour = "#333333"),
      legend.position = "bottom",
      legend.text = element_text(size = 8.5),
      plot.margin = margin(1, 1, 1, 1)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 1))

  p_top <- ggplot(z, aes(epithelial_score)) +
    geom_density(fill = "#5AB4D0", colour = "#2C7FB8", linewidth = 0.55,
                 alpha = 0.65, adjust = 0.85) +
    geom_vline(xintercept = epi_cut, linetype = "dashed", linewidth = 0.5) +
    coord_cartesian(xlim = global_x, expand = FALSE) +
    theme_void()

  p_right <- ggplot(z, aes(hallmark_EMT_score)) +
    geom_density(fill = "#55BDAE", colour = "#238B7E", linewidth = 0.55,
                 alpha = 0.65, adjust = 0.85) +
    geom_vline(xintercept = emt_cut, linetype = "dashed", linewidth = 0.5) +
    coord_flip(xlim = global_y, expand = FALSE) +
    theme_void()

  design <- "AAAB\nCCCD\nCCCD\nCCCD"
  p <- wrap_plots(A = p_top, B = plot_spacer(), C = p_main, D = p_right,
                  design = design, guides = "collect") +
    plot_annotation(
      title = sprintf("%s: Cook-derived EMT landscape", cohort_name),
      subtitle = sprintf(
        "%s malignant epithelial cells; continuous program scores with pooled global median guides",
        format(nrow(z), big.mark = ",")
      ),
      theme = theme(
        plot.title = element_text(size = 14.5, face = "bold"),
        plot.subtitle = element_text(size = 9.2, colour = "#555555"),
        plot.margin = margin(7, 8, 5, 7)
      )
    ) & theme(legend.position = "bottom")

  stem <- file.path(review_dir, unname(figure_stems[cohort_name]))
  ggsave(paste0(stem, ".pdf"), p, width = 7.2, height = 6.8,
         device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 7.2, height = 6.8,
         dpi = 300, bg = "white")
}

invisible(lapply(cohort_levels, plot_one))

audit <- data.table(
  epithelial_program = "Cook 2020 conserved downregulated epithelial-associated genes",
  epithelial_genes_total = length(cook_epithelial),
  epithelial_genes_available = length(epi_features),
  EMT_program = "MSigDB HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  EMT_genes_total = length(hallmark_emt),
  EMT_genes_available = length(emt_features),
  scoring_method = "Seurat AddModuleScore",
  nbin = 24L, ctrl = 100L, seed = 1L,
  scoring_universe = "pooled V4 malignant epithelial cells across all three cohorts",
  quadrant_rule = "pooled global medians; descriptive visualization only",
  reference = paste(
    "Cook DP, Vanderhyden BC. Nat Commun. 2020;11:2142.",
    "doi:10.1038/s41467-020-16066-2; PMID:32358524"
  )
)
fwrite(audit, file.path(table_dir, "S5_Cook2020_two_axis_EMT_method_audit.csv"))

cat("\n=== Cook-derived two-axis EMT landscape complete ===\n")
print(summary)
cat("Epithelial genes available:", length(epi_features), "/ 17\n")
cat("Hallmark EMT genes available:", length(emt_features), "/ 200\n")
