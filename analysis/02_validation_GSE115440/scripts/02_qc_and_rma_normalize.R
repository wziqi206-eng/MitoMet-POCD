# ==============================================================================
# Script: 02_qc_and_rma_normalize.R
# Project: MitoMet-POCD / NeuroMitoMap
# Dataset: GSE115440 (Tier 2 validation, Surgery vs Control + Maresin1 rescue)
# Platform: GPL11533 Affymetrix Mouse Gene 1.1 ST Array (MoGene-1_1-st-v1)
# Purpose: QC and RMA normalization of 9 CEL.gz files
# Note   : This script does NOT run limma. That is Script 03.
#          This script does NOT touch GSE95426 (Tier 1 anchor, frozen).
# ==============================================================================

# ---- 0. Setup ----------------------------------------------------------------
cat("==============================================================\n")
cat(" GSE115440: QC + RMA normalization (Script 02)\n")
cat("==============================================================\n\n")

set.seed(42)

# ---- 1. Package management ---------------------------------------------------
cat(">>> Step 1: Checking required packages...\n")

cran_pkgs <- c("ggplot2", "pheatmap", "RColorBrewer", "reshape2")
bioc_pkgs <- c("oligo", "Biobase", "pd.mogene.1.1.st.v1")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("    Installing CRAN package: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("    Installing Bioconductor package: %s\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(oligo)
  library(Biobase)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
})

cat("    All packages loaded.\n\n")

# ---- 2. Define paths ---------------------------------------------------------
cat(">>> Step 2: Setting up paths...\n")

cel_dir       <- "data/raw/GSE115440/CEL"
metadata_path <- "data/metadata/GSE115440_metadata_verified.csv"

results_dir   <- "analysis/02_validation_GSE115440/results"
tables_dir    <- file.path(results_dir, "tables")
figures_dir   <- file.path(results_dir, "figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("    CEL directory:  %s\n", cel_dir))
cat(sprintf("    Metadata file:  %s\n", metadata_path))
cat(sprintf("    Tables output:  %s\n", tables_dir))
cat(sprintf("    Figures output: %s\n\n", figures_dir))

# ---- 3. Load and validate metadata -------------------------------------------
cat(">>> Step 3: Loading sample metadata...\n")

if (!file.exists(metadata_path)) {
  stop(sprintf("Metadata file not found: %s", metadata_path))
}

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

cat(sprintf("    Loaded metadata: %d rows x %d cols\n", nrow(metadata), ncol(metadata)))
cat("    Columns: ", paste(colnames(metadata), collapse = ", "), "\n", sep = "")

gsm_col <- NULL
for (col in colnames(metadata)) {
  if (any(grepl("^GSM[0-9]+$", as.character(metadata[[col]])))) {
    gsm_col <- col
    break
  }
}

if (is.null(gsm_col)) {
  stop("Could not find a GSM column in metadata.")
}

if (any(is.na(metadata[[gsm_col]]))) {
  stop("Metadata GSM column contains NA.")
}

if (any(duplicated(metadata[[gsm_col]]))) {
  stop("Metadata GSM column has duplicates.")
}

cat(sprintf("    GSM column detected:   %s\n", gsm_col))

group_candidates <- c("group", "Group", "condition", "Condition", "treatment", "Treatment")
group_col <- NULL

for (col in group_candidates) {
  if (col %in% colnames(metadata)) {
    group_col <- col
    break
  }
}

if (is.null(group_col)) {
  warning("No standard group column found; using second non-GSM column as fallback.")
  group_col <- setdiff(colnames(metadata), gsm_col)[1]
}

cat(sprintf("    Group column detected: %s\n\n", group_col))

# ---- 4. Locate and match CEL files -------------------------------------------
cat(">>> Step 4: Locating CEL files...\n")

if (!dir.exists(cel_dir)) {
  stop(sprintf("CEL directory not found: %s", cel_dir))
}

cel_files <- list.files(
  cel_dir,
  pattern = "\\.CEL\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

cat(sprintf("    Found %d CEL.gz files.\n", length(cel_files)))

if (length(cel_files) != 9) {
  stop(sprintf("Expected 9 CEL.gz files, found %d in %s", length(cel_files), cel_dir))
}

cel_basenames <- basename(cel_files)
cel_gsm <- regmatches(cel_basenames, regexpr("GSM[0-9]+", cel_basenames))

if (length(cel_gsm) != length(cel_files)) {
  stop("Some CEL filenames do not contain a GSM ID.")
}

if (any(duplicated(cel_gsm))) {
  stop("Duplicate GSM IDs detected across CEL filenames.")
}

cel_map <- data.frame(
  cel_file = cel_files,
  cel_base = cel_basenames,
  gsm = cel_gsm,
  stringsAsFactors = FALSE
)

missing_in_cel <- setdiff(metadata[[gsm_col]], cel_map$gsm)
missing_in_meta <- setdiff(cel_map$gsm, metadata[[gsm_col]])

if (length(missing_in_cel) > 0) {
  stop(sprintf("Metadata GSMs missing from CEL files: %s", paste(missing_in_cel, collapse = ", ")))
}

if (length(missing_in_meta) > 0) {
  stop(sprintf("CEL GSMs missing from metadata: %s", paste(missing_in_meta, collapse = ", ")))
}

cel_map <- cel_map[match(metadata[[gsm_col]], cel_map$gsm), ]
cel_files_ordered <- cel_map$cel_file

cat("    CEL <-> GSM matching successful:\n")
print(data.frame(
  GSM = cel_map$gsm,
  group = metadata[[group_col]],
  file = cel_map$cel_base,
  stringsAsFactors = FALSE
))
cat("\n")

# ---- 5. Read CEL files with oligo --------------------------------------------
cat(">>> Step 5: Reading CEL files with oligo::read.celfiles...\n")

raw_data <- tryCatch(
  oligo::read.celfiles(cel_files_ordered, verbose = FALSE),
  error = function(e) {
    stop(sprintf("Failed to read CEL files: %s", e$message))
  }
)

sampleNames(raw_data) <- cel_map$gsm
pData(raw_data)$gsm <- cel_map$gsm
pData(raw_data)$group <- metadata[[group_col]]

cat(sprintf("    Raw data: %d probes/features x %d samples\n", nrow(raw_data), ncol(raw_data)))
cat(sprintf("    Annotation used by oligo: %s\n\n", annotation(raw_data)))

# ---- 6. Color setup ----------------------------------------------------------
group_factor <- factor(pData(raw_data)$group)
n_groups <- length(levels(group_factor))
group_colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
sample_colors <- group_colors[as.integer(group_factor)]

# ---- 7. Boxplot BEFORE RMA ---------------------------------------------------
cat(">>> Step 6: Boxplot BEFORE RMA...\n")

png(
  file.path(figures_dir, "GSE115440_boxplot_before_rma.png"),
  width = 1800,
  height = 1200,
  res = 200
)

par(mar = c(8, 5, 4, 2))

tryCatch({
  oligo::boxplot(
    raw_data,
    main = "GSE115440: raw intensities before RMA",
    las = 2,
    col = sample_colors,
    ylab = "raw intensity"
  )
  legend("topright", legend = levels(group_factor), fill = group_colors, bty = "n", cex = 0.8)
}, error = function(e) {
  plot.new()
  title("Pre-RMA boxplot failed; see log.")
  cat(sprintf("    Warning: pre-RMA boxplot failed: %s\n", e$message))
})

dev.off()

cat("    Saved: GSE115440_boxplot_before_rma.png\n\n")

# ---- 8. RMA normalization ----------------------------------------------------
cat(">>> Step 7: Running RMA normalization with oligo::rma...\n")

eset <- oligo::rma(raw_data)

cat(sprintf("    RMA-normalized: %d probes/features x %d samples\n\n", nrow(eset), ncol(eset)))

# ---- 9. Export RMA expression matrix -----------------------------------------
cat(">>> Step 8: Exporting RMA expression matrix...\n")

expr_mat <- exprs(eset)

expr_df <- data.frame(
  probe_id = rownames(expr_mat),
  expr_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_out <- file.path(tables_dir, "GSE115440_expression_matrix_rma_probe_level.csv")
write.csv(expr_df, expr_out, row.names = FALSE)

cat(sprintf("    Saved: %s\n", expr_out))
cat(sprintf("    Dimensions: %d probes/features x %d samples\n\n", nrow(expr_mat), ncol(expr_mat)))

# ---- 10. Boxplot AFTER RMA ---------------------------------------------------
cat(">>> Step 9: Boxplot AFTER RMA...\n")

png(
  file.path(figures_dir, "GSE115440_boxplot_after_rma.png"),
  width = 1800,
  height = 1200,
  res = 200
)

par(mar = c(8, 5, 4, 2))

boxplot(
  as.data.frame(expr_mat),
  main = "GSE115440: RMA-normalized intensities",
  las = 2,
  col = sample_colors,
  ylab = "log2 intensity (RMA)"
)

legend("topright", legend = levels(group_factor), fill = group_colors, bty = "n", cex = 0.8)

dev.off()

cat("    Saved: GSE115440_boxplot_after_rma.png\n\n")

# ---- 11. Density plot AFTER RMA ----------------------------------------------
cat(">>> Step 10: Density plot AFTER RMA...\n")

density_df <- reshape2::melt(
  as.data.frame(expr_mat),
  variable.name = "sample",
  value.name = "expression"
)

density_df$group <- pData(eset)$group[match(density_df$sample, sampleNames(eset))]

p_density <- ggplot(density_df, aes(x = expression, color = sample, linetype = group)) +
  geom_density(alpha = 0.7) +
  theme_bw(base_size = 12) +
  labs(
    title = "GSE115440: density of RMA-normalized expression",
    x = "log2 intensity (RMA)",
    y = "density"
  )

ggsave(
  file.path(figures_dir, "GSE115440_density_after_rma.png"),
  p_density,
  width = 9,
  height = 6,
  dpi = 200
)

cat("    Saved: GSE115440_density_after_rma.png\n\n")

# ---- 12. PCA AFTER RMA -------------------------------------------------------
cat(">>> Step 11: PCA AFTER RMA...\n")

pca <- prcomp(t(expr_mat), scale. = FALSE)
var_explained <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)

pca_df <- data.frame(
  sample = colnames(expr_mat),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = pData(eset)$group,
  stringsAsFactors = FALSE
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "GSE115440: PCA of RMA-normalized samples",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  scale_color_brewer(palette = "Set2")

ggsave(
  file.path(figures_dir, "GSE115440_pca_after_rma.png"),
  p_pca,
  width = 8,
  height = 6,
  dpi = 200
)

cat("    Saved: GSE115440_pca_after_rma.png\n\n")

# ---- 13. Sample correlation heatmap ------------------------------------------
cat(">>> Step 12: Sample correlation heatmap...\n")

cor_mat <- cor(expr_mat, method = "pearson")

annot_col <- data.frame(group = pData(eset)$group)
rownames(annot_col) <- colnames(cor_mat)

ann_colors <- list(group = setNames(group_colors, levels(group_factor)))

pheatmap::pheatmap(
  cor_mat,
  annotation_col = annot_col,
  annotation_row = annot_col,
  annotation_colors = ann_colors,
  main = "GSE115440: sample-sample Pearson correlation after RMA",
  display_numbers = TRUE,
  number_format = "%.3f",
  fontsize_number = 7,
  filename = file.path(figures_dir, "GSE115440_sample_correlation_heatmap.png"),
  width = 9,
  height = 8
)

cat("    Saved: GSE115440_sample_correlation_heatmap.png\n\n")

# ---- 14. QC summary table ----------------------------------------------------
cat(">>> Step 13: QC summary table...\n")

qc_summary <- data.frame(
  GSM = colnames(expr_mat),
  group = pData(eset)$group,
  cel_file = cel_map$cel_base,
  n_probes = nrow(expr_mat),
  mean_expr = round(colMeans(expr_mat), 4),
  median_expr = round(apply(expr_mat, 2, median), 4),
  sd_expr = round(apply(expr_mat, 2, sd), 4),
  min_expr = round(apply(expr_mat, 2, min), 4),
  max_expr = round(apply(expr_mat, 2, max), 4),
  PC1 = round(pca_df$PC1, 4),
  PC2 = round(pca_df$PC2, 4),
  stringsAsFactors = FALSE
)

qc_out <- file.path(tables_dir, "GSE115440_qc_summary.csv")
write.csv(qc_summary, qc_out, row.names = FALSE)

cat(sprintf("    Saved: %s\n", qc_out))
print(qc_summary)
cat("\n")

# ---- 15. Wrap-up -------------------------------------------------------------
cat("==============================================================\n")
cat(" Script 02 finished successfully.\n")
cat("==============================================================\n")
cat(" Outputs:\n")
cat(sprintf("  - %s\n", expr_out))
cat(sprintf("  - %s\n", qc_out))
cat(sprintf("  - %s/GSE115440_boxplot_before_rma.png\n", figures_dir))
cat(sprintf("  - %s/GSE115440_boxplot_after_rma.png\n", figures_dir))
cat(sprintf("  - %s/GSE115440_density_after_rma.png\n", figures_dir))
cat(sprintf("  - %s/GSE115440_pca_after_rma.png\n", figures_dir))
cat(sprintf("  - %s/GSE115440_sample_correlation_heatmap.png\n", figures_dir))
cat("\n Next suggested script:\n")
cat("   03_deg_limma_surgery_vs_control.R\n")
cat("==============================================================\n")
