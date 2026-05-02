# ============================================================================
# Script: 03_deg_limma_surgery_vs_control.R
# Project: MitoMet-POCD / NeuroMitoMap
# Dataset: GSE115440 (Validation Tier 2)
# Platform: GPL11533 Affymetrix Mouse Gene 1.1 ST Array
# Purpose: Probe-level differential expression analysis with limma.
# Primary contrast: Surgery vs Control.
# Surgery+Maresin1 samples are EXCLUDED from this analysis.
#
# IMPORTANT NOTES:
# - All DEG results are PROBE-LEVEL, not gene-level.
# - Sample size is n = 3 vs 3, which is small.
# - DEG results are SECONDARY supporting evidence only.
# - Main project evidence will come from module-level mitochondrial /
#   inflammatory scoring in later scripts.
# - GSE95426 is frozen and is NOT touched here.
# ============================================================================

cat("=========================================================\n")
cat(" GSE115440 | Script 03 | limma DEG: Surgery vs Control \n")
cat(" Probe-level analysis, n = 3 vs 3, secondary evidence \n")
cat("=========================================================\n\n")

# ---- 0. User library path ----------------------------------------------------
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
  Sys.setenv(R_LIBS_USER = user_lib)
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
cat("[INFO] R library paths:\n")
print(.libPaths())
cat("\n")

# ---- 1. Package installation checks ------------------------------------------
required_cran <- c("ggplot2", "pheatmap", "RColorBrewer")
required_bioc <- c("limma")

install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[INFO] Installing CRAN package: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[INFO] Installing Bioconductor package: %s\n", pkg))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

invisible(lapply(required_cran, install_if_missing_cran))
invisible(lapply(required_bioc, install_if_missing_bioc))

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

cat("[OK] All required packages loaded.\n\n")

# ---- 2. Path configuration ---------------------------------------------------
project_root <- getwd()
cat(sprintf("[INFO] Project root: %s\n", project_root))

dir_base <- file.path("analysis", "02_validation_GSE115440")
dir_tables <- file.path(dir_base, "results", "tables")
dir_figs <- file.path(dir_base, "results", "figures")

dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)

path_expr <- file.path(
  dir_tables,
  "GSE115440_expression_matrix_rma_probe_level.csv"
)

path_meta <- file.path(
  "data",
  "metadata",
  "GSE115440_metadata_verified.csv"
)

path_full <- file.path(
  dir_tables,
  "GSE115440_limma_surgery_vs_control_probe_level.csv"
)

path_top50 <- file.path(
  dir_tables,
  "GSE115440_limma_top50_surgery_vs_control_probe_level.csv"
)

path_summary <- file.path(
  dir_tables,
  "GSE115440_deg_summary_surgery_vs_control.csv"
)

path_volcano <- file.path(
  dir_figs,
  "GSE115440_volcano_surgery_vs_control_probe_level.png"
)

path_ma <- file.path(
  dir_figs,
  "GSE115440_MAplot_surgery_vs_control_probe_level.png"
)

path_pvalhist <- file.path(
  dir_figs,
  "GSE115440_pvalue_histogram_surgery_vs_control_probe_level.png"
)

path_heatmap <- file.path(
  dir_figs,
  "GSE115440_top50_heatmap_surgery_vs_control_probe_level.png"
)

# ---- 3. Load expression matrix ----------------------------------------------
cat("[STEP 1] Loading RMA expression matrix...\n")

if (!file.exists(path_expr)) {
  stop(sprintf("Expression matrix not found at: %s", path_expr))
}

expr_df <- read.csv(
  path_expr,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_mat <- as.matrix(expr_df)
storage.mode(expr_mat) <- "numeric"

cat(sprintf(" - Probes/features: %d\n", nrow(expr_mat)))
cat(sprintf(" - Samples in matrix: %d\n", ncol(expr_mat)))
cat(sprintf(" - First sample columns: %s\n", paste(head(colnames(expr_mat), 6), collapse = ", ")))
cat("\n")

# ---- 4. Load metadata --------------------------------------------------------
cat("[STEP 2] Loading verified metadata...\n")

if (!file.exists(path_meta)) {
  stop(sprintf("Metadata file not found at: %s", path_meta))
}

meta <- read.csv(
  path_meta,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(sprintf(" - Metadata rows: %d\n", nrow(meta)))
cat(sprintf(" - Metadata columns: %s\n", paste(colnames(meta), collapse = ", ")))

gsm_col_candidates <- c(
  "sample_geo_accession",
  "geo_accession",
  "GSM",
  "gsm",
  "gsm_id",
  "GSM_ID",
  "sample_id",
  "Sample",
  "sample",
  "accession"
)

gsm_col <- NA_character_

for (cand in gsm_col_candidates) {
  if (cand %in% colnames(meta)) {
    gsm_col <- cand
    break
  }
}

if (is.na(gsm_col)) {
  is_gsm_col <- vapply(
    meta,
    function(x) any(grepl("^GSM\\d+", as.character(x))),
    logical(1)
  )
  if (any(is_gsm_col)) {
    gsm_col <- names(is_gsm_col)[which(is_gsm_col)[1]]
  } else {
    stop("Could not auto-detect a GSM column in metadata.")
  }
}

cat(sprintf(" - Using GSM column: '%s'\n", gsm_col))

group_col_candidates <- c("group", "Group", "condition", "Condition", "treatment", "Treatment")
group_col <- NA_character_

for (cand in group_col_candidates) {
  if (cand %in% colnames(meta)) {
    group_col <- cand
    break
  }
}

if (is.na(group_col)) {
  stop(
    "Could not find a group/condition column in metadata. Expected one of: ",
    paste(group_col_candidates, collapse = ", ")
  )
}

cat(sprintf(" - Using group column: '%s'\n\n", group_col))

# ---- 5. Match expression columns to metadata GSM IDs --------------------------
cat("[STEP 3] Matching expression columns to metadata GSM IDs...\n")

extract_gsm <- function(x) {
  m <- regmatches(x, regexpr("GSM\\d+", x))
  if (length(m) == 0 || nchar(m) == 0) {
    x
  } else {
    m
  }
}

expr_cols_raw <- colnames(expr_mat)
expr_cols_gsm <- vapply(expr_cols_raw, extract_gsm, character(1))
colnames(expr_mat) <- expr_cols_gsm

meta_gsm <- as.character(meta[[gsm_col]])

common_gsm <- intersect(expr_cols_gsm, meta_gsm)

cat(sprintf(" - GSM IDs in expression matrix: %d\n", length(expr_cols_gsm)))
cat(sprintf(" - GSM IDs in metadata: %d\n", length(meta_gsm)))
cat(sprintf(" - Overlap: %d\n", length(common_gsm)))

if (length(common_gsm) == 0) {
  stop("No overlap between expression matrix columns and metadata GSM IDs.")
}

# Use expression-matrix order for matched samples.
meta_sub <- meta[match(common_gsm, meta_gsm), , drop = FALSE]
expr_sub <- expr_mat[, common_gsm, drop = FALSE]

cat("\n")

# ---- 6. Keep only Control and Surgery samples --------------------------------
cat("[STEP 4] Selecting Control and Surgery samples, excluding Surgery+Maresin1...\n")

control_gsm_known <- c("GSM3178275", "GSM3178276", "GSM3178277")
surgery_gsm_known <- c("GSM3178278", "GSM3178279", "GSM3178280")
maresin_gsm_known <- c("GSM3178281", "GSM3178282", "GSM3178283")

group_raw <- as.character(meta_sub[[group_col]])

is_maresin <- grepl("maresin", group_raw, ignore.case = TRUE)
is_control_lbl <- grepl("control", group_raw, ignore.case = TRUE) & !is_maresin
is_surgery_lbl <- grepl("surgery", group_raw, ignore.case = TRUE) & !is_maresin

keep_idx <- is_control_lbl | is_surgery_lbl

if (sum(keep_idx) < 6) {
  cat("[WARN] Group-label parsing did not yield 6 samples. Falling back to verified GSM lists.\n")
  is_control_lbl <- meta_sub[[gsm_col]] %in% control_gsm_known
  is_surgery_lbl <- meta_sub[[gsm_col]] %in% surgery_gsm_known
  is_maresin <- meta_sub[[gsm_col]] %in% maresin_gsm_known
  keep_idx <- is_control_lbl | is_surgery_lbl
}

meta_kept <- meta_sub[keep_idx, , drop = FALSE]
expr_kept <- expr_sub[, keep_idx, drop = FALSE]

g_clean <- ifelse(
  meta_kept[[gsm_col]] %in% control_gsm_known,
  "Control",
  ifelse(
    meta_kept[[gsm_col]] %in% surgery_gsm_known,
    "Surgery",
    ifelse(
      grepl("control", as.character(meta_kept[[group_col]]), ignore.case = TRUE),
      "Control",
      "Surgery"
    )
  )
)

group_factor <- factor(g_clean, levels = c("Control", "Surgery"))

cat(sprintf(" - Samples retained: %d\n", ncol(expr_kept)))
cat(sprintf(
  " - Control: %d -> %s\n",
  sum(group_factor == "Control"),
  paste(meta_kept[[gsm_col]][group_factor == "Control"], collapse = ", ")
))
cat(sprintf(
  " - Surgery: %d -> %s\n",
  sum(group_factor == "Surgery"),
  paste(meta_kept[[gsm_col]][group_factor == "Surgery"], collapse = ", ")
))
cat(sprintf(" - Excluded as Surgery+Maresin1: %d\n", sum(is_maresin)))
cat("\n")

if (length(unique(group_factor)) != 2) {
  stop("Group factor must have exactly 2 levels: Control and Surgery.")
}

if (sum(group_factor == "Control") != 3 || sum(group_factor == "Surgery") != 3) {
  stop("Expected exactly 3 Control and 3 Surgery samples for primary limma contrast.")
}

# ---- 7. limma analysis -------------------------------------------------------
cat("[STEP 5] Running limma, contrast: Surgery - Control...\n")

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)
rownames(design) <- colnames(expr_kept)

cat(" Design matrix:\n")
print(design)

contrast_mat <- makeContrasts(
  SurgeryVsControl = Surgery - Control,
  levels = design
)

cat(" Contrast matrix:\n")
print(contrast_mat)

fit <- lmFit(expr_kept, design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

tt_full <- topTable(
  fit2,
  coef = "SurgeryVsControl",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

tt_full$ProbeID <- rownames(tt_full)

tt_full <- tt_full[, c(
  "ProbeID",
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B"
)]

cat(sprintf(" - Probes/features tested: %d\n\n", nrow(tt_full)))

# ---- 8. Save result tables ---------------------------------------------------
cat("[STEP 6] Writing result tables...\n")

write.csv(tt_full, path_full, row.names = FALSE)
cat(sprintf(" [OK] Full limma table -> %s\n", path_full))

tt_top50 <- head(tt_full, 50)
write.csv(tt_top50, path_top50, row.names = FALSE)
cat(sprintf(" [OK] Top-50 table -> %s\n", path_top50))

# ---- 9. DEG summary table ----------------------------------------------------
cat("[STEP 7] Building DEG summary...\n")

fc_thresh <- log2(1.5)

count_de <- function(p_col, p_thr) {
  up <- sum(
    tt_full[[p_col]] < p_thr & tt_full$logFC >= fc_thresh,
    na.rm = TRUE
  )
  down <- sum(
    tt_full[[p_col]] < p_thr & tt_full$logFC <= -fc_thresh,
    na.rm = TRUE
  )
  data.frame(up = up, down = down, total = up + down)
}

c_raw_05 <- count_de("P.Value", 0.05)
c_adj_05 <- count_de("adj.P.Val", 0.05)
c_adj_10 <- count_de("adj.P.Val", 0.10)

summary_df <- data.frame(
  contrast = "Surgery_vs_Control",
  level = "probe",
  n_control = sum(group_factor == "Control"),
  n_surgery = sum(group_factor == "Surgery"),
  n_probes_tested = nrow(tt_full),
  threshold_logFC = fc_thresh,
  raw_p_lt_0.05_up = c_raw_05$up,
  raw_p_lt_0.05_down = c_raw_05$down,
  raw_p_lt_0.05_total = c_raw_05$total,
  adj_p_lt_0.05_up = c_adj_05$up,
  adj_p_lt_0.05_down = c_adj_05$down,
  adj_p_lt_0.05_total = c_adj_05$total,
  adj_p_lt_0.10_up = c_adj_10$up,
  adj_p_lt_0.10_down = c_adj_10$down,
  adj_p_lt_0.10_total = c_adj_10$total,
  note = "n=3 vs 3; PROBE-LEVEL; secondary supporting evidence only.",
  stringsAsFactors = FALSE
)

write.csv(summary_df, path_summary, row.names = FALSE)
cat(sprintf(" [OK] DEG summary -> %s\n", path_summary))
cat("\n DEG summary preview:\n")
print(summary_df)
cat("\n")

# ---- 10. Figures -------------------------------------------------------------
cat("[STEP 8] Generating figures...\n")

# Volcano plot
volcano_df <- tt_full
volcano_df$neglog10P <- -log10(volcano_df$P.Value)
volcano_df$Direction <- "NS"
volcano_df$Direction[
  volcano_df$P.Value < 0.05 & volcano_df$logFC >= fc_thresh
] <- "Up"
volcano_df$Direction[
  volcano_df$P.Value < 0.05 & volcano_df$logFC <= -fc_thresh
] <- "Down"

volcano_df$Direction <- factor(volcano_df$Direction, levels = c("Down", "NS", "Up"))

p_volcano <- ggplot(
  volcano_df,
  aes(x = logFC, y = neglog10P, color = Direction)
) +
  geom_point(alpha = 0.6, size = 1.1) +
  scale_color_manual(values = c("Down" = "#3B82F6", "NS" = "grey70", "Up" = "#EF4444")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    title = "GSE115440: Surgery vs Control, probe-level, n=3 vs 3",
    subtitle = "Secondary evidence only; interpret with caution",
    x = "log2 fold change",
    y = "-log10(raw p-value)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(path_volcano, p_volcano, width = 7, height = 6, dpi = 300)
cat(sprintf(" [OK] Volcano plot -> %s\n", path_volcano))

# MA plot
ma_df <- tt_full
ma_df$Direction <- "NS"
ma_df$Direction[
  ma_df$P.Value < 0.05 & ma_df$logFC >= fc_thresh
] <- "Up"
ma_df$Direction[
  ma_df$P.Value < 0.05 & ma_df$logFC <= -fc_thresh
] <- "Down"

ma_df$Direction <- factor(ma_df$Direction, levels = c("Down", "NS", "Up"))

p_ma <- ggplot(ma_df, aes(x = AveExpr, y = logFC, color = Direction)) +
  geom_point(alpha = 0.6, size = 1.1) +
  scale_color_manual(values = c("Down" = "#3B82F6", "NS" = "grey70", "Up" = "#EF4444")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = c(-fc_thresh, fc_thresh), linetype = "dotted", color = "grey50") +
  labs(
    title = "GSE115440: MA plot, Surgery vs Control, probe-level",
    subtitle = "Secondary evidence only; n=3 vs 3",
    x = "Average expression, log2",
    y = "log2 fold change"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(path_ma, p_ma, width = 7, height = 6, dpi = 300)
cat(sprintf(" [OK] MA plot -> %s\n", path_ma))

# P-value histogram
p_hist <- ggplot(tt_full, aes(x = P.Value)) +
  geom_histogram(
    breaks = seq(0, 1, by = 0.05),
    fill = "#60A5FA",
    color = "white"
  ) +
  labs(
    title = "GSE115440: p-value distribution, Surgery vs Control",
    subtitle = "Probe-level, n=3 vs 3",
    x = "Raw p-value",
    y = "Number of probes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(path_pvalhist, p_hist, width = 7, height = 5, dpi = 300)
cat(sprintf(" [OK] P-value histogram -> %s\n", path_pvalhist))

# Top-50 heatmap
top50_probes <- tt_top50$ProbeID
top50_probes <- top50_probes[top50_probes %in% rownames(expr_kept)]

if (length(top50_probes) >= 2) {
  heat_mat <- expr_kept[top50_probes, , drop = FALSE]
  row_sd <- apply(heat_mat, 1, sd, na.rm = TRUE)
  heat_mat <- heat_mat[!is.na(row_sd) & row_sd > 0, , drop = FALSE]

  if (nrow(heat_mat) >= 2) {
    heat_scaled <- t(scale(t(heat_mat)))

    annot_col <- data.frame(Group = group_factor)
    rownames(annot_col) <- colnames(heat_scaled)

    annot_colors <- list(
      Group = c(Control = "#10B981", Surgery = "#EF4444")
    )

    png(path_heatmap, width = 1800, height = 2400, res = 220)
    pheatmap(
      heat_scaled,
      annotation_col = annot_col,
      annotation_colors = annot_colors,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 7,
      fontsize_col = 9,
      main = "GSE115440 top-50 probes: Surgery vs Control, probe-level"
    )
    dev.off()

    cat(sprintf(" [OK] Top-50 heatmap -> %s\n", path_heatmap))
  } else {
    cat(" [WARN] After variance filtering, fewer than 2 probes remained; heatmap skipped.\n")
  }
} else {
  cat(" [WARN] Could not draw top-50 heatmap; insufficient probes after intersection.\n")
}

cat("\n")

# ---- 11. Final notes ---------------------------------------------------------
cat("=========================================================\n")
cat(" Script 03 finished.\n")
cat(" REMINDERS:\n")
cat(" - All DEG outputs are PROBE-LEVEL, not gene-level.\n")
cat(" - Sample size n = 3 vs 3: DEG findings are SECONDARY\n")
cat("   supporting evidence only.\n")
cat(" - Do NOT claim strong DEG discovery from this step.\n")
cat(" - Main project evidence will come from module-level\n")
cat("   mitochondrial / inflammatory scoring downstream.\n")
cat(" - GSE95426 was NOT touched.\n")
cat("\n")
cat(" NEXT SUGGESTED SCRIPT:\n")
cat(" 04_module_scoring_GSVA_ssGSEA.R\n")
cat("=========================================================\n")
