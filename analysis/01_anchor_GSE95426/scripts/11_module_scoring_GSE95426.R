#!/usr/bin/env Rscript

# Script 11: GSE95426 module scoring and POCD vs Control contrast
#
# Pre-mortem:
# 1. Most likely failure mode: sample labels or group order are wrong.
# 2. How to detect it: assert manifest has exactly 6 Control and 6 POCD samples, all present in expression matrix.
# 3. First thing to inspect if wrong: data/metadata/GSE95426_manifest_frozen.csv and expression matrix column names.
#
# Decision rules declared before seeing results:
# abs(Cohen's d) >= 0.8: strong
# 0.5 <= abs(Cohen's d) < 0.8: moderate
# 0.2 <= abs(Cohen's d) < 0.5: weak
# abs(Cohen's d) < 0.2: null_effect
#
# Scope:
# - GSE95426 only.
# - Predefined module-level scoring only.
# - No cross-dataset raw expression merging.
# - p-values are descriptive only.
# - DEG/single-gene claims are not made here.
# - PASS_WITH_CAVEAT context: GPL22782 annotation rescue supports predefined module-level scoring.

options(stringsAsFactors = FALSE)

log_dir <- "analysis/01_anchor_GSE95426/logs"
table_dir <- "analysis/01_anchor_GSE95426/results/tables"
figure_dir <- "analysis/01_anchor_GSE95426/results/figures"

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "GSE95426_script11_module_scoring.log")
sink(log_file, split = TRUE)
on.exit({
  cat("\n[sessionInfo]\n")
  print(sessionInfo())
  sink()
}, add = TRUE)

message_ts <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required file missing: %s", path), call. = FALSE)
  }
}

cohens_d <- function(x_case, x_ctrl) {
  x_case <- as.numeric(x_case)
  x_ctrl <- as.numeric(x_ctrl)
  n_case <- length(x_case)
  n_ctrl <- length(x_ctrl)
  sd_case <- stats::sd(x_case)
  sd_ctrl <- stats::sd(x_ctrl)

  if (n_case < 2 || n_ctrl < 2) {
    return(NA_real_)
  }

  pooled_var <- ((n_case - 1) * sd_case^2 + (n_ctrl - 1) * sd_ctrl^2) /
    (n_case + n_ctrl - 2)

  if (!is.finite(pooled_var) || pooled_var <= 0) {
    return(NA_real_)
  }

  (mean(x_case) - mean(x_ctrl)) / sqrt(pooled_var)
}

effect_size_label <- function(d) {
  ad <- abs(d)
  if (!is.finite(ad)) return("not_evaluable")
  if (ad >= 0.8) return("strong")
  if (ad >= 0.5) return("moderate")
  if (ad >= 0.2) return("weak")
  "null_effect"
}

direction_label <- function(delta, d) {
  if (!is.finite(delta) || !is.finite(d)) return("not_evaluable")
  if (abs(d) < 0.2) return("near_null")
  if (delta > 0) return("up_in_pocd")
  if (delta < 0) return("down_in_pocd")
  "near_null"
}

safe_t_test_p <- function(x_case, x_ctrl) {
  out <- tryCatch({
    if (length(unique(c(x_case, x_ctrl))) < 2) return(NA_real_)
    stats::t.test(x_case, x_ctrl)$p.value
  }, error = function(e) NA_real_)
  as.numeric(out)
}

message_ts("Checking input paths")

expr_path <- "data/processed/GSE95426_expression_matrix_gene_level.csv"
module_path <- "data/modules/stage_B5_target_modules.csv"
manifest_path <- "data/metadata/GSE95426_manifest_frozen.csv"

stop_if_missing(expr_path)
stop_if_missing(module_path)
stop_if_missing(manifest_path)

message_ts("Reading inputs")

expr_df <- read.csv(expr_path, check.names = FALSE)
modules_df <- read.csv(module_path, check.names = FALSE)
manifest_df <- read.csv(manifest_path, check.names = FALSE)

required_expr_cols <- c("gene_name", "selected_probe_id")
required_module_cols <- c("module", "gene_symbol", "module_class")
required_manifest_cols <- c("sample_id", "group", "include_primary_analysis")

missing_expr <- setdiff(required_expr_cols, colnames(expr_df))
missing_module <- setdiff(required_module_cols, colnames(modules_df))
missing_manifest <- setdiff(required_manifest_cols, colnames(manifest_df))

if (length(missing_expr) > 0) stop("Expression matrix missing columns: ", paste(missing_expr, collapse = ", "))
if (length(missing_module) > 0) stop("Module file missing columns: ", paste(missing_module, collapse = ", "))
if (length(missing_manifest) > 0) stop("Manifest missing columns: ", paste(missing_manifest, collapse = ", "))

message_ts("Validating frozen manifest")

manifest_use <- manifest_df[manifest_df$include_primary_analysis == "yes", , drop = FALSE]
manifest_use$group <- as.character(manifest_use$group)
manifest_use$sample_id <- as.character(manifest_use$sample_id)

control_samples <- manifest_use$sample_id[manifest_use$group == "Control"]
pocd_samples <- manifest_use$sample_id[manifest_use$group == "POCD"]

if (length(control_samples) != 6) {
  stop(sprintf("Expected 6 Control samples, found %d", length(control_samples)))
}
if (length(pocd_samples) != 6) {
  stop(sprintf("Expected 6 POCD samples, found %d", length(pocd_samples)))
}

sample_ids <- c(control_samples, pocd_samples)
missing_samples <- setdiff(sample_ids, colnames(expr_df))
if (length(missing_samples) > 0) {
  stop("Manifest samples missing from expression matrix: ", paste(missing_samples, collapse = ", "))
}

message_ts("Preparing gene-level expression matrix")

if (anyDuplicated(expr_df$gene_name) > 0) {
  stop("gene_name contains duplicated rows; expected already-collapsed gene-level matrix.")
}

expr_mat <- as.matrix(expr_df[, sample_ids, drop = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr_df$gene_name

if (any(!is.finite(expr_mat))) {
  stop("Expression matrix contains non-finite values.")
}

message_ts(sprintf("Expression matrix genes: %d", nrow(expr_mat)))
message_ts(sprintf("Expression matrix samples: %d", ncol(expr_mat)))
message_ts(sprintf("Control samples: %s", paste(control_samples, collapse = ", ")))
message_ts(sprintf("POCD samples: %s", paste(pocd_samples, collapse = ", ")))

message_ts("Computing per-gene z-scores across samples")

gene_means <- rowMeans(expr_mat)
gene_sds <- apply(expr_mat, 1, stats::sd)

z_mat <- sweep(expr_mat, 1, gene_means, "-")
nonzero_sd <- is.finite(gene_sds) & gene_sds > 0
z_mat[nonzero_sd, ] <- sweep(z_mat[nonzero_sd, , drop = FALSE], 1, gene_sds[nonzero_sd], "/")
z_mat[!nonzero_sd, ] <- 0

message_ts(sprintf("Zero-variance genes assigned z-score 0: %d", sum(!nonzero_sd)))

message_ts("Computing module scores")

modules_df$module <- as.character(modules_df$module)
modules_df$gene_symbol <- as.character(modules_df$gene_symbol)
modules_df$module_class <- as.character(modules_df$module_class)

# Add a reproducible negative-control module for Script 11 group-separation QC.
# This is not a biological module and must not be interpreted mechanistically.
# It is generated after input validation from the GSE95426 gene-level universe with a fixed seed.
set.seed(95426)
target_genes_all <- unique(modules_df$gene_symbol)
negative_pool <- setdiff(rownames(z_mat), target_genes_all)
if (length(negative_pool) < 200) {
  stop(sprintf("Negative-control pool too small: %d genes", length(negative_pool)))
}
negative_genes <- sample(negative_pool, size = 200, replace = FALSE)
negative_df <- data.frame(
  module = "RANDOM_CONTROL_200",
  gene_symbol = negative_genes,
  module_class = "negative_control",
  stringsAsFactors = FALSE
)
modules_df <- rbind(modules_df, negative_df)

module_names <- unique(modules_df$module)

score_list <- list()
overlap_list <- list()
contrast_list <- list()
summary_list <- list()

for (module_name in module_names) {
  module_rows <- modules_df[modules_df$module == module_name, , drop = FALSE]
  module_class <- unique(module_rows$module_class)
  if (length(module_class) != 1) {
    module_class <- paste(module_class, collapse = ";")
  }

  module_genes <- unique(module_rows$gene_symbol)
  detected_genes <- intersect(module_genes, rownames(z_mat))
  missing_genes <- setdiff(module_genes, rownames(z_mat))

  total_genes <- length(module_genes)
  detected_n <- length(detected_genes)
  coverage_pct <- ifelse(total_genes > 0, 100 * detected_n / total_genes, NA_real_)

  if (detected_n == 0) {
    module_scores <- rep(NA_real_, length(sample_ids))
    names(module_scores) <- sample_ids
  } else {
    module_scores <- colMeans(z_mat[detected_genes, , drop = FALSE], na.rm = TRUE)
  }

  score_df <- data.frame(
    module = module_name,
    module_class = module_class,
    sample_id = sample_ids,
    group = manifest_use$group[match(sample_ids, manifest_use$sample_id)],
    module_score_mean_zscore = as.numeric(module_scores[sample_ids]),
    detected_genes = detected_n,
    total_module_genes = total_genes,
    coverage_pct = coverage_pct,
    stringsAsFactors = FALSE
  )
  score_list[[module_name]] <- score_df

  overlap_list[[module_name]] <- data.frame(
    module = module_name,
    module_class = module_class,
    total_module_genes = total_genes,
    detected_genes = detected_n,
    coverage_pct = coverage_pct,
    detected_gene_symbols = paste(detected_genes, collapse = ";"),
    missing_gene_symbols = paste(missing_genes, collapse = ";"),
    stringsAsFactors = FALSE
  )

  ctrl_values <- score_df$module_score_mean_zscore[score_df$group == "Control"]
  pocd_values <- score_df$module_score_mean_zscore[score_df$group == "POCD"]

  mean_control <- mean(ctrl_values, na.rm = TRUE)
  mean_pocd <- mean(pocd_values, na.rm = TRUE)
  delta <- mean_pocd - mean_control
  d <- cohens_d(pocd_values, ctrl_values)
  p_val <- safe_t_test_p(pocd_values, ctrl_values)

  contrast_list[[module_name]] <- data.frame(
    module = module_name,
    module_class = module_class,
    n_control = length(ctrl_values),
    n_pocd_or_surgery = length(pocd_values),
    mean_control = mean_control,
    mean_pocd_or_surgery = mean_pocd,
    delta_pocd_vs_ctrl = delta,
    cohens_d = d,
    effect_size_category = effect_size_label(d),
    direction = direction_label(delta, d),
    descriptive_ttest_p = p_val,
    detected_genes = detected_n,
    total_module_genes = total_genes,
    coverage_pct = coverage_pct,
    caveat = "GSE95426_PASS_WITH_CAVEAT_module_level_only_p_values_descriptive",
    stringsAsFactors = FALSE
  )

  for (grp in c("Control", "POCD")) {
    vals <- score_df$module_score_mean_zscore[score_df$group == grp]
    summary_list[[paste(module_name, grp, sep = "__")]] <- data.frame(
      module = module_name,
      module_class = module_class,
      group = grp,
      n = length(vals),
      mean_module_score = mean(vals, na.rm = TRUE),
      sd_module_score = stats::sd(vals, na.rm = TRUE),
      median_module_score = stats::median(vals, na.rm = TRUE),
      min_module_score = min(vals, na.rm = TRUE),
      max_module_score = max(vals, na.rm = TRUE),
      detected_genes = detected_n,
      total_module_genes = total_genes,
      coverage_pct = coverage_pct,
      stringsAsFactors = FALSE
    )
  }
}

scores_all <- do.call(rbind, score_list)
overlap_all <- do.call(rbind, overlap_list)
contrast_all <- do.call(rbind, contrast_list)
summary_all <- do.call(rbind, summary_list)

directionality_summary <- contrast_all[, c(
  "module",
  "module_class",
  "delta_pocd_vs_ctrl",
  "cohens_d",
  "effect_size_category",
  "direction",
  "detected_genes",
  "total_module_genes",
  "coverage_pct",
  "caveat"
)]

message_ts("Writing tables")

write.csv(
  scores_all,
  file.path(table_dir, "GSE95426_module_scores_mean_zscore.csv"),
  row.names = FALSE
)

write.csv(
  summary_all,
  file.path(table_dir, "GSE95426_module_group_summary.csv"),
  row.names = FALSE
)

write.csv(
  contrast_all,
  file.path(table_dir, "GSE95426_module_pocd_vs_control_results.csv"),
  row.names = FALSE
)

write.csv(
  directionality_summary,
  file.path(table_dir, "GSE95426_module_directionality_summary.csv"),
  row.names = FALSE
)

write.csv(
  overlap_all,
  file.path(table_dir, "GSE95426_module_gene_overlap.csv"),
  row.names = FALSE
)

message_ts("Computing random-control sensitivity QC")

set.seed(95426)
random_sensitivity_list <- list()

for (i in seq_len(50)) {
  random_genes_i <- sample(negative_pool, size = 200, replace = FALSE)
  random_scores_i <- colMeans(z_mat[random_genes_i, sample_ids, drop = FALSE], na.rm = TRUE)

  ctrl_i <- random_scores_i[control_samples]
  pocd_i <- random_scores_i[pocd_samples]

  delta_i <- mean(pocd_i, na.rm = TRUE) - mean(ctrl_i, na.rm = TRUE)
  d_i <- cohens_d(pocd_i, ctrl_i)

  random_sensitivity_list[[i]] <- data.frame(
    random_module = sprintf("RANDOM_CONTROL_200_%02d", i),
    delta_pocd_vs_ctrl = delta_i,
    cohens_d = d_i,
    abs_cohens_d = abs(d_i),
    effect_size_category = effect_size_label(d_i),
    direction = direction_label(delta_i, d_i),
    stringsAsFactors = FALSE
  )
}

random_sensitivity_all <- do.call(rbind, random_sensitivity_list)

write.csv(
  random_sensitivity_all,
  file.path(table_dir, "GSE95426_random_control_sensitivity_QC.csv"),
  row.names = FALSE
)

message_ts(sprintf(
  "Random-control sensitivity: abs(d) >= 0.5 in %d/50; abs(d) >= 0.8 in %d/50",
  sum(random_sensitivity_all$abs_cohens_d >= 0.5, na.rm = TRUE),
  sum(random_sensitivity_all$abs_cohens_d >= 0.8, na.rm = TRUE)
))

message_ts("Writing figures with base R")

png(file.path(figure_dir, "GSE95426_module_delta_barplot.png"), width = 1200, height = 800, res = 120)
op <- par(mar = c(12, 5, 4, 2))
barplot(
  contrast_all$delta_pocd_vs_ctrl,
  names.arg = contrast_all$module,
  las = 2,
  ylab = "Delta module score: POCD - Control",
  main = "GSE95426 module-score contrast"
)
abline(h = 0, lty = 2)
par(op)
dev.off()

png(file.path(figure_dir, "GSE95426_module_score_boxplots.png"), width = 1400, height = 900, res = 120)
op <- par(mar = c(12, 5, 4, 2))
boxplot(
  module_score_mean_zscore ~ module + group,
  data = scores_all,
  las = 2,
  ylab = "Module score: mean gene z-score",
  main = "GSE95426 per-sample module scores"
)
abline(h = 0, lty = 2)
par(op)
dev.off()

score_wide <- reshape(
  scores_all[, c("module", "sample_id", "module_score_mean_zscore")],
  idvar = "module",
  timevar = "sample_id",
  direction = "wide"
)
rownames(score_wide) <- score_wide$module
heat_mat <- as.matrix(score_wide[, setdiff(colnames(score_wide), "module"), drop = FALSE])
colnames(heat_mat) <- sub("^module_score_mean_zscore\\.", "", colnames(heat_mat))

png(file.path(figure_dir, "GSE95426_module_score_heatmap.png"), width = 1200, height = 700, res = 120)
op <- par(mar = c(8, 14, 4, 2))
image(
  x = seq_len(ncol(heat_mat)),
  y = seq_len(nrow(heat_mat)),
  z = t(heat_mat[nrow(heat_mat):1, , drop = FALSE]),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "GSE95426 module scores"
)
axis(1, at = seq_len(ncol(heat_mat)), labels = colnames(heat_mat), las = 2)
axis(2, at = seq_len(nrow(heat_mat)), labels = rev(rownames(heat_mat)), las = 2)
box()
par(op)
dev.off()

message_ts("Completed Script 11")
message_ts("Key output: analysis/01_anchor_GSE95426/results/tables/GSE95426_module_pocd_vs_control_results.csv")
message_ts("Interpretation remains conservative: module-level only; p-values descriptive only.")
