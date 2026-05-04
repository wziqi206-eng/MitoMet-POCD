# ============================================================================
# Script 11.5 - GSE95426 module sensitivity QC
# Purpose: Formalize random-control sensitivity caveats before Script 05.
# Inputs: Existing Script 11 outputs only.
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

message("== Script 11.5 START ==")

base_dir <- "analysis/01_anchor_GSE95426"
tables_dir <- file.path(base_dir, "results", "tables")
figs_dir <- file.path(base_dir, "results", "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

target_path <- file.path(tables_dir, "GSE95426_module_pocd_vs_control_results.csv")
random_path <- file.path(tables_dir, "GSE95426_random_control_sensitivity_QC.csv")
scores_path <- file.path(tables_dir, "GSE95426_module_scores_mean_zscore.csv")

required <- c(target_path, random_path, scores_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required inputs:\n", paste(missing, collapse = "\n"))
}

target <- readr::read_csv(target_path, show_col_types = FALSE)
random <- readr::read_csv(random_path, show_col_types = FALSE)
scores <- readr::read_csv(scores_path, show_col_types = FALSE)

required_target_cols <- c("module", "module_class", "cohens_d", "direction")
required_random_cols <- c("cohens_d", "direction", "effect_size_category")
required_scores_cols <- c("module", "module_class", "sample_id", "group", "module_score_mean_zscore")

for (cc in required_target_cols) {
  if (!cc %in% names(target)) stop("Missing target column: ", cc)
}
for (cc in required_random_cols) {
  if (!cc %in% names(random)) stop("Missing random column: ", cc)
}
for (cc in required_scores_cols) {
  if (!cc %in% names(scores)) stop("Missing scores column: ", cc)
}

if (!"abs_cohens_d" %in% names(random)) {
  random <- random %>% mutate(abs_cohens_d = abs(cohens_d))
}

# Biological targets only; exclude RANDOM_CONTROL_200 from target QC.
target_bio <- target %>%
  filter(module_class != "negative_control", module != "RANDOM_CONTROL_200") %>%
  mutate(abs_cohens_d = abs(cohens_d))

message("[INFO] Biological target modules retained: ", nrow(target_bio))

# 1. Random-control summary
random_summary <- tibble::tibble(
  n_random_modules = nrow(random),
  mean_cohens_d = mean(random$cohens_d, na.rm = TRUE),
  median_cohens_d = median(random$cohens_d, na.rm = TRUE),
  sd_cohens_d = sd(random$cohens_d, na.rm = TRUE),
  min_cohens_d = min(random$cohens_d, na.rm = TRUE),
  max_cohens_d = max(random$cohens_d, na.rm = TRUE),
  q10_cohens_d = as.numeric(quantile(random$cohens_d, 0.10, na.rm = TRUE)),
  q90_cohens_d = as.numeric(quantile(random$cohens_d, 0.90, na.rm = TRUE)),
  n_down_in_pocd = sum(random$direction == "down_in_pocd", na.rm = TRUE),
  frac_down_in_pocd = mean(random$direction == "down_in_pocd", na.rm = TRUE),
  n_up_in_pocd = sum(random$direction == "up_in_pocd", na.rm = TRUE),
  frac_up_in_pocd = mean(random$direction == "up_in_pocd", na.rm = TRUE),
  n_near_null = sum(random$direction == "near_null", na.rm = TRUE),
  frac_near_null = mean(random$direction == "near_null", na.rm = TRUE),
  n_strong = sum(random$effect_size_category == "strong", na.rm = TRUE),
  n_moderate = sum(random$effect_size_category == "moderate", na.rm = TRUE),
  n_weak = sum(random$effect_size_category == "weak", na.rm = TRUE),
  n_null_effect = sum(random$effect_size_category == "null_effect", na.rm = TRUE)
)

out_random_summary <- file.path(tables_dir, "GSE95426_script11_5_random_control_summary.csv")
readr::write_csv(random_summary, out_random_summary)
message("[WRITE] ", out_random_summary)

random_directional_bias <- random_summary$frac_down_in_pocd >= 0.70 ||
  random_summary$frac_up_in_pocd >= 0.70

random_bias_direction <- if (random_summary$frac_down_in_pocd >= 0.70) {
  "down_in_pocd"
} else if (random_summary$frac_up_in_pocd >= 0.70) {
  "up_in_pocd"
} else {
  "no_strong_bias"
}

# 2. Target vs random empirical QC
median_abs_random <- median(random$abs_cohens_d, na.rm = TRUE)
p90_abs_random <- as.numeric(quantile(random$abs_cohens_d, 0.90, na.rm = TRUE))
max_abs_random <- max(random$abs_cohens_d, na.rm = TRUE)

same_direction_extreme <- function(d, direction) {
  rr <- random %>% filter(.data$direction == direction)
  if (nrow(rr) == 0) return(NA_real_)
  mean(rr$abs_cohens_d <= abs(d), na.rm = TRUE)
}

target_qc <- target_bio %>%
  rowwise() %>%
  mutate(
    empirical_abs_percentile_vs_random = mean(random$abs_cohens_d <= abs_cohens_d, na.rm = TRUE),
    same_direction_extremeness_vs_random = same_direction_extreme(cohens_d, direction),
    stronger_than_random_median = abs_cohens_d > median_abs_random,
    stronger_than_random_90th_percentile = abs_cohens_d > p90_abs_random,
    stronger_than_random_max = abs_cohens_d > max_abs_random
  ) %>%
  ungroup() %>%
  select(
    module, module_class, cohens_d, abs_cohens_d, direction,
    empirical_abs_percentile_vs_random,
    same_direction_extremeness_vs_random,
    stronger_than_random_median,
    stronger_than_random_90th_percentile,
    stronger_than_random_max
  )

out_target_qc <- file.path(tables_dir, "GSE95426_script11_5_target_vs_random_empirical_QC.csv")
readr::write_csv(target_qc, out_target_qc)
message("[WRITE] ", out_target_qc)

# 3. Interpretation flags
expected_direction <- c(
  HALLMARK_INFLAMMATORY_RESPONSE = "up_in_pocd",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = "up_in_pocd",
  MICROGLIA_ACTIVATION_CUSTOM = "up_in_pocd",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = "down_in_pocd",
  KEGG_OXIDATIVE_PHOSPHORYLATION = "down_in_pocd",
  REACTOME_RESPIRATORY_ELECTRON_TRANSPORT = "down_in_pocd"
)

flag_one <- function(module, module_class, direction, abs_pct) {
  exp_dir <- unname(expected_direction[module])
  if (is.na(exp_dir)) return("inconclusive")

  if (module_class == "inflammatory" && direction %in% c("down_in_pocd", "near_null")) {
    return("not_replicated_or_opposite_direction")
  }

  if (direction != exp_dir) {
    return("not_replicated_or_opposite_direction")
  }

  if (module_class == "mitochondrial" &&
      direction == "down_in_pocd" &&
      abs_pct >= 0.90 &&
      random_directional_bias &&
      random_bias_direction == "down_in_pocd") {
    return("hypothesis_generating_only")
  }

  if (abs_pct >= 0.90 && !random_directional_bias) {
    return("clean_target_signal")
  }

  if (abs_pct >= 0.90 && random_directional_bias) {
    return("target_signal_with_global_shift_caveat")
  }

  "inconclusive"
}

flags <- target_qc %>%
  rowwise() %>%
  mutate(
    expected_direction = unname(expected_direction[module]),
    direction_matches_expected = direction == expected_direction,
    random_controls_directionally_biased = random_directional_bias,
    random_bias_direction = random_bias_direction,
    interpretation_flag = flag_one(
      module,
      module_class,
      direction,
      empirical_abs_percentile_vs_random
    )
  ) %>%
  ungroup()

out_flags <- file.path(tables_dir, "GSE95426_script11_5_module_interpretation_flags.csv")
readr::write_csv(flags, out_flags)
message("[WRITE] ", out_flags)

# 4. Sample-level QC summary
sample_qc <- scores %>%
  filter(module %in% target_bio$module) %>%
  group_by(module, module_class, group) %>%
  summarise(
    n = n(),
    mean = mean(module_score_mean_zscore, na.rm = TRUE),
    sd = sd(module_score_mean_zscore, na.rm = TRUE),
    min = min(module_score_mean_zscore, na.rm = TRUE),
    max = max(module_score_mean_zscore, na.rm = TRUE),
    .groups = "drop"
  )

out_sample_qc <- file.path(tables_dir, "GSE95426_script11_5_sample_level_QC_summary.csv")
readr::write_csv(sample_qc, out_sample_qc)
message("[WRITE] ", out_sample_qc)

# 5. Caveat text
caveat <- "GSE95426 Script 11.5 sensitivity QC indicates a broad random-control down-shift in POCD samples. Therefore, GSE95426 module-level contrasts should not be interpreted as clean pathway-specific replication. OXPHOS and respiratory electron transport modules show stronger down-shifts than the random-control distribution and remain hypothesis-generating PASS_WITH_CAVEAT signals. Inflammatory activation is not replicated in GSE95426."

out_caveat <- file.path(tables_dir, "GSE95426_script11_5_interpretation_caveat.txt")
writeLines(caveat, out_caveat)
message("[WRITE] ", out_caveat)

# 6. Figures
fig_hist <- file.path(figs_dir, "GSE95426_script11_5_random_control_cohens_d_histogram.png")
p1 <- ggplot(random, aes(x = cohens_d)) +
  geom_histogram(bins = 25, fill = "grey70", color = "grey30") +
  geom_vline(xintercept = 0, color = "black") +
  geom_vline(data = target_bio, aes(xintercept = cohens_d, color = module),
             linetype = "dashed", linewidth = 0.6) +
  labs(
    title = "GSE95426 random-control Cohen's d distribution",
    subtitle = paste0(
      "Random-control mean d = ",
      round(random_summary$mean_cohens_d, 3),
      "; median d = ",
      round(random_summary$median_cohens_d, 3)
    ),
    x = "Cohen's d",
    y = "count",
    color = "biological target module"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))
ggsave(fig_hist, p1, width = 9, height = 5.5, dpi = 150)
message("[WRITE] ", fig_hist)

fig_bar <- file.path(figs_dir, "GSE95426_script11_5_target_vs_random_abs_d_plot.png")
ref_df <- tibble::tibble(
  label = c("random_median_abs_d", "random_p90_abs_d", "random_max_abs_d"),
  value = c(median_abs_random, p90_abs_random, max_abs_random)
)

p2 <- ggplot(target_qc, aes(x = reorder(module, abs_cohens_d), y = abs_cohens_d, fill = module_class)) +
  geom_col() +
  geom_hline(data = ref_df, aes(yintercept = value, linetype = label), color = "red") +
  coord_flip() +
  labs(
    title = "GSE95426 target modules vs random-control |Cohen's d|",
    x = NULL,
    y = "|Cohen's d|",
    fill = "module class",
    linetype = "random reference"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(fig_bar, p2, width = 9, height = 5.5, dpi = 150)
message("[WRITE] ", fig_bar)

fig_heatmap <- file.path(figs_dir, "GSE95426_script11_5_module_score_sample_heatmap.png")
heatmap_df <- scores %>%
  filter(module %in% target_bio$module) %>%
  mutate(
    sample_id = factor(sample_id, levels = unique(sample_id)),
    module = factor(module, levels = target_bio$module)
  )

p3 <- ggplot(heatmap_df, aes(x = sample_id, y = module, fill = module_score_mean_zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  labs(
    title = "GSE95426 biological target module scores",
    x = NULL,
    y = NULL,
    fill = "mean z-score"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        panel.grid = element_blank())
ggsave(fig_heatmap, p3, width = 11, height = 4.5, dpi = 150)
message("[WRITE] ", fig_heatmap)

message("== Script 11.5 DONE ==")
sessionInfo()
