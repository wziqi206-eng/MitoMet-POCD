# ============================================================
# Script 05: Conservative cross-dataset module concordance
# Project: MitoMet-POCD / NeuroMitoMap
# Scope: Partial-concordance, NOT clean replication.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

options(readr.show_col_types = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

# ---------- Paths ----------
in_gse95426_contrast <- "analysis/01_anchor_GSE95426/results/tables/GSE95426_module_pocd_vs_control_results.csv"
in_gse95426_flags <- "analysis/01_anchor_GSE95426/results/tables/GSE95426_script11_5_module_interpretation_flags.csv"
in_gse115440_contrast <- "analysis/02_validation_GSE115440/results/tables/GSE115440_module_surgery_vs_control_results.csv"

out_dir_tab <- "analysis/03_cross_dataset_concordance/results/tables"
out_dir_fig <- "analysis/03_cross_dataset_concordance/results/figures"

dir.create(out_dir_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)

out_concordance_csv <- file.path(out_dir_tab, "cross_dataset_module_concordance_table.csv")
out_directionality_csv <- file.path(out_dir_tab, "cross_dataset_directionality_summary.csv")
out_flags_csv <- file.path(out_dir_tab, "cross_dataset_interpretation_flags.csv")
out_audit_csv <- file.path(out_dir_tab, "cross_dataset_input_audit.csv")
out_summary_txt <- file.path(out_dir_tab, "cross_dataset_conservative_summary.txt")

out_heatmap_png <- file.path(out_dir_fig, "cross_dataset_module_effect_size_heatmap.png")
out_tileplot_png <- file.path(out_dir_fig, "cross_dataset_directionality_tileplot.png")

# ---------- Strict input checks ----------
required_inputs <- c(in_gse95426_contrast, in_gse95426_flags, in_gse115440_contrast)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0) {
  stop("Missing required input(s):\n  - ", paste(missing_inputs, collapse = "\n  - "))
}

log_msg("All required inputs found.")

# ---------- Helpers ----------
resolve_col <- function(df, candidates, required = TRUE, label = "") {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) {
      stop(sprintf(
        "Required column not found for [%s]. Tried: %s. Available: %s",
        label,
        paste(candidates, collapse = ", "),
        paste(colnames(df), collapse = ", ")
      ))
    }
    return(NA_character_)
  }
  hit[1]
}

harmonize_direction <- function(x) {
  x_l <- tolower(as.character(x))
  dplyr::case_when(
    x_l %in% c("up_in_pocd", "up_in_surgery", "up", "upregulated", "up_regulated", "increased") ~ "up",
    x_l %in% c("down_in_pocd", "down_in_surgery", "down", "downregulated", "down_regulated", "decreased") ~ "down",
    x_l %in% c("near_null", "null", "ns", "near-null", "nearnull", "none") ~ "near_null",
    is.na(x_l) ~ NA_character_,
    TRUE ~ "near_null"
  )
}

is_random_or_negctrl <- function(mod, cls) {
  mod_l <- as.character(mod)
  cls_l <- tolower(as.character(cls))

  mod_l[is.na(mod_l)] <- ""
  cls_l[is.na(cls_l)] <- ""

  cls_l == "negative_control" |
    mod_l == "RANDOM_CONTROL_200" |
    grepl("^RANDOM_CONTROL", mod_l, ignore.case = TRUE)
}

classify_module_family <- function(cls, mod) {
  cls_l <- tolower(as.character(cls))
  mod_l <- tolower(as.character(mod))

  if (is.na(cls_l)) cls_l <- ""
  if (is.na(mod_l)) mod_l <- ""

  if (grepl("inflam|immune|cytokine|nfkb", cls_l, ignore.case = TRUE)) return("inflammatory")
  if (grepl("mito|oxphos|respir|electron|ret|tca|ox_phos", cls_l, ignore.case = TRUE)) return("mitochondrial")
  if (grepl("inflam|immune|cytokine|nfkb", mod_l, ignore.case = TRUE)) return("inflammatory")
  if (grepl("mito|oxphos|respir|electron|ret|tca", mod_l, ignore.case = TRUE)) return("mitochondrial")

  "other"
}

classify_concordance <- function(d_g95, d_g115) {
  if (is.na(d_g95) || is.na(d_g115)) return("partial_or_near_null")
  if (d_g95 == "near_null" || d_g115 == "near_null") return("partial_or_near_null")
  if (d_g95 == d_g115) return("concordant_same_direction")
  "discordant_opposite_direction"
}

final_interp <- function(dir_class, flag_g95) {
  fl <- tolower(as.character(flag_g95))

  if (!is.na(fl) && fl == "hypothesis_generating_only") {
    return("caveated_hypothesis_generating")
  }

  if (!is.na(fl) && fl == "not_replicated_or_opposite_direction") {
    if (dir_class == "concordant_same_direction") return("partial_or_near_null")
    return(dir_class)
  }

  dir_class
}

# ---------- Load ----------
log_msg("Loading GSE95426 contrast table.")
g95 <- readr::read_csv(in_gse95426_contrast)

log_msg("Loading GSE95426 Script 11.5 flag table.")
g95_flags <- readr::read_csv(in_gse95426_flags)

log_msg("Loading GSE115440 contrast table.")
g115 <- readr::read_csv(in_gse115440_contrast)

# ---------- Input audit ----------
audit_df <- bind_rows(
  tibble(file = basename(in_gse95426_contrast), column = colnames(g95)),
  tibble(file = basename(in_gse95426_flags), column = colnames(g95_flags)),
  tibble(file = basename(in_gse115440_contrast), column = colnames(g115))
)

readr::write_csv(audit_df, out_audit_csv)
log_msg("Wrote input audit: ", out_audit_csv)

# ---------- Resolve GSE95426 contrast columns ----------
c_g95_module <- resolve_col(g95, c("module"), TRUE, "GSE95426 module")
c_g95_class <- resolve_col(g95, c("module_class", "class"), TRUE, "GSE95426 module_class")
c_g95_delta <- resolve_col(g95, c("delta_pocd_vs_ctrl", "delta", "mean_diff"), TRUE, "GSE95426 delta")
c_g95_d <- resolve_col(g95, c("cohens_d", "cohen_d", "d"), TRUE, "GSE95426 Cohen's d")
c_g95_dir <- resolve_col(g95, c("direction"), TRUE, "GSE95426 direction")
c_g95_cov <- resolve_col(g95, c("coverage_pct", "coverage", "detected_pct"), FALSE, "GSE95426 coverage")
c_g95_caveat <- resolve_col(g95, c("caveat", "note", "caveat_note"), FALSE, "GSE95426 caveat")

# ---------- Resolve GSE95426 Script 11.5 flag columns ----------
c_f_module <- resolve_col(g95_flags, c("module"), TRUE, "GSE95426 flags module")
c_f_interp <- resolve_col(g95_flags, c("interpretation_flag"), TRUE, "GSE95426 flags interpretation_flag")
c_f_expected <- resolve_col(g95_flags, c("expected_direction"), FALSE, "GSE95426 expected_direction")
c_f_match <- resolve_col(g95_flags, c("direction_matches_expected"), FALSE, "GSE95426 direction_matches_expected")
c_f_pct <- resolve_col(g95_flags, c("empirical_abs_percentile_vs_random"), FALSE, "GSE95426 empirical percentile")
c_f_stronger <- resolve_col(g95_flags, c("stronger_than_random_max"), FALSE, "GSE95426 stronger_than_random_max")
c_f_bias <- resolve_col(g95_flags, c("random_controls_directionally_biased"), FALSE, "GSE95426 random bias")
c_f_bias_dir <- resolve_col(g95_flags, c("random_bias_direction"), FALSE, "GSE95426 random bias direction")

# ---------- Resolve GSE115440 columns ----------
c_g115_module <- resolve_col(g115, c("module"), TRUE, "GSE115440 module")
c_g115_class <- resolve_col(g115, c("module_class", "class"), TRUE, "GSE115440 module_class")
c_g115_delta <- resolve_col(
  g115,
  c("delta_surgery_vs_ctrl", "delta_surgery_vs_control", "delta", "mean_diff"),
  TRUE,
  "GSE115440 delta"
)
c_g115_d <- resolve_col(g115, c("cohens_d", "cohen_d", "d"), TRUE, "GSE115440 Cohen's d")
c_g115_dir <- resolve_col(g115, c("direction"), TRUE, "GSE115440 direction")
c_g115_rescue <- resolve_col(g115, c("rescue_consistency", "rescue_flag"), FALSE, "GSE115440 rescue")
c_g115_cov <- resolve_col(g115, c("coverage_pct", "coverage", "detected_pct", "detected_genes"), FALSE, "GSE115440 coverage")

# ---------- Standardize GSE95426 ----------
g95_std <- tibble(
  dataset = "GSE95426",
  module = g95[[c_g95_module]],
  module_class = g95[[c_g95_class]],
  disease_delta = g95[[c_g95_delta]],
  disease_cohens_d = g95[[c_g95_d]],
  original_direction = g95[[c_g95_dir]],
  disease_direction = harmonize_direction(g95[[c_g95_dir]]),
  caveat = if (!is.na(c_g95_caveat)) g95[[c_g95_caveat]] else NA_character_,
  coverage = if (!is.na(c_g95_cov)) g95[[c_g95_cov]] else NA_real_
)

g95_flags_std <- tibble(
  module = g95_flags[[c_f_module]],
  interpretation_flag = g95_flags[[c_f_interp]],
  expected_direction = if (!is.na(c_f_expected)) g95_flags[[c_f_expected]] else NA_character_,
  direction_matches_expected = if (!is.na(c_f_match)) g95_flags[[c_f_match]] else NA,
  empirical_abs_percentile_vs_random = if (!is.na(c_f_pct)) g95_flags[[c_f_pct]] else NA_real_,
  stronger_than_random_max = if (!is.na(c_f_stronger)) g95_flags[[c_f_stronger]] else NA,
  random_controls_directionally_biased = if (!is.na(c_f_bias)) g95_flags[[c_f_bias]] else NA,
  random_bias_direction = if (!is.na(c_f_bias_dir)) g95_flags[[c_f_bias_dir]] else NA_character_
)

g95_std <- g95_std %>%
  left_join(g95_flags_std, by = "module")

# ---------- Standardize GSE115440 ----------
g115_std <- tibble(
  dataset = "GSE115440",
  module = g115[[c_g115_module]],
  module_class = g115[[c_g115_class]],
  disease_delta = g115[[c_g115_delta]],
  disease_cohens_d = g115[[c_g115_d]],
  original_direction = g115[[c_g115_dir]],
  disease_direction = harmonize_direction(g115[[c_g115_dir]]),
  caveat = NA_character_,
  coverage = if (!is.na(c_g115_cov)) g115[[c_g115_cov]] else NA_real_,
  rescue_consistency = if (!is.na(c_g115_rescue)) g115[[c_g115_rescue]] else NA_character_
)

# ---------- Exclude negative controls ----------
g95_bio <- g95_std %>%
  filter(!is_random_or_negctrl(module, module_class))

g115_bio <- g115_std %>%
  filter(!is_random_or_negctrl(module, module_class))

log_msg("GSE95426 biological modules: ", nrow(g95_bio), " of ", nrow(g95_std))
log_msg("GSE115440 biological modules: ", nrow(g115_bio), " of ", nrow(g115_std))

# ---------- Join ----------
joined <- inner_join(
  g95_bio %>% rename_with(~ paste0(., "_g95"), -c(module, module_class)),
  g115_bio %>% rename_with(~ paste0(., "_g115"), -c(module, module_class)),
  by = c("module", "module_class")
)

if (nrow(joined) == 0) {
  stop("No overlapping biological modules between GSE95426 and GSE115440.")
}

log_msg("Overlapping biological modules: ", nrow(joined))

# ---------- Concordance ----------
concordance <- joined %>%
  rowwise() %>%
  mutate(
    direction_concordance = classify_concordance(disease_direction_g95, disease_direction_g115),
    final_interpretation = final_interp(direction_concordance, interpretation_flag_g95),
    module_family = classify_module_family(module_class, module)
  ) %>%
  ungroup()

concordance_out <- concordance %>%
  transmute(
    module,
    module_class,
    module_family,
    g95_disease_delta = disease_delta_g95,
    g95_disease_cohens_d = disease_cohens_d_g95,
    g95_original_direction = original_direction_g95,
    g95_disease_direction = disease_direction_g95,
    g95_caveat = caveat_g95,
    g95_interpretation_flag = interpretation_flag_g95,
    g95_expected_direction = expected_direction_g95,
    g95_direction_matches_expected = direction_matches_expected_g95,
    g95_empirical_abs_percentile_vs_random = empirical_abs_percentile_vs_random_g95,
    g95_stronger_than_random_max = stronger_than_random_max_g95,
    g95_random_controls_directionally_biased = random_controls_directionally_biased_g95,
    g95_random_bias_direction = random_bias_direction_g95,
    g115_disease_delta = disease_delta_g115,
    g115_disease_cohens_d = disease_cohens_d_g115,
    g115_original_direction = original_direction_g115,
    g115_disease_direction = disease_direction_g115,
    g115_rescue_consistency = rescue_consistency_g115,
    direction_concordance,
    final_interpretation
  ) %>%
  arrange(module_family, module)

readr::write_csv(concordance_out, out_concordance_csv)
log_msg("Wrote concordance table: ", out_concordance_csv)

# ---------- Summaries ----------
dir_summary <- concordance_out %>%
  count(module_family, final_interpretation, name = "n_modules") %>%
  arrange(module_family, desc(n_modules))

readr::write_csv(dir_summary, out_directionality_csv)
log_msg("Wrote directionality summary: ", out_directionality_csv)

flags_out <- concordance_out %>%
  transmute(
    module,
    module_class,
    module_family,
    g95_disease_direction,
    g115_disease_direction,
    g95_interpretation_flag,
    direction_concordance,
    final_interpretation
  )

readr::write_csv(flags_out, out_flags_csv)
log_msg("Wrote interpretation flags: ", out_flags_csv)

# ---------- Conservative summary ----------
infl_rows <- concordance_out %>%
  filter(module_family == "inflammatory")

mito_rows <- concordance_out %>%
  filter(module_family == "mitochondrial")

infl_g115_up <- sum(infl_rows$g115_disease_direction == "up", na.rm = TRUE)
infl_g95_up <- sum(infl_rows$g95_disease_direction == "up", na.rm = TRUE)
mito_g95_down <- sum(mito_rows$g95_disease_direction == "down", na.rm = TRUE)
mito_caveated <- sum(mito_rows$final_interpretation == "caveated_hypothesis_generating", na.rm = TRUE)

summary_lines <- c(
  "Cross-dataset module concordance: conservative summary",
  "=======================================================",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Scope:",
  "- Anchor: GSE95426, POCD vs Control.",
  "- Validation: GSE115440, Surgery vs Control.",
  "- Negative controls and RANDOM_CONTROL_200 excluded from biological concordance.",
  "",
  "Key findings:",
  "1. GSE115440 provides the clearer inflammatory activation evidence.",
  sprintf("   Inflammatory modules up in GSE115440: %d", infl_g115_up),
  sprintf("   Inflammatory modules up in GSE95426: %d", infl_g95_up),
  "   Therefore, GSE95426 does not replicate inflammatory activation.",
  "",
  "2. GSE95426 shows OXPHOS / respiratory electron transport down-shift,",
  "   but Script 11.5 random-control QC indicates broad downward bias.",
  sprintf("   Mitochondrial modules down in GSE95426: %d", mito_g95_down),
  sprintf("   Mitochondrial modules flagged caveated_hypothesis_generating: %d", mito_caveated),
  "   Therefore, mitochondrial down-shift is hypothesis-generating PASS_WITH_CAVEAT.",
  "",
  "3. Overall cross-dataset concordance is partial, with informative discordance",
  "   and dataset-specific caveats.",
  "",
  "4. This supports a conservative technical-report narrative, not a preprint-level",
  "   clean replicated mechanism claim.",
  "",
  "Hard rules respected:",
  "- No raw expression merged across datasets.",
  "- DEG not used as primary evidence.",
  "- Descriptive p-values not used as primary evidence.",
  "- No causal claims.",
  "- No AI/ML claims.",
  "- GSE174412 not touched.",
  "- Surgery+Maresin1 not used as primary validation."
)

writeLines(summary_lines, out_summary_txt)
log_msg("Wrote conservative summary: ", out_summary_txt)

# ---------- Figure 1: effect size heatmap ----------
heat_df <- concordance_out %>%
  select(module, module_family, g95_disease_cohens_d, g115_disease_cohens_d) %>%
  pivot_longer(
    cols = c(g95_disease_cohens_d, g115_disease_cohens_d),
    names_to = "dataset",
    values_to = "cohens_d"
  ) %>%
  mutate(
    dataset = recode(
      dataset,
      g95_disease_cohens_d = "GSE95426 (POCD)",
      g115_disease_cohens_d = "GSE115440 (Surgery)"
    )
  )

p_heat <- ggplot(
  heat_df,
  aes(x = dataset, y = module, fill = cohens_d)
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    na.value = "grey85",
    name = "Cohen's d"
  ) +
  facet_grid(module_family ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Cross-dataset module effect sizes",
    subtitle = "Disease vs control; conservative partial-concordance analysis",
    x = NULL,
    y = "Module"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0)
  )

ggsave(out_heatmap_png, p_heat, width = 7, height = max(4, 0.35 * nrow(concordance_out) + 2), dpi = 200, limitsize = FALSE)
log_msg("Wrote heatmap: ", out_heatmap_png)

# ---------- Figure 2: direction tileplot ----------
tile_df <- concordance_out %>%
  select(
    module,
    module_family,
    `GSE95426 (POCD)` = g95_disease_direction,
    `GSE115440 (Surgery)` = g115_disease_direction,
    final_interpretation
  ) %>%
  pivot_longer(
    cols = c(`GSE95426 (POCD)`, `GSE115440 (Surgery)`),
    names_to = "dataset",
    values_to = "disease_direction"
  ) %>%
  mutate(
    disease_direction = factor(disease_direction, levels = c("up", "near_null", "down"))
  )

p_tile <- ggplot(tile_df, aes(x = dataset, y = module, fill = disease_direction)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(up = "#b2182b", near_null = "grey85", down = "#2166ac"),
    na.value = "grey95",
    name = "Direction"
  ) +
  facet_grid(module_family ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Cross-dataset module direction",
    subtitle = "Harmonized disease direction: GSE95426 POCD and GSE115440 Surgery",
    x = NULL,
    y = "Module"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0)
  )

ggsave(out_tileplot_png, p_tile, width = 7, height = max(4, 0.35 * nrow(concordance_out) + 2), dpi = 200, limitsize = FALSE)
log_msg("Wrote tileplot: ", out_tileplot_png)

# ---------- Sanity assertions ----------
stopifnot(!any(concordance_out$module == "RANDOM_CONTROL_200"))
stopifnot(!any(tolower(as.character(concordance_out$module_class)) == "negative_control"))

log_msg("Sanity checks passed: no RANDOM_CONTROL_200 or negative_control rows in concordance table.")

log_msg("Session info:")
print(sessionInfo())

log_msg("Script 05 complete.")
