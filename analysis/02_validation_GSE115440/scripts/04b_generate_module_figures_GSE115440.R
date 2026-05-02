suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

proj_root <- "."
val_root <- file.path(proj_root, "analysis", "02_validation_GSE115440")

tables_dir <- file.path(val_root, "results", "tables")
figures_dir <- file.path(val_root, "results", "figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

in_meta <- file.path(proj_root, "data", "metadata", "GSE115440_metadata_verified.csv")
in_scores <- file.path(tables_dir, "GSE115440_module_scores_mean_zscore.csv")
in_overlap <- file.path(tables_dir, "GSE115440_module_gene_overlap.csv")
in_primary <- file.path(tables_dir, "GSE115440_module_surgery_vs_control_results.csv")

if (!file.exists(in_meta)) stop("Missing metadata: ", in_meta)
if (!file.exists(in_scores)) stop("Missing module scores: ", in_scores)
if (!file.exists(in_overlap)) stop("Missing overlap table: ", in_overlap)
if (!file.exists(in_primary)) stop("Missing primary result table: ", in_primary)

message(">> Loading module scoring outputs")

meta_raw <- read_csv(in_meta, show_col_types = FALSE)
scores_tbl <- read_csv(in_scores, show_col_types = FALSE)
overlap_tbl <- read_csv(in_overlap, show_col_types = FALSE)
primary_res <- read_csv(in_primary, show_col_types = FALSE)

sid_col <- intersect(
  c("sample_geo_accession", "geo_accession", "sample_id", "sample", "GSM", "Sample"),
  colnames(meta_raw)
)

grp_col <- intersect(
  c("group", "condition", "treatment", "Group", "Condition"),
  colnames(meta_raw)
)

if (length(sid_col) == 0 || length(grp_col) == 0) {
  stop(
    "Cannot identify metadata columns. Found columns: ",
    paste(colnames(meta_raw), collapse = ", ")
  )
}

sid_col <- sid_col[1]
grp_col <- grp_col[1]

meta <- meta_raw %>%
  transmute(
    sample_id = as.character(.data[[sid_col]]),
    group = as.character(.data[[grp_col]])
  )

preferred_order <- c("Control", "Surgery", "Surgery+Maresin1")

meta <- meta %>%
  mutate(
    group = factor(
      group,
      levels = c(
        intersect(preferred_order, unique(group)),
        setdiff(unique(group), preferred_order)
      )
    )
  )

sample_cols <- meta$sample_id

message(">> Generating heatmap")

pheat_mat <- scores_tbl %>%
  column_to_rownames("module") %>%
  as.matrix()

mode(pheat_mat) <- "numeric"

existing_samples <- intersect(sample_cols, colnames(pheat_mat))

if (length(existing_samples) == 0) {
  stop("No metadata samples found in module score table.")
}

pheat_mat <- pheat_mat[, existing_samples, drop = FALSE]

ann_col <- meta %>%
  filter(sample_id %in% existing_samples) %>%
  column_to_rownames("sample_id") %>%
  select(group)

ann_row <- overlap_tbl %>%
  select(module, module_class) %>%
  filter(module %in% rownames(pheat_mat)) %>%
  column_to_rownames("module")

ann_row <- ann_row[rownames(pheat_mat), , drop = FALSE]

png(
  filename = file.path(figures_dir, "GSE115440_module_score_heatmap.png"),
  width = 2400,
  height = max(1400, 220 * nrow(pheat_mat)),
  res = 300
)

pheatmap(
  pheat_mat,
  annotation_col = ann_col,
  annotation_row = ann_row,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 7,
  fontsize_col = 8,
  main = "GSE115440 module scores"
)

dev.off()

message(">> Generating boxplots")

scores_long <- scores_tbl %>%
  pivot_longer(-module, names_to = "sample_id", values_to = "score") %>%
  left_join(meta, by = "sample_id") %>%
  left_join(
    overlap_tbl %>% select(module, module_class),
    by = "module"
  )

plot_long <- scores_long %>%
  mutate(module = factor(module, levels = primary_res$module))

p_box <- ggplot(plot_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.9) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 7)
  ) +
  labs(
    x = NULL,
    y = "Module score",
    title = "GSE115440 module scores by group",
    subtitle = "n = 3 per group; directionality-focused exploratory analysis"
  )

ggsave(
  filename = file.path(figures_dir, "GSE115440_module_score_boxplots.png"),
  plot = p_box,
  width = 10,
  height = max(5, 1.5 * ceiling(length(unique(plot_long$module)) / 3)),
  dpi = 200
)

message(">> Generating delta barplot")

plot_delta <- primary_res %>%
  mutate(
    module = factor(
      module,
      levels = module[order(delta_surgery_vs_ctrl)]
    )
  )

p_bar <- ggplot(
  plot_delta,
  aes(
    x = module,
    y = delta_surgery_vs_ctrl,
    fill = module_class
  )
) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 10) +
  labs(
    x = NULL,
    y = "Delta module score: Surgery - Control",
    fill = "Module class",
    title = "GSE115440 Surgery vs Control module shift",
    subtitle = "Effect direction only; n = 3 vs 3"
  )

ggsave(
  filename = file.path(figures_dir, "GSE115440_module_delta_barplot.png"),
  plot = p_bar,
  width = 9,
  height = max(4, 0.3 * nrow(plot_delta) + 2),
  dpi = 200
)

message("Done. Figures written to: ", figures_dir)
