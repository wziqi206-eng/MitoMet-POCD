# ============================================================
# Script 04 -- Module-level scoring for GSE115440
# Project : MitoMet-POCD / NeuroMitoMap
# Dataset : GSE115440 (mouse hippocampus)
# Platform : GPL11533 -- Affymetrix Mouse Gene 1.1 ST
# Samples : n = 3 per group (Control / Surgery / Surgery+Maresin1)
# Primary : Surgery vs Control (module-level directionality)
# Auxiliary : Surgery+Maresin1 rescue direction (NOT primary evidence)
#
# Design principles
# - Module-level directionality is the primary readout
# - DEG / single-gene evidence is treated as secondary
# - Honest small-n language; no causal overclaim
# - No raw merging with GSE95426
# - Includes random negative-control modules with fixed seed
# ============================================================

required_cran <- c(
  "dplyr", "tidyr", "readr", "tibble", "stringr",
  "ggplot2", "pheatmap", "RColorBrewer", "matrixStats", "msigdbr"
)

missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

required_bioc <- c("AnnotationDbi", "mogene11sttranscriptcluster.db")
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(20251215)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

proj_root <- "."

val_root <- file.path(proj_root, "analysis", "02_validation_GSE115440")

in_meta <- file.path(
  proj_root,
  "data",
  "metadata",
  "GSE115440_metadata_verified.csv"
)

in_expr <- file.path(
  val_root,
  "results",
  "tables",
  "GSE115440_expression_matrix_rma_probe_level.csv"
)

tables_dir <- file.path(val_root, "results", "tables")
figures_dir <- file.path(val_root, "results", "figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Input checks
# ------------------------------------------------------------

if (!file.exists(in_meta)) {
  stop("Missing metadata file: ", in_meta)
}

if (!file.exists(in_expr)) {
  stop("Missing probe-level expression file: ", in_expr)
}

# ------------------------------------------------------------
# 2. Load inputs
# ------------------------------------------------------------

message(">> [1/9] Loading metadata and probe-level expression matrix")

meta_raw <- read_csv(in_meta, show_col_types = FALSE)
expr_raw <- read_csv(in_expr, show_col_types = FALSE)

sid_col <- intersect(
  c("sample_geo_accession", "geo_accession", "sample_id", "sample", "GSM", "Sample"),
  colnames(meta_raw)
)

sid_col <- if (length(sid_col)) sid_col[1] else NA_character_

grp_col <- intersect(
  c("group", "condition", "treatment", "Group", "Condition"),
  colnames(meta_raw)
)

grp_col <- if (length(grp_col)) grp_col[1] else NA_character_

if (is.na(sid_col) || is.na(grp_col)) {
  stop(
    "Cannot find sample-id / group columns in metadata. ",
    "Looked for sample id in {geo_accession, sample_id, sample, GSM, Sample} ",
    "and group in {group, condition, treatment, Group, Condition}. ",
    "Found columns: ",
    paste(colnames(meta_raw), collapse = ", ")
  )
}

meta <- meta_raw %>%
  transmute(
    sample_id = as.character(.data[[sid_col]]),
    group = as.character(.data[[grp_col]])
  )

preferred_order <- c(
  "Control", "control",
  "Surgery", "surgery",
  "Surgery+Maresin1", "Surgery_Maresin1",
  "surgery+maresin1", "Maresin1", "MaR1"
)

ordered_levels <- c(
  intersect(preferred_order, unique(meta$group)),
  setdiff(unique(meta$group), preferred_order)
)

meta <- meta %>%
  mutate(group = factor(group, levels = ordered_levels))

message(" samples in metadata: ", nrow(meta))
message(" group levels: ", paste(levels(meta$group), collapse = " | "))

probe_col <- colnames(expr_raw)[1]
message(" probe id column: ", probe_col)

sample_cols <- as.character(meta$sample_id)

missing_samp <- setdiff(sample_cols, colnames(expr_raw))

if (length(missing_samp) > 0) {
  stop(
    "These samples are listed in metadata but missing from expression matrix: ",
    paste(missing_samp, collapse = ", ")
  )
}

expr_probe <- expr_raw %>%
  dplyr::select(all_of(probe_col), all_of(sample_cols))

colnames(expr_probe)[1] <- "probe_id"

expr_probe <- expr_probe %>%
  mutate(probe_id = as.character(probe_id))

# ------------------------------------------------------------
# 3. Probe -> gene symbol annotation
# ------------------------------------------------------------

message(">> [2/9] Mapping GPL11533 probes to mouse gene symbols")

annot_pkg <- "mogene11sttranscriptcluster.db"

have_annot <- requireNamespace(annot_pkg, quietly = TRUE) &&
  requireNamespace("AnnotationDbi", quietly = TRUE)

if (have_annot) {
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(mogene11sttranscriptcluster.db)
  })

  probe_ids <- expr_probe$probe_id

  map_df <- suppressMessages(
    AnnotationDbi::select(
      mogene11sttranscriptcluster.db,
      keys = probe_ids,
      columns = c("SYMBOL"),
      keytype = "PROBEID"
    )
  ) %>%
    as_tibble() %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(PROBEID, SYMBOL)

  message(
    " annotated probes: ",
    n_distinct(map_df$PROBEID),
    " / ",
    length(probe_ids)
  )
} else {
  warning(
    "Annotation package '",
    annot_pkg,
    "' not installed. Falling back to identity mapping."
  )

  map_df <- tibble(
    PROBEID = expr_probe$probe_id,
    SYMBOL = expr_probe$probe_id
  )
}

# ------------------------------------------------------------
# 4. Collapse to gene-level
# ------------------------------------------------------------

message(">> [3/9] Collapsing to gene-level by highest-variance probe per gene")

expr_long <- expr_probe %>%
  inner_join(map_df, by = c("probe_id" = "PROBEID"))

mat_for_var <- as.matrix(expr_long[, sample_cols, drop = FALSE])
mode(mat_for_var) <- "numeric"

expr_long$.var <- matrixStats::rowVars(mat_for_var, na.rm = TRUE)

best_probe <- expr_long %>%
  group_by(SYMBOL) %>%
  slice_max(.var, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_gene <- best_probe %>%
  dplyr::select(gene_symbol = SYMBOL, all_of(sample_cols))

write_csv(
  expr_gene,
  file.path(tables_dir, "GSE115440_expression_matrix_rma_gene_level.csv")
)

message(
  " gene-level matrix: ",
  nrow(expr_gene),
  " genes x ",
  length(sample_cols),
  " samples"
)

# ------------------------------------------------------------
# 5. Define modules
# ------------------------------------------------------------

message(">> [4/9] Defining mitochondrial / inflammatory / negative-control modules")

msig_all <- if (requireNamespace("msigdbr", quietly = TRUE)) {
  tryCatch(
    msigdbr::msigdbr(species = "Mus musculus"),
    error = function(e) {
      message(" msigdbr load failed: ", e$message)
      NULL
    }
  )
} else {
  NULL
}

get_msig <- function(name) {
  if (is.null(msig_all)) {
    return(NULL)
  }

  syms <- msig_all %>%
    filter(.data$gs_name == name) %>%
    pull(.data$gene_symbol) %>%
    unique()

  if (length(syms) == 0) {
    NULL
  } else {
    syms
  }
}

mod_oxphos_h <- get_msig("HALLMARK_OXIDATIVE_PHOSPHORYLATION")
mod_oxphos_kegg <- get_msig("KEGG_OXIDATIVE_PHOSPHORYLATION")
mod_etc_react <- get_msig("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")
mod_inflam_h <- get_msig("HALLMARK_INFLAMMATORY_RESPONSE")
mod_tnfa_h <- get_msig("HALLMARK_TNFA_SIGNALING_VIA_NFKB")

mod_microglia_custom <- c(
  "Aif1", "Itgam", "Cd68", "Cd86", "Tlr2", "Tlr4",
  "Tnf", "Il1b", "Il6", "Nos2", "Nfkb1",
  "Ccl2", "Ccl3", "Ccl4", "Cxcl10", "Ptgs2", "Cybb",
  "Trem2", "Tyrobp", "Csf1r", "C1qa", "C1qb", "C1qc"
)

modules <- list()

add_mod <- function(lst, key, genes) {
  if (!is.null(genes) && length(genes) >= 5) {
    lst[[key]] <- genes
  }
  lst
}

modules <- add_mod(modules, "HALLMARK_OXIDATIVE_PHOSPHORYLATION", mod_oxphos_h)
modules <- add_mod(modules, "KEGG_OXIDATIVE_PHOSPHORYLATION", mod_oxphos_kegg)
modules <- add_mod(modules, "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT", mod_etc_react)
modules <- add_mod(modules, "HALLMARK_INFLAMMATORY_RESPONSE", mod_inflam_h)
modules <- add_mod(modules, "HALLMARK_TNFA_SIGNALING_VIA_NFKB", mod_tnfa_h)
modules <- add_mod(modules, "MICROGLIA_ACTIVATION_CUSTOM", mod_microglia_custom)

if (length(modules) <= 1) {
  warning(
    "msigdbr unavailable or empty; using built-in fallback modules. ",
    "Install msigdbr for full module coverage."
  )

  modules$HALLMARK_OXIDATIVE_PHOSPHORYLATION_fallback <- c(
    "Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5f1", "Atp5g1", "Atp5g2",
    "Atp5g3", "Atp5h", "Atp5j", "Atp5l", "Atp5o", "Cox4i1", "Cox5a", "Cox5b",
    "Cox6a1", "Cox6b1", "Cox6c", "Cox7a2", "Cox7b", "Cox7c", "Cox8a",
    "Ndufa1", "Ndufa2", "Ndufa3", "Ndufa4", "Ndufa5", "Ndufa6", "Ndufa7",
    "Ndufa8", "Ndufb1", "Ndufb2", "Ndufb3", "Ndufb4", "Ndufb5", "Ndufb6",
    "Ndufb7", "Ndufb8", "Ndufb9", "Ndufb10", "Ndufs1", "Ndufs2", "Ndufs3",
    "Ndufs4", "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8", "Ndufv1", "Ndufv2",
    "Ndufv3", "Sdha", "Sdhb", "Sdhc", "Sdhd", "Uqcrc1", "Uqcrc2", "Uqcrfs1",
    "Uqcrh", "Uqcrq", "Cycs"
  )

  modules$HALLMARK_INFLAMMATORY_RESPONSE_fallback <- c(
    "Il1a", "Il1b", "Il6", "Tnf", "Tlr1", "Tlr2", "Tlr4", "Tlr6",
    "Nfkb1", "Nfkb2", "Rela", "Relb", "Ifit1", "Ifit2", "Ifit3",
    "Cxcl1", "Cxcl2", "Cxcl10", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
    "Ptgs2", "Nos2", "Nlrp3", "Casp1", "Stat1", "Stat3", "Mmp9",
    "Cd14", "Cd40", "Cd86", "Icam1", "Vcam1"
  )
}

expressed_genes <- expr_gene$gene_symbol

sizes_pool <- vapply(modules, length, integer(1))
target_size <- as.integer(round(median(sizes_pool)))

n_neg <- 3

for (i in seq_len(n_neg)) {
  modules[[paste0("RANDOM_CONTROL_", i)]] <- sample(
    expressed_genes,
    size = min(target_size, length(expressed_genes))
  )
}

message(" total modules defined: ", length(modules))

# ------------------------------------------------------------
# 6. Module / dataset gene overlap
# ------------------------------------------------------------

message(">> [5/9] Computing module x dataset gene overlap")

overlap_tbl <- tibble(
  module = names(modules),
  module_size = vapply(modules, length, integer(1)),
  detected = vapply(
    modules,
    function(g) length(intersect(g, expressed_genes)),
    integer(1)
  )
) %>%
  mutate(
    detected_frac = round(detected / module_size, 3),
    module_class = case_when(
      grepl("OXIDATIVE|RESPIRATORY|OXPHOS", module) ~ "mitochondrial",
      grepl("INFLAMMATORY|TNFA|MICROGLIA", module) ~ "inflammatory",
      grepl("RANDOM_CONTROL", module) ~ "negative_control",
      TRUE ~ "other"
    )
  )

write_csv(
  overlap_tbl,
  file.path(tables_dir, "GSE115440_module_gene_overlap.csv")
)

keep_mod <- overlap_tbl %>%
  filter(detected >= 5) %>%
  pull(module)

if (length(keep_mod) == 0) {
  stop("No module has >=5 detected genes; check annotation / mapping.")
}

modules <- modules[keep_mod]

overlap_tbl <- overlap_tbl %>%
  filter(module %in% keep_mod)

message(" modules retained with >=5 detected genes: ", length(modules))

# ------------------------------------------------------------
# 7. Per-sample module scores
# ------------------------------------------------------------

message(">> [6/9] Computing mean-z-score module scores")

expr_mat <- expr_gene %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

mode(expr_mat) <- "numeric"

gene_sd <- matrixStats::rowSds(expr_mat, na.rm = TRUE)

expr_mat <- expr_mat[gene_sd > 0, , drop = FALSE]

expr_z <- t(scale(t(expr_mat)))

score_one <- function(genes) {
  g <- intersect(genes, rownames(expr_z))

  if (length(g) < 5) {
    return(rep(NA_real_, ncol(expr_z)))
  }

  colMeans(expr_z[g, , drop = FALSE], na.rm = TRUE)
}

score_mat <- t(vapply(modules, score_one, numeric(ncol(expr_z))))

colnames(score_mat) <- colnames(expr_z)

scores_tbl <- as.data.frame(score_mat) %>%
  rownames_to_column("module") %>%
  as_tibble()

write_csv(
  scores_tbl,
  file.path(tables_dir, "GSE115440_module_scores_mean_zscore.csv")
)

scores_long <- scores_tbl %>%
  pivot_longer(-module, names_to = "sample_id", values_to = "score") %>%
  left_join(meta, by = "sample_id") %>%
  left_join(
    overlap_tbl %>% dplyr::select(module, module_class),
    by = "module"
  )

# ------------------------------------------------------------
# 8. Group summary and Surgery vs Control
# ------------------------------------------------------------

message(">> [7/9] Group summary and Surgery-vs-Control comparison")

group_summary <- scores_long %>%
  group_by(module, module_class, group) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = sum(!is.na(score)),
    .groups = "drop"
  )

write_csv(
  group_summary,
  file.path(tables_dir, "GSE115440_module_group_summary.csv")
)

glabels <- as.character(unique(meta$group))

ctrl_lab <- glabels[grepl("^(Control|control|CTRL|ctrl|Sham|sham)", glabels)]
ctrl_lab <- if (length(ctrl_lab)) ctrl_lab[1] else NA_character_

res_lab <- glabels[grepl("Maresin", glabels, ignore.case = TRUE)]
res_lab <- if (length(res_lab)) res_lab[1] else NA_character_

surg_lab <- glabels[
  grepl("Surgery|surgery", glabels) &
    !grepl("Maresin", glabels, ignore.case = TRUE)
]
surg_lab <- if (length(surg_lab)) surg_lab[1] else NA_character_

if (is.na(ctrl_lab) || is.na(surg_lab)) {
  stop(
    "Cannot identify Control / Surgery groups from metadata. Found: ",
    paste(glabels, collapse = " | ")
  )
}

message(" Control = ", ctrl_lab)
message(" Surgery = ", surg_lab)
message(" Rescue = ", ifelse(is.na(res_lab), "(none)", res_lab))

primary_test <- function(mod_name) {
  s_ctrl <- scores_long %>%
    filter(module == mod_name, group == ctrl_lab) %>%
    pull(score)

  s_surg <- scores_long %>%
    filter(module == mod_name, group == surg_lab) %>%
    pull(score)

  s_res <- if (!is.na(res_lab)) {
    scores_long %>%
      filter(module == mod_name, group == res_lab) %>%
      pull(score)
  } else {
    numeric(0)
  }

  delta_sc <- mean(s_surg, na.rm = TRUE) - mean(s_ctrl, na.rm = TRUE)

  tt <- tryCatch(
    t.test(s_surg, s_ctrl, var.equal = FALSE),
    error = function(e) NULL
  )

  pval <- if (!is.null(tt)) tt$p.value else NA_real_

  pooled_sd <- sqrt(
    (
      stats::var(s_surg, na.rm = TRUE) +
        stats::var(s_ctrl, na.rm = TRUE)
    ) / 2
  )

  cohens_d <- if (is.finite(pooled_sd) && pooled_sd > 0) {
    delta_sc / pooled_sd
  } else {
    NA_real_
  }

  delta_rs <- if (length(s_res) > 0) {
    mean(s_res, na.rm = TRUE) - mean(s_surg, na.rm = TRUE)
  } else {
    NA_real_
  }

  rescue_dir <- if (is.na(delta_rs)) {
    NA_character_
  } else if (delta_rs * delta_sc < 0) {
    "consistent_rescue"
  } else {
    "not_consistent"
  }

  tibble(
    module = mod_name,
    mean_control = mean(s_ctrl, na.rm = TRUE),
    mean_surgery = mean(s_surg, na.rm = TRUE),
    delta_surgery_vs_ctrl = delta_sc,
    cohens_d = cohens_d,
    welch_t_pvalue = pval,
    direction = ifelse(delta_sc > 0, "up_in_surgery", "down_in_surgery"),
    mean_rescue = if (length(s_res) > 0) mean(s_res, na.rm = TRUE) else NA_real_,
    delta_rescue_vs_surg = delta_rs,
    rescue_consistency = rescue_dir
  )
}

primary_res <- bind_rows(lapply(names(modules), primary_test)) %>%
  left_join(
    overlap_tbl %>% dplyr::select(module, module_class, detected, module_size),
    by = "module"
  ) %>%
  mutate(bh_qvalue = p.adjust(welch_t_pvalue, method = "BH")) %>%
  dplyr::select(
    module,
    module_class,
    module_size,
    detected,
    mean_control,
    mean_surgery,
    delta_surgery_vs_ctrl,
    cohens_d,
    welch_t_pvalue,
    bh_qvalue,
    direction,
    mean_rescue,
    delta_rescue_vs_surg,
    rescue_consistency
  ) %>%
  arrange(module_class, desc(abs(cohens_d)))

write_csv(
  primary_res,
  file.path(tables_dir, "GSE115440_module_surgery_vs_control_results.csv")
)

direction_summary <- primary_res %>%
  group_by(module_class) %>%
  summarise(
    n_modules = n(),
    n_up_in_surgery = sum(direction == "up_in_surgery", na.rm = TRUE),
    n_down_in_surgery = sum(direction == "down_in_surgery", na.rm = TRUE),
    mean_delta = mean(delta_surgery_vs_ctrl, na.rm = TRUE),
    median_delta = median(delta_surgery_vs_ctrl, na.rm = TRUE),
    mean_cohens_d = mean(cohens_d, na.rm = TRUE),
    n_rescue_consistent = sum(
      rescue_consistency == "consistent_rescue",
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  direction_summary,
  file.path(tables_dir, "GSE115440_module_directionality_summary.csv")
)

# ------------------------------------------------------------
# 9. Figures
# ------------------------------------------------------------

message(">> [8/9] Generating figures")

pheat_mat <- as.matrix(scores_tbl[, -1, drop = FALSE])
rownames(pheat_mat) <- scores_tbl$module
pheat_mat <- pheat_mat[, sample_cols, drop = FALSE]

ann_col <- data.frame(
  group = meta$group,
  row.names = meta$sample_id
)

ann_row <- data.frame(
  class = overlap_tbl$module_class[
    match(rownames(pheat_mat), overlap_tbl$module)
  ],
  row.names = rownames(pheat_mat)
)

row_can_cluster <- nrow(pheat_mat) >= 2 &&
  all(matrixStats::rowSds(pheat_mat, na.rm = TRUE) > 0, na.rm = TRUE) &&
  all(is.finite(pheat_mat))

pheatmap(
  pheat_mat,
  annotation_col = ann_col,
  annotation_row = ann_row,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
  cluster_rows = row_can_cluster,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "GSE115440 module scores (mean z-score)",
  filename = file.path(figures_dir, "GSE115440_module_score_heatmap.png"),
  width = 8,
  height = max(4, 0.32 * nrow(pheat_mat) + 2)
)

plot_long <- scores_long %>%
  mutate(module = factor(module, levels = primary_res$module))

p_box <- ggplot(plot_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.15, size = 1.3, alpha = 0.9) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 7)
  ) +
  labs(
    x = NULL,
    y = "Module score (mean z-score)",
    title = "GSE115440 -- module scores by group (n = 3 per group)",
    subtitle = "Exploratory; small-n. Treat directionality as primary readout."
  )

ggsave(
  file.path(figures_dir, "GSE115440_module_score_boxplots.png"),
  p_box,
  width = 10,
  height = max(5, 1.6 * ceiling(length(modules) / 3)),
  dpi = 200
)

plot_delta <- primary_res %>%
  mutate(module = factor(module, levels = module[order(delta_surgery_vs_ctrl)]))

p_bar <- ggplot(
  plot_delta,
  aes(x = module, y = delta_surgery_vs_ctrl, fill = module_class)
) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 10) +
  labs(
    x = NULL,
    y = "Delta module score (Surgery - Control)",
    fill = "module class",
    title = "GSE115440 -- Surgery vs Control module shift",
    subtitle = "n = 3 vs 3; bars show effect direction, not significance"
  )

ggsave(
  file.path(figures_dir, "GSE115440_module_delta_barplot.png"),
  p_bar,
  width = 9,
  height = max(4, 0.28 * nrow(plot_delta) + 2),
  dpi = 200
)

# ------------------------------------------------------------
# 10. Console summary
# ------------------------------------------------------------

message(">> [9/9] Module directionality summary")
message("---------------------------------------------------")
message("GSE115440 -- Surgery vs Control module-level summary")
message("Sample size: n = 3 per group; exploratory only")
message("---------------------------------------------------")

print(direction_summary, n = Inf)

message("")
message("Per-module results sorted within class by |Cohen's d|:")

print(
  primary_res %>%
    dplyr::select(
      module,
      module_class,
      delta_surgery_vs_ctrl,
      cohens_d,
      welch_t_pvalue,
      bh_qvalue,
      direction,
      rescue_consistency
    ),
  n = Inf
)

message("")
message("INTERPRETATION CAVEATS")
message(" - n = 3 vs 3: p-values are unstable; use module-level direction as primary.")
message(" - Random-control modules are reference baselines, not formal tests.")
message(" - Surgery+Maresin1 columns are auxiliary rescue evidence, not validation.")
message(" - Cross-dataset confirmation with GSE95426 happens in a downstream step.")
message("Done. Outputs written under: ", val_root)
