#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- normalizePath(sub(file_arg, "", args[grep(file_arg, args)]))
script_dir <- dirname(script_path)

find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir)
  repeat {
    if (file.exists(file.path(cur, ".git"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Could not locate project root (.git).", call. = FALSE)
    cur <- parent
  }
}

project_root <- find_project_root(script_dir)
analysis_root <- file.path(project_root, "analysis", "01_anchor_GSE95426")
results_tables <- file.path(analysis_root, "results", "tables")
results_figures <- file.path(analysis_root, "results", "figures")
logs_dir <- file.path(analysis_root, "logs")
dir.create(results_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(results_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(logs_dir, "GSE95426_stage_B5_module_coverage_QC.log")
log_con <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  writeLines(msg, con = log_con)
}

required_pkgs <- c("ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install_cmd <- sprintf('install.packages(c(%s))', paste(sprintf('"%s"', missing_pkgs), collapse = ", "))
  stop(
    paste0(
      "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
      ". Please install manually in R with: ", install_cmd
    ),
    call. = FALSE
  )
}

required_modules <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "MICROGLIA_ACTIVATION_CUSTOM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"
)
optional_modules <- c("MITOCARTA_3.0_MOUSE", "MITOCARTA_3_MOUSE", "MITOCARTA_MOUSE")

clean_symbols <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x) | x == ""] <- NA_character_
  unique(x[!is.na(x)])
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  modules <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) next
    modules[[parts[1]]] <- clean_symbols(parts[-c(1, 2)])
  }
  modules
}

extract_modules_from_table <- function(df) {
  nms <- names(df)
  nms_upper <- toupper(nms)
  module_col <- which(nms_upper %in% c("MODULE", "GENESET", "GENE_SET", "SET", "PATHWAY"))[1]
  gene_col <- which(nms_upper %in% c("GENE", "GENE_SYMBOL", "SYMBOL", "GENES"))[1]
  out <- list()

  if (!is.na(module_col) && !is.na(gene_col)) {
    mod_vals <- trimws(as.character(df[[module_col]]))
    gene_vals <- trimws(as.character(df[[gene_col]]))
    keep <- !(is.na(mod_vals) | mod_vals == "" | is.na(gene_vals) | gene_vals == "")
    if (any(keep)) {
      split_genes <- split(gene_vals[keep], mod_vals[keep])
      out <- lapply(split_genes, clean_symbols)
    }
    return(out)
  }

  for (col in names(df)) {
    vals <- clean_symbols(df[[col]])
    if (length(vals) > 0) out[[col]] <- vals
  }
  out
}

read_module_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  log_msg("Reading module file: ", path)
  if (ext == "gmt") return(read_gmt(path))

  sep <- if (ext == "csv") "," else "\t"
  df <- tryCatch(
    read.table(path, header = TRUE, sep = sep, quote = '"', comment.char = "", stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) {
    df <- read.table(path, header = TRUE, sep = "", quote = '"', comment.char = "", stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE)
  }
  extract_modules_from_table(df)
}

module_search_dirs <- file.path(project_root, c(
  "data/modules",
  "data/gene_sets",
  "analysis/01_anchor_GSE95426/modules",
  "analysis/01_anchor_GSE95426/data",
  "resources",
  "config"
))

candidate_files <- c()
for (d in module_search_dirs) {
  if (dir.exists(d)) {
    candidate_files <- c(candidate_files, list.files(d, pattern = "\\.(csv|tsv|txt|gmt)$", full.names = TRUE, ignore.case = TRUE))
  }
}
candidate_files <- unique(candidate_files)

preferred_module_file <- file.path(project_root, "data", "modules", "stage_B5_target_modules.csv")
if (file.exists(preferred_module_file)) {
  candidate_files <- c(preferred_module_file, setdiff(candidate_files, preferred_module_file))
  log_msg("Preferred module file found and prioritized: ", preferred_module_file)
}

if (length(candidate_files) == 0) {
  stop(
    paste0(
      "No module files found. Place a module file (.csv/.tsv/.txt/.gmt) under one of: ",
      paste(module_search_dirs, collapse = ", ")
    ),
    call. = FALSE
  )
}

best_modules <- NULL
best_file <- NULL
best_match <- -1
for (f in candidate_files) {
  mods <- tryCatch(read_module_file(f), error = function(e) {
    log_msg("Skipping unreadable module file: ", f, " | ", conditionMessage(e))
    NULL
  })
  if (is.null(mods) || length(mods) == 0) next
  mod_names_upper <- toupper(names(mods))
  req_match <- sum(required_modules %in% mod_names_upper)
  if (req_match > best_match) {
    best_match <- req_match
    best_modules <- mods
    best_file <- f
  }
}

if (is.null(best_modules)) {
  stop("No readable module files found in search locations.", call. = FALSE)
}

names(best_modules) <- toupper(names(best_modules))
log_msg("Selected module file: ", best_file, " (required modules matched: ", best_match, "/", length(required_modules), ")")

if (best_match == 0) {
  stop(
    paste0(
      "Stage B.5 cannot run because standard target module gene lists are not present. ",
      "Provide an authoritative module file at ", preferred_module_file,
      " with columns module,gene_symbol,module_class for required Hallmark/KEGG/Reactome/Microglia modules."
    ),
    call. = FALSE
  )
}

missing_required <- setdiff(required_modules, names(best_modules))
if (length(missing_required) > 0) {
  stop(
    paste0(
      "Required module(s) missing from selected module file: ", paste(missing_required, collapse = ", "),
      ". Add these modules to a file under: ", paste(module_search_dirs, collapse = ", "),
      ". Selected file was: ", best_file
    ),
    call. = FALSE
  )
}

expr_candidates <- c(
  file.path(analysis_root, "results", "tables", "GSE95426_expression_matrix_rescued_protein_coding_gene_level.csv"),
  file.path(analysis_root, "results", "tables", "GSE95426_rescued_protein_coding_gene_level_matrix.csv"),
  file.path(project_root, "data", "processed", "GSE95426_expression_matrix_rescued_protein_coding_gene_level.csv"),
  file.path(project_root, "data", "processed", "GSE95426_expression_matrix_gene_level_protein_coding.csv"),
  file.path(project_root, "data", "processed", "GSE95426_expression_matrix_gene_level.csv")
)
expr_path <- expr_candidates[file.exists(expr_candidates)][1]
if (is.na(expr_path)) {
  stop(
    paste0(
      "Rescued protein-coding gene-level matrix not found. Expected one of: ",
      paste(expr_candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

log_msg("Reading expression matrix: ", expr_path)
expr_df <- read.csv(expr_path, check.names = FALSE, stringsAsFactors = FALSE)
col_upper <- toupper(names(expr_df))
gene_col_idx <- which(col_upper %in% c("GENE_SYMBOL", "GENE", "SYMBOL", "GENESYMBOL", "GENEID"))[1]
if (is.na(gene_col_idx)) gene_col_idx <- 1
universe <- clean_symbols(expr_df[[gene_col_idx]])
if (length(universe) == 0) stop("Gene universe extraction failed (no valid symbols).", call. = FALSE)
log_msg("Gene universe size: ", length(universe))

modules_eval <- best_modules[required_modules]
mito_opt_name <- optional_modules[optional_modules %in% names(best_modules)][1]
if (!is.na(mito_opt_name)) {
  modules_eval[[mito_opt_name]] <- best_modules[[mito_opt_name]]
  log_msg("Optional module detected and included: ", mito_opt_name)
}
set.seed(20260504)
modules_eval[["RANDOM_CONTROL_200"]] <- sample(universe, size = min(200, length(universe)), replace = FALSE)

module_class <- function(m) {
  if (grepl("RANDOM_CONTROL", m)) return("negative_control")
  if (grepl("INFLAMMATORY|TNFA|MICROGLIA", m)) return("inflammation")
  if (grepl("OXIDATIVE|RESPIRATORY|MITOCARTA|ETC", m)) return("mitochondrial")
  "other"
}

status_from_percent <- function(p) {
  if (p >= 85) return("PASS_STRONG")
  if (p >= 70) return("PASS_MIN")
  if (p >= 60) return("EXCLUDED")
  "FORBIDDEN"
}

rows <- list(); miss_rows <- list()
for (m in names(modules_eval)) {
  genes <- clean_symbols(modules_eval[[m]])
  recovered <- intersect(genes, universe)
  missing <- setdiff(genes, universe)
  pct <- if (length(genes) == 0) 0 else round(100 * length(recovered) / length(genes), 2)
  st <- status_from_percent(pct)
  interp <- if (st %in% c("PASS_STRONG", "PASS_MIN")) "YES" else "NO"
  note <- if (m == "RANDOM_CONTROL_200") "Coverage reference only; no group separation used." else ""
  rows[[length(rows)+1]] <- data.frame(
    module = m,
    module_class = module_class(m),
    source_gene_count = length(genes),
    recovered_gene_count = length(recovered),
    missing_gene_count = length(missing),
    recovery_percent = pct,
    status = st,
    interpretation_allowed = interp,
    notes = note,
    stringsAsFactors = FALSE
  )
  if (length(missing) > 0) {
    miss_rows[[length(miss_rows)+1]] <- data.frame(module = m, missing_gene_symbol = missing, stringsAsFactors = FALSE)
  }
}
coverage_tbl <- do.call(rbind, rows)
missing_tbl <- if (length(miss_rows) > 0) do.call(rbind, miss_rows) else data.frame(module=character(), missing_gene_symbol=character())

key_inflam <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
key_strong <- all(coverage_tbl$recovery_percent[match(key_inflam, coverage_tbl$module)] >= 85)
key_lt70 <- any(coverage_tbl$recovery_percent[match(key_inflam, coverage_tbl$module)] < 70)

target_tbl <- coverage_tbl[coverage_tbl$module %in% c(required_modules, optional_modules), , drop = FALSE]
n_target <- nrow(target_tbl)
n_ge70 <- sum(target_tbl$recovery_percent >= 70)
n_lt60 <- sum(target_tbl$recovery_percent < 60)
mito_mask <- grepl("OXIDATIVE|RESPIRATORY|MITOCARTA|ETC", target_tbl$module)
n_mito_lt70 <- sum(target_tbl$recovery_percent[mito_mask] < 70)
n_target_lt70 <- sum(target_tbl$recovery_percent < 70)

verdict <- "PASS_WITH_CAVEAT"
if (key_lt70 || n_lt60 >= 3) {
  verdict <- "FAIL_DO_NOT_SCORE"
} else if (key_strong && n_mito_lt70 >= 1) {
  verdict <- "INFLAMMATION_ONLY_PASS"
} else if ((!key_strong) && n_target_lt70 >= 2) {
  verdict <- "REBUILD_ANNOTATION"
} else if (!(n_target >= 6 && n_ge70 >= 6 && key_strong)) {
  verdict <- "PASS_WITH_CAVEAT"
}

neg_status <- coverage_tbl$status[coverage_tbl$module == "RANDOM_CONTROL_200"]
if (length(neg_status) == 0) neg_status <- "NOT_GENERATED"
next_action <- switch(verdict,
  "FAIL_DO_NOT_SCORE" = "Do not score modules; fix mapping/module definitions first.",
  "REBUILD_ANNOTATION" = "Rebuild annotation rescue and remap to protein-coding genes.",
  "INFLAMMATION_ONLY_PASS" = "Proceed only with inflammation modules; exclude low-coverage mito modules.",
  "PASS_WITH_CAVEAT" = "Proceed with modules >=70% only, document exclusions.",
  "Manual review required."
)

summary_tbl <- data.frame(
  final_verdict = verdict,
  key_inflammation_pass_strong = key_strong,
  n_target_modules_total = n_target,
  n_modules_ge_70 = n_ge70,
  n_modules_lt_60 = n_lt60,
  n_mito_modules_lt_70 = n_mito_lt70,
  negative_control_status = neg_status,
  next_action = next_action,
  stringsAsFactors = FALSE
)

write.csv(coverage_tbl, file.path(results_tables, "GSE95426_stage_B5_module_coverage_QC.csv"), row.names = FALSE)
write.csv(missing_tbl, file.path(results_tables, "GSE95426_stage_B5_missing_genes_by_module.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(results_tables, "GSE95426_stage_B5_decision_summary.csv"), row.names = FALSE)

plot_tbl <- coverage_tbl
plot_tbl$module <- factor(plot_tbl$module, levels = plot_tbl$module[order(plot_tbl$recovery_percent)])

p <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = recovery_percent, y = module, fill = status)) +
  ggplot2::geom_col() +
  ggplot2::geom_vline(xintercept = c(60, 70, 85), linetype = "dashed", color = c("#b22222", "#e67e22", "#2e8b57")) +
  ggplot2::scale_fill_manual(values = c(PASS_STRONG="#2e8b57", PASS_MIN="#4c9f70", EXCLUDED="#e67e22", FORBIDDEN="#b22222")) +
  ggplot2::labs(
    title = "GSE95426 Stage B.5 Module Coverage QC",
    x = "Recovery percent (%)",
    y = "Module",
    fill = "Status"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  filename = file.path(results_figures, "GSE95426_stage_B5_module_coverage_barplot.png"),
  plot = p, width = 10, height = 6, dpi = 300
)

log_msg("Coverage QC completed.")
log_msg("Final verdict: ", verdict)
print(summary_tbl)
