#!/usr/bin/env Rscript

# Stage B.5 target module generator
# Output: data/modules/stage_B5_target_modules.csv

required_pkgs <- c("msigdbr", "dplyr", "readr", "stringr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  cmds <- paste0('install.packages("', missing_pkgs, '")')
  stop(
    paste0(
      "Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\n",
      "Install with:\n", paste(cmds, collapse = "\n")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

find_project_root <- function(start_dir = getwd()) {
  cur <- normalizePath(start_dir)
  repeat {
    if (file.exists(file.path(cur, ".git"))) return(cur)
    parent <- dirname(cur)
    if (identical(cur, parent)) stop("Could not locate project root (.git).", call. = FALSE)
    cur <- parent
  }
}

project_root <- find_project_root(getwd())
out_dir <- file.path(project_root, "data", "modules")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "stage_B5_target_modules.csv")

required_msig_modules <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"
)

module_class_map <- c(
  HALLMARK_INFLAMMATORY_RESPONSE = "inflammatory",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = "inflammatory",
  MICROGLIA_ACTIVATION_CUSTOM = "inflammatory",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = "mitochondrial",
  KEGG_OXIDATIVE_PHOSPHORYLATION = "mitochondrial",
  REACTOME_RESPIRATORY_ELECTRON_TRANSPORT = "mitochondrial"
)

message("Loading MSigDB sets via msigdbr for Mus musculus...")
msig_all <- msigdbr::msigdbr(species = "Mus musculus")

msig_tbl <- msig_all %>%
  filter(.data$gs_name %in% required_msig_modules) %>%
  transmute(
    module = .data$gs_name,
    gene_symbol = str_trim(.data$gene_symbol)
  ) %>%
  filter(!is.na(.data$gene_symbol), .data$gene_symbol != "") %>%
  distinct(.data$module, .data$gene_symbol)

found_modules <- sort(unique(msig_tbl$module))
missing_modules <- setdiff(required_msig_modules, found_modules)
if (length(missing_modules) > 0) {
  stop(
    paste0(
      "msigdbr did not return required module(s): ",
      paste(missing_modules, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

# Source copied exactly from:
# analysis/02_validation_GSE115440/scripts/04_module_scoring_GSE115440.R lines 300-305
microglia_custom <- c(
  "Aif1", "Itgam", "Cd68", "Cd86", "Tlr2", "Tlr4",
  "Tnf", "Il1b", "Il6", "Nos2", "Nfkb1",
  "Ccl2", "Ccl3", "Ccl4", "Cxcl10", "Ptgs2", "Cybb",
  "Trem2", "Tyrobp", "Csf1r", "C1qa", "C1qb", "C1qc"
)

microglia_tbl <- tibble::tibble(
  module = "MICROGLIA_ACTIVATION_CUSTOM",
  gene_symbol = microglia_custom
)

out_tbl <- bind_rows(msig_tbl, microglia_tbl) %>%
  mutate(
    module_class = unname(module_class_map[.data$module]),
    gene_symbol = str_trim(.data$gene_symbol)
  ) %>%
  filter(!is.na(.data$gene_symbol), .data$gene_symbol != "") %>%
  distinct(.data$module, .data$gene_symbol, .keep_all = TRUE) %>%
  arrange(.data$module, .data$gene_symbol) %>%
  select(.data$module, .data$gene_symbol, .data$module_class)

missing_class <- out_tbl %>% filter(is.na(.data$module_class)) %>% pull(.data$module) %>% unique()
if (length(missing_class) > 0) {
  stop(paste0("Missing module_class mapping for: ", paste(missing_class, collapse = ", ")), call. = FALSE)
}

write_csv(out_tbl, out_file)
message("Wrote: ", out_file)
print(out_tbl %>% count(.data$module, .data$module_class, name = "n_genes"))
