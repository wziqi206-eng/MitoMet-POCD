#!/usr/bin/env Rscript

# Script 06: Rank-based cross-dataset concordance sensitivity analysis
#
# Purpose:
#   Conservative final analysis for the MitoMet-POCD / NeuroMitoMap packaging phase.
#   Uses only Script 05 output and only the six predefined biological target modules.
#
# Inputs:
#   analysis/03_cross_dataset_concordance/results/tables/cross_dataset_module_concordance_table.csv
#
# Outputs:
#   analysis/03_cross_dataset_concordance/results/tables/cross_dataset_rank_concordance.csv
#   analysis/03_cross_dataset_concordance/results/tables/cross_dataset_rank_concordance.txt
#
# Guardrails:
#   - No raw expression merging.
#   - No new datasets.
#   - No DEG-based claim.
#   - No inflammation-first / mechanism / discovery language.
#   - Concordance is descriptive and conservative.

options(stringsAsFactors = FALSE)

set.seed(20260504)

input_file <- "analysis/03_cross_dataset_concordance/results/tables/cross_dataset_module_concordance_table.csv"
out_csv <- "analysis/03_cross_dataset_concordance/results/tables/cross_dataset_rank_concordance.csv"
out_txt <- "analysis/03_cross_dataset_concordance/results/tables/cross_dataset_rank_concordance.txt"

target_modules <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "MICROGLIA_ACTIVATION_CUSTOM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"
)

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

df <- read.csv(input_file, check.names = FALSE)

required_base <- "module"
if (!required_base %in% names(df)) {
  stop("Input table must contain a 'module' column.")
}

possible_gse95426_d <- c(
  "cohens_d_GSE95426",
  "GSE95426_cohens_d",
  "cohens_d_95426",
  "d_GSE95426",
  "g95_disease_cohens_d"
)

possible_gse115440_d <- c(
  "cohens_d_GSE115440",
  "GSE115440_cohens_d",
  "cohens_d_115440",
  "d_GSE115440",
  "g115_disease_cohens_d"
)

find_col <- function(possible_names, available_names, label) {
  hit <- intersect(possible_names, available_names)
  if (length(hit) == 0) {
    stop(
      "Could not find column for ", label, ". Tried: ",
      paste(possible_names, collapse = ", "),
      "\nAvailable columns: ",
      paste(available_names, collapse = ", ")
    )
  }
  hit[[1]]
}

d1_col <- find_col(possible_gse95426_d, names(df), "GSE95426 Cohen's d")
d2_col <- find_col(possible_gse115440_d, names(df), "GSE115440 Cohen's d")

bio <- df[df$module %in% target_modules, , drop = FALSE]

missing_modules <- setdiff(target_modules, bio$module)
if (length(missing_modules) > 0) {
  stop("Missing target modules from Script 05 table: ", paste(missing_modules, collapse = ", "))
}

if (nrow(bio) != 6) {
  stop("Expected exactly 6 biological target modules, found: ", nrow(bio))
}

bio <- bio[match(target_modules, bio$module), , drop = FALSE]

d_gse95426 <- as.numeric(bio[[d1_col]])
d_gse115440 <- as.numeric(bio[[d2_col]])

if (anyNA(d_gse95426) || anyNA(d_gse115440)) {
  stop("Cohen's d columns contain NA after numeric conversion.")
}

safe_spearman <- function(x, y) {
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, method = "spearman"))
}

spearman_rho <- safe_spearman(d_gse95426, d_gse115440)

boot_n <- 1000
boot_rho <- numeric(boot_n)

for (i in seq_len(boot_n)) {
  idx <- sample(seq_along(d_gse95426), size = length(d_gse95426), replace = TRUE)
  boot_rho[[i]] <- safe_spearman(d_gse95426[idx], d_gse115440[idx])
}

valid_boot <- boot_rho[!is.na(boot_rho)]
if (length(valid_boot) == 0) {
  ci_low <- NA_real_
  ci_high <- NA_real_
} else {
  ci <- quantile(valid_boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  ci_low <- ci[[1]]
  ci_high <- ci[[2]]
}

direction_gse95426 <- sign(d_gse95426)
direction_gse115440 <- sign(d_gse115440)

direction_agree <- direction_gse95426 == direction_gse115440
direction_agree[direction_gse95426 == 0 | direction_gse115440 == 0] <- FALSE

n_agree <- sum(direction_agree)
n_modules <- length(target_modules)

sign_test <- binom.test(n_agree, n_modules, p = 0.5, alternative = "two.sided")

result <- data.frame(
  analysis = "rank_concordance_sign_test",
  input_file = input_file,
  module_n = n_modules,
  spearman_rho = spearman_rho,
  bootstrap_n = boot_n,
  bootstrap_valid_n = length(valid_boot),
  spearman_bootstrap_ci_low = ci_low,
  spearman_bootstrap_ci_high = ci_high,
  direction_agree_n = n_agree,
  direction_agree_fraction = n_agree / n_modules,
  sign_test_p_value = sign_test$p.value,
  interpretation = "descriptive_partial_concordance_only",
  stringsAsFactors = FALSE
)

write.csv(result, out_csv, row.names = FALSE)

txt <- c(
  "Script 06: Rank-based concordance and sign test",
  "================================================",
  "",
  paste0("Input file: ", input_file),
  paste0("GSE95426 effect-size column: ", d1_col),
  paste0("GSE115440 effect-size column: ", d2_col),
  paste0("Biological target modules used: ", n_modules),
  "",
  "Modules:",
  paste0("  - ", bio$module),
  "",
  "Effect sizes:",
  paste0(
    "  - ", bio$module,
    ": GSE95426 d = ", round(d_gse95426, 4),
    "; GSE115440 d = ", round(d_gse115440, 4),
    "; direction_agree = ", direction_agree
  ),
  "",
  paste0("Spearman rho: ", round(spearman_rho, 4)),
  paste0(
    "Bootstrap 95% CI: [",
    round(ci_low, 4), ", ",
    round(ci_high, 4), "]"
  ),
  paste0("Bootstrap valid replicates: ", length(valid_boot), " / ", boot_n),
  "",
  paste0("Direction agreement: ", n_agree, " / ", n_modules),
  paste0("Binomial sign test p-value: ", signif(sign_test$p.value, 4)),
  "",
  "Conservative interpretation:",
  "  This analysis is descriptive. It tests whether the six predefined biological",
  "  modules show rank-level and direction-level concordance across the two",
  "  independently processed datasets. It does not support raw-expression merging,",
  "  mechanism discovery, AI/ML claims, or clean cross-dataset validation language.",
  "",
  "sessionInfo():",
  paste(capture.output(sessionInfo()), collapse = "\n")
)

writeLines(txt, out_txt)

cat("Script 06 complete.\n")
cat("Wrote: ", out_csv, "\n", sep = "")
cat("Wrote: ", out_txt, "\n", sep = "")
