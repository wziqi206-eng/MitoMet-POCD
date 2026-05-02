# ============================================================
# Script: 01_download_GSE115440_raw.R
# Project: MitoMet-POCD / NeuroMitoMap
# Stage  : Validation cohort (Tier 2) - GSE115440
# Purpose: Download GSE115440 supplementary RAW tar from GEO
#          and extract CEL files. NO normalization here.
#
# IMPORTANT: This script must NOT touch GSE95426 (Tier 1, frozen).
# ============================================================

cat("============================================================\n")
cat(" GSE115440 :: Step 01 - Download raw CEL files\n")
cat("============================================================\n\n")

# ---- 0. Sanity checks -------------------------------------------------
if (!dir.exists("analysis")) {
  warning(
    "Working directory does not contain 'analysis/'. ",
    "Please run this script from the PROJECT ROOT.\n",
    "  Current wd: ", getwd()
  )
}

# ---- 1. Configuration -------------------------------------------------
gse_id <- "GSE115440"
gse_stub <- "GSE115nnn"
expected_n_samples <- 9L

raw_dir <- file.path("data", "raw", gse_id)
cel_dir <- file.path(raw_dir, "CEL")
tar_filename <- paste0(gse_id, "_RAW.tar")
tar_path <- file.path(raw_dir, tar_filename)

geo_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  gse_stub, "/", gse_id, "/suppl/", tar_filename
)

# Hard guard: refuse to operate on GSE95426 paths
if (any(grepl("GSE95426", c(raw_dir, cel_dir, tar_path)))) {
  stop("Refusing to run: this script must not modify GSE95426.")
}

# ---- 2. Create directories -------------------------------------------
cat("[1/5] Preparing directories ...\n")
for (d in c(raw_dir, cel_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    cat("    created : ", d, "\n", sep = "")
  } else {
    cat("    exists  : ", d, "\n", sep = "")
  }
}
cat("\n")

# ---- 3. Download tar archive -----------------------------------------
cat("[2/5] Downloading raw archive from GEO ...\n")
cat("    URL  : ", geo_url, "\n", sep = "")
cat("    Save : ", tar_path, "\n", sep = "")

options(timeout = max(3600, getOption("timeout")))

download_ok <- FALSE

if (file.exists(tar_path) && file.info(tar_path)$size > 1024) {
  cat("    NOTE : tar file already exists (",
      round(file.info(tar_path)$size / 1024 / 1024, 1),
      " MB). Skipping download.\n", sep = "")
  download_ok <- TRUE
} else {
  download_ok <- tryCatch({
    res <- download.file(
      url = geo_url,
      destfile = tar_path,
      mode = "wb",
      quiet = FALSE
    )
    isTRUE(res == 0) || file.exists(tar_path)
  }, error = function(e) {
    cat("    download.file ERROR  : ", conditionMessage(e), "\n", sep = "")
    FALSE
  }, warning = function(w) {
    cat("    download.file WARNING: ", conditionMessage(w), "\n", sep = "")
    file.exists(tar_path) && file.info(tar_path)$size > 1024
  })

  if (!isTRUE(download_ok)) {
    cat("    Trying fallback via GEOquery::getGEOSuppFiles ...\n")
    if (requireNamespace("GEOquery", quietly = TRUE)) {
      download_ok <- tryCatch({
        GEOquery::getGEOSuppFiles(
          GEO = gse_id,
          baseDir = file.path("data", "raw"),
          makeDirectory = TRUE,
          fetch_files = TRUE
        )
        file.exists(tar_path) && file.info(tar_path)$size > 1024
      }, error = function(e) {
        cat("    GEOquery fallback failed: ", conditionMessage(e), "\n", sep = "")
        FALSE
      })
    } else {
      cat("    GEOquery is not installed. To enable fallback, run:\n")
      cat("      if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n")
      cat("      BiocManager::install('GEOquery')\n")
    }
  }
}

if (!isTRUE(download_ok) || !file.exists(tar_path)) {
  stop(
    "Failed to download GSE115440 raw archive.\n",
    "Manual URL:\n  ", geo_url, "\n",
    "Place it at:\n  ", tar_path
  )
}

tar_size_mb <- round(file.info(tar_path)$size / 1024 / 1024, 1)
cat("    OK  : downloaded/found (", tar_size_mb, " MB)\n", sep = "")

if (file.info(tar_path)$size < 1024) {
  stop("Downloaded file is smaller than 1 KB. This is likely an HTML error page.")
}
cat("\n")

# ---- 4. Extract CEL files --------------------------------------------
cat("[3/5] Extracting CEL files into ", cel_dir, " ...\n", sep = "")

extract_ok <- tryCatch({
  untar(tarfile = tar_path, exdir = cel_dir)
  TRUE
}, error = function(e) {
  cat("    untar ERROR: ", conditionMessage(e), "\n", sep = "")
  FALSE
}, warning = function(w) {
  cat("    untar WARNING: ", conditionMessage(w), "\n", sep = "")
  TRUE
})

if (!isTRUE(extract_ok)) {
  stop("Extraction failed. The tar archive may be corrupted.")
}
cat("    OK  : extraction finished\n\n")

# ---- 5. Inventory CEL files ------------------------------------------
cat("[4/5] Counting extracted CEL files ...\n")

cel_files <- list.files(
  cel_dir,
  pattern = "\\.cel(\\.gz)?$",
  ignore.case = TRUE,
  recursive = TRUE,
  full.names = TRUE
)

n_cel <- length(cel_files)
cat("    Found ", n_cel, " CEL file(s):\n", sep = "")
for (f in cel_files) {
  cat("      - ", basename(f), "\n", sep = "")
}

if (n_cel == 0L) {
  stop("No CEL files were extracted. Please verify the archive.")
}

if (n_cel != expected_n_samples) {
  cat("    WARNING: expected ", expected_n_samples,
      " CEL files but found ", n_cel, ".\n", sep = "")
  cat("    Please re-check data/metadata/GSE115440_metadata_verified.csv\n")
} else {
  cat("    OK  : count matches verified metadata (", expected_n_samples, ")\n", sep = "")
}
cat("\n")

# ---- 6. Done ----------------------------------------------------------
cat("[5/5] Step 01 complete.\n")
cat("============================================================\n")
cat(" Next suggested script:\n")
cat("   analysis/02_validation_GSE115440/scripts/02_qc_and_rma_normalize.R\n")
cat("============================================================\n")
