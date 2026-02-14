###############################################################################
# run_pipeline.R — Master orchestration script
#
# Runs the COMPLETE biological aging analysis pipeline:
#   1. Train all clocks on NHANES III -> project to NHANES IV
#   2. Assemble unified dataset
#   3. QC diagnostics
#   4. Residualize and Z-score
#   5. Survival analysis (Cox models)
#   6. Correlation / independence analysis
#   7. Risk stratification (KM curves)
#   8. Export deployment bundle for clinical use
#
# Usage (from R console):
#   source("/full/path/to/run_pipeline.R")
#
# Or set working directory first:
#   setwd("/path/to/project")
#   source("run_pipeline.R")
#
# To modify biomarker panels, edit R/config.R before running.
###############################################################################

cat("==========================================================\n")
cat("  BIOLOGICAL AGING - MULTI-CLOCK ANALYSIS PIPELINE\n")
cat("  Based on BioAge (dayoonkwon) + PhenoAge (Levine 2018)\n")
cat("==========================================================\n\n")

# ======================================================================
# 0. AUTO-DETECT PROJECT ROOT AND SET WORKING DIRECTORY
# ======================================================================
# This makes source() calls work regardless of where R was launched.

PROJECT_ROOT <- tryCatch({
  # When called via source("/path/to/run_pipeline.R")
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  # Fallback: assume we are already in the project root
  getwd()
})

PROJECT_ROOT <- normalizePath(PROJECT_ROOT)

# Verify R/ folder exists
if (!dir.exists(file.path(PROJECT_ROOT, "R")) ||
    !file.exists(file.path(PROJECT_ROOT, "R", "config.R"))) {
  stop(
    "Cannot find the R/ modules folder.\n",
    "  Looked in: ", PROJECT_ROOT, "\n",
    "  Please set your working directory to the project root first:\n",
    "    setwd('/home/observer/Escritorio/R/phenoage')  # <-- your path\n",
    "    source('run_pipeline.R')"
  )
}

setwd(PROJECT_ROOT)
message(">> Project root: ", PROJECT_ROOT)

# ======================================================================
# 0b. LOAD R PACKAGES
# ======================================================================

message(">> Checking required packages...")

needed_pkgs <- c("BioAge", "dplyr", "purrr", "survival", "ggplot2",
                 "tidyr", "tibble")

missing <- needed_pkgs[!vapply(needed_pkgs, requireNamespace,
                               logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "), "\n",
    "  Install with: install.packages(c('",
    paste(setdiff(missing, "BioAge"), collapse = "','"), "'))\n",
    "  For BioAge:   devtools::install_github('dayoonkwon/BioAge')"
  )
}

suppressPackageStartupMessages({
  library(BioAge)
  library(dplyr)
  library(purrr)
  library(survival)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

has_corrplot <- requireNamespace("corrplot", quietly = TRUE)
if (has_corrplot) library(corrplot)
message(">> All packages loaded OK")

# ======================================================================
# 0c. SOURCE ALL TOOLKIT MODULES
# ======================================================================

message(">> Loading toolkit modules...")

module_files <- c(
  "R/config.R",
  "R/train.R",
  "R/assemble.R",
  "R/qc.R",
  "R/residualize.R",
  "R/survival.R",
  "R/correlation.R",
  "R/stratify.R",
  "R/export.R"
)

for (mf in module_files) {
  fp <- file.path(PROJECT_ROOT, mf)
  if (!file.exists(fp)) {
    stop("Module not found: ", fp)
  }
  source(fp, local = FALSE)
}

message(">> All modules loaded OK\n")

# ======================================================================
# STEP 1. TRAIN CLOCKS
# ======================================================================

message(strrep("=", 60))
message("STEP 1: Training biological age clocks")
message(strrep("=", 60))

clocks <- train_all_clocks(
  panels_sub = PANELS_SUB,
  train_kdm  = TRUE,
  train_hd   = TRUE,
  verbose    = TRUE
)

# ======================================================================
# STEP 2. ASSEMBLE DATASET
# ======================================================================

message("\n", strrep("=", 60))
message("STEP 2: Assembling NHANES IV dataset")
message(strrep("=", 60))

data <- assemble_data(clocks)

# ======================================================================
# STEP 3. QUALITY CONTROL
# ======================================================================

message("\n", strrep("=", 60))
message("STEP 3: Quality control diagnostics")
message(strrep("=", 60))

# PhenoAge vs age sanity check (expect ~0.94)
qc_pheno_age_cor(data)

# QC table for all advancement + HD columns
qc <- qc_table(data)
cat("\n--- QC Summary ---\n")
print(qc, row.names = FALSE)

# Determine valid clocks (enough data + non-zero variance)
valid <- valid_clocks(qc)
message("\n>> Valid clocks for survival analysis: ", length(valid))
for (v in valid) cat("    - ", v, "\n")

# Overlap analysis
overlap <- qc_overlap(data)
cat("\n--- Sample Overlap ---\n")
print(overlap, row.names = FALSE)

# ======================================================================
# STEP 4. RESIDUALIZE & Z-SCORE
# ======================================================================

message("\n", strrep("=", 60))
message("STEP 4: Residualization and Z-scoring")
message(strrep("=", 60))

# Z-scores on raw BA columns (for raw Cox models)
ba_cols <- list_ba_columns(data)
message(">> BA columns: ", paste(ba_cols, collapse = ", "))
data <- add_zscores(data, ba_cols)

# Residualize advancement columns against age, then Z-score
data <- add_residuals(data, valid)

# ======================================================================
# STEP 5. SURVIVAL ANALYSIS
# ======================================================================

message("\n", strrep("=", 60))
message("STEP 5: Survival analysis (Cox proportional hazards)")
message(strrep("=", 60))

has_mortality <- all(c("time", "status") %in% names(data)) &&
  !all(is.na(data$status))

hr_all <- NULL

if (has_mortality) {

  message(">> Mortality data detected: ",
          sum(data$status == 1, na.rm = TRUE), " deaths / ",
          sum(!is.na(data$status)), " subjects")

  # Raw BA Cox models
  cat("\n--- Raw BA models (Surv ~ age + z_clock, per +1 SD) ---\n")
  hr_raw <- cox_raw(data, ba_cols, adjust_gender = FALSE)
  if (!is.null(hr_raw) && nrow(hr_raw) > 0) {
    print(hr_raw, row.names = FALSE)
  } else {
    message("   (no valid raw models)")
  }

  # Residual age Cox models
  cat("\n--- Residual age models (Surv ~ age + resid_z, per +1 SD) ---\n")
  hr_resid <- cox_residual(data, valid, adjust_gender = FALSE)
  if (!is.null(hr_resid) && nrow(hr_resid) > 0) {
    print(hr_resid, row.names = FALSE)
  } else {
    message("   (no valid residual models)")
  }

  # Combined table
  hr_all <- hr_table(data, valid, adjust_gender = FALSE)

  # Forest plot
  if (!is.null(hr_resid) && nrow(hr_resid) > 0) {
    p_forest <- forest_plot(
      hr_resid,
      title    = "Residual Biological Age and Mortality Risk",
      subtitle = "Hazard Ratio per +1 SD (adjusted for chronological age)"
    )
    print(p_forest)
    message(">> Forest plot displayed")
  }

} else {
  message(">> No mortality data detected - skipping survival analysis")
  message("   (This is expected for clinical deployment datasets)")
}

# ======================================================================
# STEP 6. CORRELATION & INDEPENDENCE ANALYSIS
# ======================================================================

message("\n", strrep("=", 60))
message("STEP 6: Correlation and independence analysis")
message(strrep("=", 60))

# BA vs chronological age
cat("\n--- BA vs Chronological Age ---\n")
cor_age_result <- cor_ba_vs_age(data)
print(cor_age_result, row.names = FALSE)

# Global vs sub-clocks
cat("\n--- PhenoAge Global vs Sub-clocks ---\n")
cor_gs <- cor_global_vs_subs(data)
print(cor_gs, row.names = FALSE)

# Residualized sub-clock correlation matrix + heatmap
resid_cols <- list_resid_columns(data)
if (length(resid_cols) >= 2) {
  cat("\n--- Residual sub-clock correlation matrix ---\n")
  mat_resid <- cor_matrix(data, resid_cols)
  print(round(mat_resid, 2))

  if (has_corrplot) {
    plot_heatmap_corrplot(mat_resid,
                          title = "Residual Sub-Clocks - Clinical Independence")
  } else {
    print(plot_heatmap_gg(mat_resid, "Residual Sub-Clocks"))
  }
}

# Advancement correlation matrix + heatmap
adv_cols <- list_advance_columns(data)
if (length(adv_cols) >= 2) {
  cat("\n--- Advancement (BA - CA) correlation matrix ---\n")
  mat_adv <- cor_matrix(data, adv_cols)
  print(round(mat_adv, 2))

  if (has_corrplot) {
    plot_heatmap_corrplot(mat_adv,
                          title = "Biological Age Advancement - Concordance")
  }
}

# Discordance analysis
cat("\n--- Pairwise system discordance ---\n")
disc <- discordance_analysis(data)
if (!is.null(disc)) print(disc, row.names = FALSE)

# ======================================================================
# STEP 7. RISK STRATIFICATION
# ======================================================================

message("\n", strrep("=", 60))
message("STEP 7: Risk stratification")
message(strrep("=", 60))

if (has_mortality) {

  global_resid <- intersect(
    c("phenoage_orig_advance_resid", "phenoage_global_advance_resid"),
    names(data)
  )

  if (length(global_resid) > 0) {
    message(">> Stratifying on: ", global_resid[1])
    data <- add_risk_groups(data, global_resid[1])

    # Kaplan-Meier
    plot_km(data)

    # Categorical Cox
    cat("\n--- Categorical Cox (risk tertiles vs mortality) ---\n")
    cat_cox <- cox_categorical(data)
    print(cat_cox)
  }

  # Multi-clock stratification
  cat("\n--- Multi-clock risk stratification ---\n")
  strat <- stratify_multi(data, resid_cols)
  if (!is.null(strat) && nrow(strat) > 0) print(strat, row.names = FALSE)

  # Density plots
  if (length(resid_cols) > 0) {
    print(plot_resid_density(data, resid_cols))
  }

} else {
  message(">> Skipping risk stratification (no mortality data)")
}

# ======================================================================
# STEP 8. EXPORT DEPLOYMENT BUNDLE
# ======================================================================

message("\n", strrep("=", 60))
message("STEP 8: Exporting deployment bundle for clinical use")
message(strrep("=", 60))

bundle <- build_deployment_bundle(
  clocks      = clocks,
  data        = data,
  clock_vars  = valid,
  hr_results  = hr_all,
  qc          = qc
)

save_bundle(bundle, file.path(PROJECT_ROOT, "bioage_deployment_bundle.rds"))

# ======================================================================
# FINAL SUMMARY
# ======================================================================

cat("\n")
cat("==========================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("==========================================================\n")
cat("  Project root:   ", PROJECT_ROOT, "\n")
cat("  Subjects:       ", nrow(data), "\n")
cat("  Clocks trained: ", length(bundle$meta$clocks_trained), "\n")
cat("  Valid for HR:   ", length(valid), "\n")
if (!is.null(hr_all)) {
  cat("  HR models:      ", nrow(hr_all), "\n")
}
cat("  Bundle saved:    bioage_deployment_bundle.rds\n")
cat("==========================================================\n")

# Store everything in a single list for interactive use
pipeline_results <- list(
  data    = data,
  clocks  = clocks,
  qc      = qc,
  valid   = valid,
  hr      = hr_all,
  bundle  = bundle
)

message("\n>> All results stored in 'pipeline_results'")
message(">> Access: pipeline_results$data, pipeline_results$hr, etc.")
message(">> To score new patients later: source('run_clinical.R')")
