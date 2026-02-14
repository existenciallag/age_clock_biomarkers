###############################################################################
# run_pipeline.R — Master orchestration script
#
# Runs the full biological aging analysis pipeline:
#   1. Train all clocks on NHANES III → project to NHANES IV
#   2. Assemble unified dataset
#   3. QC diagnostics
#   4. Residualize and Z-score
#   5. Survival analysis (Cox models)
#   6. Correlation / independence analysis
#   7. Risk stratification (KM curves)
#   8. Export deployment bundle for clinical use
#
# Usage:
#   source("run_pipeline.R")
#
# To modify biomarker panels, edit R/config.R before running.
# To add a sub-clock, add an entry to PANELS_SUB in R/config.R.
###############################################################################

cat("==========================================================\n")
cat("  BIOLOGICAL AGING — MULTI-CLOCK ANALYSIS PIPELINE\n")
cat("  Based on BioAge (dayoonkwon) + PhenoAge (Levine 2018)\n")
cat("==========================================================\n\n")

# ── 0. Load dependencies ─────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(BioAge)
  library(dplyr)
  library(purrr)
  library(survival)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

# Optionally load corrplot for heatmaps
has_corrplot <- requireNamespace("corrplot", quietly = TRUE)
if (has_corrplot) library(corrplot)

# ── Load toolkit modules ──────────────────────────────────────────────────

source("R/config.R")
source("R/train.R")
source("R/assemble.R")
source("R/qc.R")
source("R/residualize.R")
source("R/survival.R")
source("R/correlation.R")
source("R/stratify.R")
source("R/export.R")

# ==========================================================================
# 1. TRAIN CLOCKS
# ==========================================================================

message("\n", strrep("=", 60))
message("STEP 1: Training biological age clocks")
message(strrep("=", 60))

clocks <- train_all_clocks(
  panels_sub = PANELS_SUB,
  train_kdm  = TRUE,
  train_hd   = TRUE,
  verbose    = TRUE
)

# ==========================================================================
# 2. ASSEMBLE DATASET
# ==========================================================================

message("\n", strrep("=", 60))
message("STEP 2: Assembling NHANES IV dataset")
message(strrep("=", 60))

data <- assemble_data(clocks)

# ==========================================================================
# 3. QUALITY CONTROL
# ==========================================================================

message("\n", strrep("=", 60))
message("STEP 3: Quality control diagnostics")
message(strrep("=", 60))

# PhenoAge vs age sanity check
qc_pheno_age_cor(data)

# QC table for all advancement + HD columns
qc <- qc_table(data)
cat("\n--- QC Summary ---\n")
print(qc, row.names = FALSE)

# Determine valid clocks
valid <- valid_clocks(qc)
message("\n>> Valid clocks for survival analysis: ", length(valid))
cat("   ", paste(valid, collapse = "\n    "), "\n")

# Overlap analysis
overlap <- qc_overlap(data)
cat("\n--- Sample Overlap ---\n")
print(overlap, row.names = FALSE)

# ==========================================================================
# 4. RESIDUALIZE & Z-SCORE
# ==========================================================================

message("\n", strrep("=", 60))
message("STEP 4: Residualization and Z-scoring")
message(strrep("=", 60))

# Add Z-scores for raw BA columns (for raw Cox models)
ba_cols <- list_ba_columns(data)
data <- add_zscores(data, ba_cols)

# Add residual + Z-scored residual for advancement columns
data <- add_residuals(data, valid)

# ==========================================================================
# 5. SURVIVAL ANALYSIS
# ==========================================================================

message("\n", strrep("=", 60))
message("STEP 5: Survival analysis (Cox proportional hazards)")
message(strrep("=", 60))

has_mortality <- all(c("time", "status") %in% names(data)) &&
  !all(is.na(data$status))

if (has_mortality) {

  # Raw BA Cox models
  cat("\n--- Raw BA models (per +1 SD) ---\n")
  hr_raw <- cox_raw(data, ba_cols, adjust_gender = FALSE)
  if (!is.null(hr_raw)) print(hr_raw, row.names = FALSE)

  # Residual age Cox models
  cat("\n--- Residual age models (per +1 SD) ---\n")
  hr_resid <- cox_residual(data, valid, adjust_gender = FALSE)
  if (!is.null(hr_resid)) print(hr_resid, row.names = FALSE)

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
  }

} else {
  message(">> No mortality data detected — skipping survival analysis")
  message("   (This is expected for clinical deployment datasets)")
  hr_all <- NULL
}

# ==========================================================================
# 6. CORRELATION & INDEPENDENCE ANALYSIS
# ==========================================================================

message("\n", strrep("=", 60))
message("STEP 6: Correlation and independence analysis")
message(strrep("=", 60))

# BA vs chronological age
cat("\n--- BA vs Chronological Age ---\n")
cor_age <- cor_ba_vs_age(data)
print(cor_age, row.names = FALSE)

# Global vs sub-clocks
cat("\n--- PhenoAge Global vs Sub-clocks ---\n")
cor_gs <- cor_global_vs_subs(data)
print(cor_gs, row.names = FALSE)

# Residualized sub-clock correlation matrix
resid_cols <- list_resid_columns(data)
if (length(resid_cols) >= 2) {
  cat("\n--- Residual sub-clock correlation ---\n")
  mat_resid <- cor_matrix(data, resid_cols)
  print(round(mat_resid, 2))

  # Heatmap
  if (has_corrplot) {
    plot_heatmap_corrplot(
      mat_resid,
      title = "Residual Sub-Clocks — Clinical Independence"
    )
  } else {
    p_heat <- plot_heatmap_gg(mat_resid, "Residual Sub-Clocks")
    print(p_heat)
  }
}

# Advancement correlation matrix
adv_cols <- list_advance_columns(data)
if (length(adv_cols) >= 2) {
  cat("\n--- Advancement (BA - CA) correlation ---\n")
  mat_adv <- cor_matrix(data, adv_cols)
  print(round(mat_adv, 2))

  if (has_corrplot) {
    plot_heatmap_corrplot(
      mat_adv,
      title = "Biological Age Advancement — System Concordance"
    )
  }
}

# Discordance analysis
cat("\n--- Pairwise system discordance ---\n")
disc <- discordance_analysis(data)
if (!is.null(disc)) print(disc, row.names = FALSE)

# ==========================================================================
# 7. RISK STRATIFICATION
# ==========================================================================

message("\n", strrep("=", 60))
message("STEP 7: Risk stratification")
message(strrep("=", 60))

if (has_mortality) {

  # Use the global PhenoAge residual for main stratification
  global_resid <- intersect(
    c("phenoage_orig_advance_resid", "phenoage_global_advance_resid"),
    names(data)
  )

  if (length(global_resid) > 0) {
    data <- add_risk_groups(data, global_resid[1])

    # Kaplan-Meier
    plot_km(data)

    # Categorical Cox
    cat("\n--- Categorical Cox (risk tertiles) ---\n")
    cat_cox <- cox_categorical(data)
    print(cat_cox)
  }

  # Multi-clock stratification
  cat("\n--- Multi-clock risk stratification ---\n")
  strat <- stratify_multi(data, resid_cols)
  if (!is.null(strat)) print(strat, row.names = FALSE)

  # Density plots of residual age
  if (length(resid_cols) > 0) {
    p_dens <- plot_resid_density(data, resid_cols)
    print(p_dens)
  }

} else {
  message(">> Skipping risk stratification (no mortality data)")
}

# ==========================================================================
# 8. EXPORT DEPLOYMENT BUNDLE
# ==========================================================================

message("\n", strrep("=", 60))
message("STEP 8: Exporting deployment bundle")
message(strrep("=", 60))

bundle <- build_deployment_bundle(
  clocks      = clocks,
  data        = data,
  clock_vars  = valid,
  hr_results  = hr_all,
  qc          = qc
)

save_bundle(bundle, "bioage_deployment_bundle.rds")

# ==========================================================================
# SUMMARY
# ==========================================================================

cat("\n")
cat("==========================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("==========================================================\n")
cat("  Subjects:      ", nrow(data), "\n")
cat("  Clocks trained:", length(bundle$meta$clocks_trained), "\n")
cat("  Valid for HR:  ", length(valid), "\n")
cat("  Bundle saved:   bioage_deployment_bundle.rds\n")
cat("==========================================================\n")

# Make key objects available in the global environment
# for interactive exploration
pipeline_results <- list(
  data    = data,
  clocks  = clocks,
  qc      = qc,
  valid   = valid,
  hr      = hr_all,
  bundle  = bundle
)

message("\n>> All results stored in 'pipeline_results' list")
message(">> Access with: pipeline_results$data, pipeline_results$hr, etc.")
