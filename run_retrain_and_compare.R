###############################################################################
# run_retrain_and_compare.R — Retrain bundle with new panels + unit fixes
#
# This script:
#   1. Retrains ALL clocks (including 4 new sub-clocks)
#   2. Rebuilds the deployment bundle
#   3. Scores the Argentine cohort with CORRECTED unit conversions
#   4. Compares old vs new results for hema_glucose
#   5. Reports diagnostics for all new clocks
#
# Usage: source("run_retrain_and_compare.R")
###############################################################################

cat("==========================================================\n")
cat("  RETRAIN + COMPARE: New panels + unit fix\n")
cat("==========================================================\n\n")

# ── 0. Setup ──────────────────────────────────────────────────────────────

PROJECT_ROOT <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) getwd())
PROJECT_ROOT <- normalizePath(PROJECT_ROOT)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(BioAge)
  library(flexsurv)
})

# Source all modules
for (mf in c("R/config.R", "R/train.R", "R/assemble.R", "R/qc.R",
             "R/residualize.R", "R/export.R")) {
  source(file.path(PROJECT_ROOT, mf), local = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 1: RETRAIN ALL CLOCKS (including new panels)
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 1: Retraining clocks with new panels")
message(strrep("=", 60))
message("Panels to train:")
for (nm in names(PANELS_SUB)) {
  message("  ", nm, ": ", paste(PANELS_SUB[[nm]], collapse = ", "))
}

clocks <- train_all_clocks(
  panels_sub = PANELS_SUB,
  train_kdm  = TRUE,
  train_hd   = TRUE,
  verbose    = TRUE
)

message("\n>> Successfully trained: ",
        paste(names(clocks$subclocks), collapse = ", "))
if (length(clocks$failed) > 0) {
  message(">> FAILED: ", paste(clocks$failed, collapse = ", "))
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 2: ASSEMBLE + BUILD BUNDLE
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 2: Assembling dataset + building bundle")
message(strrep("=", 60))

data <- assemble_data(clocks)
message(">> Assembled: ", nrow(data), " subjects, ", ncol(data), " columns")

# Identify advancement columns
adv_cols <- grep("_advance$", names(data), value = TRUE)
message(">> Advancement columns: ", length(adv_cols))

# Build bundle
bundle <- build_deployment_bundle(
  clocks     = clocks,
  data       = data,
  clock_vars = adv_cols,
  hr_results = NULL,
  qc         = NULL
)

# Save
bundle_path <- file.path(PROJECT_ROOT, "bioage_deployment_bundle.rds")
save_bundle(bundle, bundle_path)

# ══════════════════════════════════════════════════════════════════════════
# STEP 3: DIAGNOSTIC — Print model coefficients for new clocks
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 3: Model diagnostics for all sub-clocks")
message(strrep("=", 60))

for (nm in names(bundle$subclock_models)) {
  fit <- bundle$subclock_models[[nm]]
  cat(sprintf("\n=== %s ===\n", nm))
  print(fit$coef)
  cat(sprintf("  m_n=%.4f  m_d=%.6f  BA_n=%.6f  BA_d=%.6f  BA_i=%.2f\n",
              fit$m_n, fit$m_d, fit$BA_n, fit$BA_d, fit$BA_i))
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 4: SCORE ARGENTINE COHORT WITH CORRECTED UNITS
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 4: Scoring Argentine cohort (corrected units)")
message(strrep("=", 60))

# Load clinical data
f1 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE2.csv")
raw <- read.csv(f1, stringsAsFactors = FALSE)
if (file.exists(f2)) raw <- rbind(raw, read.csv(f2, stringsAsFactors = FALSE))

# Transform with CORRECT unit conversions (SI units to match BioAge models)
df <- raw %>%
  mutate(
    patient_id = Protocolo,
    age = as.numeric(age),

    # SI unit conversions (match BioAge::phenoage_nhanes internal units)
    albumin = as.numeric(albumin) * 10,               # g/dL → g/L
    alp     = as.numeric(alp),                         # U/L (no change)
    glucose = as.numeric(glucose) / 18.0,              # mg/dL → mmol/L
    crp_raw = as.numeric(crp),
    crp     = crp_raw / 10,                            # mg/L → mg/dL
    lncrp   = ifelse(is.na(crp) | crp <= 0, NA, log(crp)),
    creat   = as.numeric(creatinine) * 88.4,           # mg/dL → μmol/L
    lncreat = ifelse(is.na(creat) | creat <= 0, NA, log(creat)),
    lymph   = as.numeric(lymphocyte),
    mcv     = as.numeric(mcv),
    rdw     = as.numeric(rdw),
    wbc     = ifelse(as.numeric(wbc) > 300, as.numeric(wbc) / 1000, as.numeric(wbc)),
    rbc     = as.numeric(rbc),
    rbc     = ifelse(rbc > 100, rbc / 1e6, rbc),
    ggt     = as.numeric(ggt),
    insulin = as.numeric(insulin),
    trig    = as.numeric(triglycerides),
    totchol = as.numeric(cholesterol),
    hba1c   = as.numeric(hba1c),
    vitaminB12 = as.numeric(vitamin_b12),
    bun     = as.numeric(bun),
    uap     = as.numeric(uric_acid),
    gender  = ifelse(sex == "M", 1, 2),
    lnglucose = ifelse(is.na(glucose) | glucose <= 0, NA, log(glucose))
  )

df <- df[!is.na(df$age) & df$age >= 1 & df$age <= 120, ]
message(">> Cohort: ", nrow(df), " patients")

# ── Define all clocks to score ──
clocks_to_score <- list()

if (!is.null(bundle$models$pheno_orig)) {
  clocks_to_score$phenoage_orig <- list(
    biomarkers = BIOMARKERS_PHENO_LEVINE,
    fit        = bundle$models$pheno_orig
  )
}

for (nm in names(bundle$subclock_models)) {
  if (nm %in% names(bundle$panels$subclocks)) {
    clocks_to_score[[paste0("pheno_", nm)]] <- list(
      biomarkers = bundle$panels$subclocks[[nm]],
      fit        = bundle$subclock_models[[nm]]
    )
  }
}

# ── Score each clock ──
results <- data.frame(
  patient_id = df$patient_id,
  age        = df$age,
  sex        = df$sex,
  stringsAsFactors = FALSE
)

score_summary <- list()

for (clock_name in names(clocks_to_score)) {
  clock <- clocks_to_score[[clock_name]]
  bm    <- clock$biomarkers
  fit   <- clock$fit

  bm_plus_age <- c(bm, "age")
  missing_cols <- setdiff(bm_plus_age, names(df))

  if (length(missing_cols) > 0) {
    message(">> ", clock_name, ": SKIP — missing: ", paste(missing_cols, collapse=", "))
    next
  }

  complete_mask <- complete.cases(df[, bm_plus_age, drop = FALSE])
  n_complete <- sum(complete_mask)
  if (n_complete == 0) { message(">> ", clock_name, ": SKIP — 0 complete"); next }

  dat_subset <- df[complete_mask, bm_plus_age, drop = FALSE]

  calc_result <- tryCatch(
    phenoage_calc(data = dat_subset, biomarkers = bm, fit = fit),
    error = function(e) { message("   ERROR: ", e$message); NULL }
  )
  if (is.null(calc_result)) next

  pa_vals  <- calc_result$data$phenoage
  pa_adv   <- calc_result$data$phenoage_advance
  n_valid  <- sum(!is.na(pa_vals) & is.finite(pa_vals))

  ba_col  <- clock_name
  adv_col <- paste0(clock_name, "_advance")
  results[[ba_col]]  <- NA_real_
  results[[adv_col]] <- NA_real_
  results[[ba_col]][complete_mask]  <- pa_vals
  results[[adv_col]][complete_mask] <- pa_adv

  score_summary[[clock_name]] <- data.frame(
    clock    = clock_name,
    n_complete = n_complete,
    n_valid  = n_valid,
    pct_valid = round(100 * n_valid / n_complete, 1),
    pct_total = round(100 * n_valid / nrow(df), 1),
    stringsAsFactors = FALSE
  )

  message(sprintf(">> %-35s  %6d complete → %6d valid (%5.1f%%)",
                  clock_name, n_complete, n_valid,
                  100 * n_valid / n_complete))
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 5: COMPARISON TABLE
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 5: Results comparison")
message(strrep("=", 60))

if (length(score_summary) > 0) {
  ss <- do.call(rbind, score_summary)
  rownames(ss) <- NULL
  cat("\n--- Scoring Summary ---\n")
  print(ss, row.names = FALSE)
}

# Correlation with age
cat("\n--- Correlation with Chronological Age ---\n")
cor_results <- list()
ba_cols <- names(score_summary)

for (bc in ba_cols) {
  vals <- results[[bc]]
  ages <- results$age
  ok   <- !is.na(vals) & !is.na(ages) & is.finite(vals)
  n_ok <- sum(ok)

  if (n_ok >= 30) {
    ct <- cor.test(ages[ok], vals[ok])
    lm_fit <- lm(vals[ok] ~ ages[ok])
    s <- summary(lm_fit)

    cor_results[[bc]] <- data.frame(
      clock     = bc,
      n         = n_ok,
      r         = round(ct$estimate, 4),
      r_sq      = round(s$r.squared, 4),
      slope     = round(coef(lm_fit)[2], 4),
      intercept = round(coef(lm_fit)[1], 2),
      RMSE      = round(sqrt(mean(s$residuals^2)), 2),
      mean_adv  = round(mean(results[[paste0(bc, "_advance")]], na.rm=TRUE), 2),
      sd_adv    = round(sd(results[[paste0(bc, "_advance")]], na.rm=TRUE), 2),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %-35s r=%.4f  slope=%.3f  RMSE=%.1f  (n=%d)\n",
                bc, ct$estimate, coef(lm_fit)[2],
                sqrt(mean(s$residuals^2)), n_ok))
  }
}

if (length(cor_results) > 0) {
  cor_table <- do.call(rbind, cor_results)
  rownames(cor_table) <- NULL
  cat("\n")
  print(cor_table, row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 6: PLOTS
# ══════════════════════════════════════════════════════════════════════════

message("\n>> Generating comparison plots...")

# BA vs CA for each clock
for (bc in ba_cols) {
  vals <- results[[bc]]
  ok   <- !is.na(vals) & is.finite(vals)
  if (sum(ok) < 30) next

  r_val <- cor(results$age[ok], vals[ok])
  label <- gsub("pheno_|phenoage_", "", bc)

  p <- ggplot(results[ok, ], aes(x = age, y = .data[[bc]])) +
    geom_point(alpha = 0.08, size = 0.5, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "red", linewidth = 0.7) +
    geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("BA vs CA: ", label),
         subtitle = sprintf("r = %.3f, n = %s", r_val, format(sum(ok), big.mark=",")),
         x = "Chronological Age (years)",
         y = "Biological Age (PhenoAge)")
  print(p)
}

# ══════════════════════════════════════════════════════════════════════════
# STEP 7: SAVE
# ══════════════════════════════════════════════════════════════════════════

write.csv(results, file.path(PROJECT_ROOT, "clinical_scored_retrain.csv"),
          row.names = FALSE)
message(">> Saved: clinical_scored_retrain.csv")

if (length(cor_results) > 0) {
  write.csv(cor_table, file.path(PROJECT_ROOT, "clinical_stats_retrain.csv"),
            row.names = FALSE)
  message(">> Saved: clinical_stats_retrain.csv")
}

cat("\n==========================================================\n")
cat("  RETRAIN + COMPARE COMPLETE\n")
cat("==========================================================\n")
cat("  New panels: hema_glucose_log, metabolic_renal,\n")
cat("              hema_metabolic, pheno_max_noalb\n")
cat("  Unit fix:   glucose÷18→mmol/L, albumin×10→g/L,\n")
cat("              creatinine×88.4→μmol/L\n")
cat("==========================================================\n")
