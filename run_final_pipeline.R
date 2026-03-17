###############################################################################
# run_final_pipeline.R — Complete Biological Aging Analysis (Paper-Ready)
#
# UNIFIED PIPELINE: Train → Score → Recenter → Discovery
#
# This script runs the entire analysis in a single source() call:
#   STEP 1. Train clocks on NHANES III → validate on NHANES IV
#   STEP 2. Score Argentine clinical cohort (131K patients)
#   STEP 3. Population recentering (3 methods)
#   STEP 4. Fast/slow ager profiling & biomarker signature
#
# FOCUS CLOCKS (by Argentine data coverage):
#   - pheno_hema_integrated          (85% coverage)
#   - pheno_hema_glucose             (67% coverage)
#   - pheno_pheno_max_v2             (57% coverage — 6/9 Levine)
#   - pheno_micronutrient_methylation (1.3% coverage — low but informative)
#
# USAGE:  source("run_final_pipeline.R")
###############################################################################

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  BIOLOGICAL AGING — COMPLETE ANALYSIS PIPELINE          ║\n")
cat("║  Train (NHANES) → Score (Argentina) → Recenter → Profile║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# ══════════════════════════════════════════════════════════════════════════
# 0. SETUP
# ══════════════════════════════════════════════════════════════════════════

PROJECT_ROOT <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) getwd())
PROJECT_ROOT <- normalizePath(PROJECT_ROOT)
setwd(PROJECT_ROOT)
message(">> Project root: ", PROJECT_ROOT)

suppressPackageStartupMessages({
  library(BioAge)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(flexsurv)
})

has_corrplot <- requireNamespace("corrplot", quietly = TRUE)
if (has_corrplot) library(corrplot)

# Source modules
for (mf in c("R/config.R", "R/train.R", "R/assemble.R", "R/residualize.R",
             "R/export.R", "R/qc_recenter.R")) {
  source(file.path(PROJECT_ROOT, mf), local = FALSE)
}
message(">> All modules loaded")

# Output directory
OUTPUT_DIR <- file.path(PROJECT_ROOT, "final_output")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Focus clocks — the ones with meaningful Argentine coverage
FOCUS_CLOCKS <- c("hema_integrated", "hema_glucose", "pheno_max_v2",
                  "micronutrient_methylation")

# ══════════════════════════════════════════════════════════════════════════
# STEP 1. TRAIN CLOCKS ON NHANES III → VALIDATE ON NHANES IV
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 1: TRAIN CLOCKS ON NHANES                        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

message(">> Training all clocks (NHANES III → NHANES IV)...")

clocks <- train_all_clocks(
  panels_sub = PANELS_SUB,
  train_kdm  = TRUE,
  train_hd   = TRUE,
  verbose    = TRUE
)

# Assemble NHANES IV dataset
data_nhanes <- assemble_data(clocks)

# ── 1a. NHANES IV Performance Report ──

cat("\n--- NHANES IV Clock Performance ---\n\n")

nhanes_perf <- list()

# BA columns for all trained clocks
ba_cols_nhanes <- list_ba_columns(data_nhanes)

for (bc in ba_cols_nhanes) {
  vals <- data_nhanes[[bc]]
  ages <- data_nhanes$age
  ok   <- !is.na(vals) & !is.na(ages) & is.finite(vals)
  n_ok <- sum(ok)
  if (n_ok < 30) next

  ct     <- cor.test(ages[ok], vals[ok])
  lm_fit <- lm(vals[ok] ~ ages[ok])
  s      <- summary(lm_fit)

  nhanes_perf[[bc]] <- data.frame(
    clock     = bc,
    n         = n_ok,
    r         = round(ct$estimate, 4),
    r_sq      = round(s$r.squared, 4),
    slope     = round(coef(lm_fit)[2], 4),
    intercept = round(coef(lm_fit)[1], 2),
    RMSE      = round(sqrt(mean(s$residuals^2)), 2),
    mean_adv  = round(mean(vals[ok] - ages[ok]), 2),
    sd_adv    = round(sd(vals[ok] - ages[ok]), 2),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

if (length(nhanes_perf) > 0) {
  nhanes_table <- do.call(rbind, nhanes_perf)
  rownames(nhanes_table) <- NULL
  print(nhanes_table, row.names = FALSE)
  write.csv(nhanes_table, file.path(OUTPUT_DIR, "step1_nhanes_performance.csv"),
            row.names = FALSE)
  message(">> Saved: step1_nhanes_performance.csv")
}

# ── 1b. NHANES IV scatter plots (BA vs CA) for focus clocks ──

cat("\n--- NHANES IV: BA vs CA Scatter Plots ---\n\n")

focus_ba <- paste0("pheno_", FOCUS_CLOCKS)
focus_ba <- intersect(focus_ba, ba_cols_nhanes)

for (bc in focus_ba) {
  vals <- data_nhanes[[bc]]
  ok   <- !is.na(vals) & is.finite(vals)
  if (sum(ok) < 30) next

  r_val <- cor(data_nhanes$age[ok], vals[ok])
  label <- gsub("pheno_", "", bc)

  p <- ggplot(data_nhanes[ok, ], aes(x = age, y = .data[[bc]])) +
    geom_point(alpha = 0.08, size = 0.5, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "red", linewidth = 0.7) +
    geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("NHANES IV: BA vs CA — ", label),
         subtitle = sprintf("r = %.4f, n = %s", r_val, format(sum(ok), big.mark = ",")),
         x = "Chronological Age (years)",
         y = "Biological Age (PhenoAge)")
  print(p)

  ggsave(file.path(OUTPUT_DIR, paste0("step1_nhanes_scatter_", label, ".png")),
         p, width = 8, height = 6, dpi = 150)
}

# ── 1c. Save deployment bundle ──

# Z-score and residualize for bundle export
ba_all <- list_ba_columns(data_nhanes)
data_nhanes <- add_zscores(data_nhanes, ba_all)
adv_cols_nhanes <- list_advance_columns(data_nhanes)
data_nhanes <- add_residuals(data_nhanes, adv_cols_nhanes)

bundle <- build_deployment_bundle(
  clocks     = clocks,
  data       = data_nhanes,
  clock_vars = adv_cols_nhanes,
  hr_results = NULL,
  qc         = NULL
)
save_bundle(bundle, file.path(PROJECT_ROOT, "bioage_deployment_bundle.rds"))

message("\n>> STEP 1 COMPLETE: ", length(bundle$meta$clocks_trained), " clocks trained")


# ══════════════════════════════════════════════════════════════════════════
# STEP 2. SCORE ARGENTINE CLINICAL COHORT
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 2: SCORE ARGENTINE COHORT                        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

message(">> Loading clinical data...")
f1 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE2.csv")
cohort <- read.csv(f1, stringsAsFactors = FALSE)
if (file.exists(f2)) cohort <- rbind(cohort, read.csv(f2, stringsAsFactors = FALSE))
message(">> Rows: ", nrow(cohort))

# ── 2a. Transform biomarkers ──

df <- cohort %>%
  mutate(
    patient_id = Protocolo,
    age = as.numeric(age),
    albumin = as.numeric(albumin) * 10,         # g/dL → g/L
    alp     = as.numeric(alp),
    glucose = as.numeric(glucose) / 18.0,        # mg/dL → mmol/L
    crp_raw = as.numeric(crp),
    crp     = crp_raw / 10,                     # mg/L → mg/dL
    lncrp   = ifelse(is.na(crp) | crp <= 0, NA, log(crp)),
    creat   = as.numeric(creatinine) * 88.4,    # mg/dL → μmol/L
    lncreat = ifelse(is.na(creat) | creat <= 0, NA, log(creat)),
    lymph   = as.numeric(lymphocyte),
    mcv     = as.numeric(mcv),
    rdw     = as.numeric(rdw),
    wbc     = ifelse(as.numeric(wbc) > 300, as.numeric(wbc) / 1000, as.numeric(wbc)),
    rbc     = as.numeric(rbc),
    rbc     = ifelse(rbc > 100, rbc / 1e6, rbc),
    ggt     = as.numeric(ggt),
    insulin    = as.numeric(insulin),
    trig       = as.numeric(triglycerides),
    totchol    = as.numeric(cholesterol),
    hba1c      = as.numeric(hba1c),
    vitaminB12 = as.numeric(vitamin_b12),
    bun        = as.numeric(bun),
    uap        = as.numeric(uric_acid),
    gender     = ifelse(sex == "M", 1, 2),
    lnglucose  = ifelse(is.na(glucose) | glucose <= 0, NA, log(glucose))
  )

df <- df[!is.na(df$age) & df$age >= 18 & df$age <= 120, ]
message(">> After age filter: ", nrow(df), " rows")
message(">> Age: ", round(mean(df$age, na.rm = TRUE), 1), " ± ",
        round(sd(df$age, na.rm = TRUE), 1), " years")

# ── 2b. Biomarker availability ──

cat("\n--- Biomarker Availability ---\n\n")

all_bm <- c("albumin", "alp", "glucose", "lncrp", "lncreat",
            "lymph", "mcv", "rdw", "wbc", "rbc", "ggt", "insulin",
            "trig", "totchol", "vitaminB12", "hba1c", "bun", "uap", "lnglucose")
for (v in all_bm) {
  if (v %in% names(df)) {
    n_ok <- sum(!is.na(df[[v]]))
    cat(sprintf("  %-14s %6d / %d  (%5.1f%%)\n",
                v, n_ok, nrow(df), 100 * n_ok / nrow(df)))
  }
}

# ── 2c. Score focus clocks ──

cat("\n--- Scoring Focus Clocks ---\n\n")

clocks_to_score <- list()

# Add focus sub-clocks from bundle
for (nm in names(bundle$subclock_models)) {
  if (nm %in% names(bundle$panels$subclocks)) {
    clocks_to_score[[paste0("pheno_", nm)]] <- list(
      biomarkers = bundle$panels$subclocks[[nm]],
      fit        = bundle$subclock_models[[nm]]
    )
  }
}

# Master results
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
  if (length(missing_cols) > 0) next

  complete_mask <- complete.cases(df[, bm_plus_age, drop = FALSE])
  n_complete <- sum(complete_mask)
  if (n_complete == 0) {
    message(">> ", clock_name, ": SKIP — 0 complete cases")
    next
  }

  message(">> ", clock_name, ": ", n_complete, " complete cases (",
          round(100 * n_complete / nrow(df), 1), "%)")

  dat_subset <- df[complete_mask, bm_plus_age, drop = FALSE]

  calc_result <- tryCatch({
    phenoage_calc(data = dat_subset, biomarkers = bm, fit = fit)
  }, error = function(e) {
    message("   ERROR: ", e$message)
    NULL
  })
  if (is.null(calc_result)) next

  pa_vals <- calc_result$data$phenoage
  n_valid <- sum(!is.na(pa_vals) & is.finite(pa_vals))
  message("   -> ", n_valid, " valid PhenoAge values")

  if (n_valid > 0) {
    ba_col  <- clock_name
    adv_col <- paste0(clock_name, "_advance")
    results[[ba_col]]  <- NA_real_
    results[[adv_col]] <- NA_real_
    results[[ba_col]][complete_mask]  <- pa_vals
    results[[adv_col]][complete_mask] <- calc_result$data$phenoage_advance

    score_summary[[clock_name]] <- data.frame(
      clock    = clock_name,
      n_scored = n_valid,
      pct      = round(100 * n_valid / nrow(df), 1),
      stringsAsFactors = FALSE
    )
  }
}

# ── 2d. Argentine scoring performance ──

cat("\n--- Argentine Cohort: Clock Performance ---\n\n")

ba_cols_arg <- names(score_summary)
arg_perf <- list()

for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ages <- results$age
  ok   <- !is.na(vals) & !is.na(ages) & is.finite(vals)
  n_ok <- sum(ok)
  if (n_ok < 10) next

  ct     <- cor.test(ages[ok], vals[ok])
  lm_fit <- lm(vals[ok] ~ ages[ok])
  s      <- summary(lm_fit)

  arg_perf[[bc]] <- data.frame(
    clock     = bc,
    n         = n_ok,
    r         = round(ct$estimate, 4),
    r_sq      = round(s$r.squared, 4),
    slope     = round(coef(lm_fit)[2], 4),
    intercept = round(coef(lm_fit)[1], 2),
    RMSE      = round(sqrt(mean(s$residuals^2)), 2),
    mean_adv  = round(mean(vals[ok] - ages[ok]), 2),
    sd_adv    = round(sd(vals[ok] - ages[ok]), 2),
    stringsAsFactors = FALSE, row.names = NULL
  )

  cat(sprintf("  %-40s r = %.4f  RMSE = %.2f  n = %s\n",
              bc, ct$estimate, sqrt(mean(s$residuals^2)),
              format(n_ok, big.mark = ",")))
}

if (length(arg_perf) > 0) {
  arg_table <- do.call(rbind, arg_perf)
  rownames(arg_table) <- NULL
  cat("\n")
  print(arg_table, row.names = FALSE)
  write.csv(arg_table, file.path(OUTPUT_DIR, "step2_argentina_performance.csv"),
            row.names = FALSE)
  message(">> Saved: step2_argentina_performance.csv")
}

# ── 2d-bis. RMSE by Age Group (Argentine Cohort) ──

cat("\n--- RMSE by Age Group (Argentine Cohort) ---\n\n")

age_breaks  <- c(0, 30, 40, 50, 60, 70, 80, 120)
age_labels  <- c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
results$age_group_rmse <- cut(results$age, breaks = age_breaks, labels = age_labels,
                              right = FALSE)

rmse_by_age <- list()
for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ages <- results$age
  grp  <- results$age_group_rmse
  for (ag in age_labels) {
    ok <- !is.na(vals) & is.finite(vals) & !is.na(ages) & !is.na(grp) & grp == ag
    n_ok <- sum(ok)
    if (n_ok < 10) next
    resid <- vals[ok] - ages[ok]
    rmse_val  <- sqrt(mean(resid^2))
    mae_val   <- mean(abs(resid))
    r_val     <- cor(ages[ok], vals[ok])
    bias_val  <- mean(resid)
    rmse_by_age[[paste0(bc, "_", ag)]] <- data.frame(
      clock     = bc,
      age_group = ag,
      n         = n_ok,
      RMSE      = round(rmse_val, 2),
      MAE       = round(mae_val, 2),
      bias      = round(bias_val, 2),
      r         = round(r_val, 4),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }
}

if (length(rmse_by_age) > 0) {
  rmse_age_table <- do.call(rbind, rmse_by_age)
  rownames(rmse_age_table) <- NULL

  # Print per clock
  for (bc in ba_cols_arg) {
    sub_t <- rmse_age_table[rmse_age_table$clock == bc, ]
    if (nrow(sub_t) == 0) next
    cat(sprintf("\n  %s:\n", bc))
    cat(sprintf("    %-8s %6s %7s %7s %7s %7s\n", "Age", "n", "RMSE", "MAE", "Bias", "r"))
    for (i in seq_len(nrow(sub_t))) {
      r <- sub_t[i, ]
      cat(sprintf("    %-8s %6d %7.2f %7.2f %7.2f %7.4f\n",
                  r$age_group, r$n, r$RMSE, r$MAE, r$bias, r$r))
    }
    # Highlight worst RMSE group
    worst <- sub_t[which.max(sub_t$RMSE), ]
    best  <- sub_t[which.min(sub_t$RMSE), ]
    cat(sprintf("    >> Worst RMSE: %s (%.2f)  |  Best RMSE: %s (%.2f)\n",
                worst$age_group, worst$RMSE, best$age_group, best$RMSE))
  }

  write.csv(rmse_age_table, file.path(OUTPUT_DIR, "step2_rmse_by_age_group.csv"),
            row.names = FALSE)
  message(">> Saved: step2_rmse_by_age_group.csv")

  # Plot: RMSE by age group per clock
  for (bc in unique(rmse_age_table$clock)) {
    sub_t <- rmse_age_table[rmse_age_table$clock == bc, ]
    if (nrow(sub_t) < 3) next
    label <- gsub("pheno_", "", bc)
    sub_t$age_group <- factor(sub_t$age_group, levels = age_labels)

    p_rmse <- ggplot(sub_t, aes(x = age_group, y = RMSE)) +
      geom_col(fill = "steelblue", alpha = 0.8, width = 0.6) +
      geom_text(aes(label = sprintf("%.1f", RMSE)), vjust = -0.5, size = 3.5) +
      geom_line(aes(x = as.numeric(age_group), y = RMSE),
                colour = "darkred", linewidth = 0.8, linetype = "dashed") +
      theme_minimal(base_size = 13) +
      labs(title = paste0("RMSE by Age Group — ", label),
           subtitle = sprintf("n per group shown below bars"),
           x = "Age Group", y = "RMSE (years)") +
      geom_text(aes(label = paste0("n=", format(n, big.mark = ",")), y = 0),
                vjust = 1.5, size = 2.8, colour = "grey40")
    print(p_rmse)

    fn <- file.path(OUTPUT_DIR, paste0("step2_rmse_by_age_", label, ".png"))
    ggsave(fn, p_rmse, width = 8, height = 5, dpi = 150)
  }
  message(">> Saved: RMSE by age group plots")
}

# ── 2e. Argentine scatter plots (BA vs CA) ──

for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ok   <- !is.na(vals) & is.finite(vals)
  if (sum(ok) < 30) next

  r_val <- cor(results$age[ok], vals[ok])
  label <- gsub("pheno_", "", bc)

  p <- ggplot(results[ok, ], aes(x = age, y = .data[[bc]])) +
    geom_point(alpha = 0.05, size = 0.3, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "red", linewidth = 0.7) +
    geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Argentina: BA vs CA — ", label),
         subtitle = sprintf("r = %.4f, n = %s", r_val, format(sum(ok), big.mark = ",")),
         x = "Chronological Age (years)",
         y = "Biological Age (PhenoAge)")
  print(p)

  ggsave(file.path(OUTPUT_DIR, paste0("step2_arg_scatter_", label, ".png")),
         p, width = 8, height = 6, dpi = 150)
}

# ── 2f. BA vs CA by sex ──

for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ok <- !is.na(vals) & is.finite(vals) & results$sex %in% c("F", "M")
  if (sum(ok) < 50) next
  label <- gsub("pheno_", "", bc)

  p_sex <- ggplot(results[ok, ], aes(x = age, y = .data[[bc]], colour = sex)) +
    geom_point(alpha = 0.03, size = 0.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray40") +
    scale_colour_manual(values = c("F" = "firebrick", "M" = "steelblue")) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Argentina — BA vs CA by Sex: ", label),
         x = "Chronological Age", y = "Biological Age", colour = "Sex")
  print(p_sex)

  ggsave(file.path(OUTPUT_DIR, paste0("step2_arg_sex_", label, ".png")),
         p_sex, width = 8, height = 6, dpi = 150)
}

# Save scored data
write.csv(results, file.path(PROJECT_ROOT, "clinical_scored_lite.csv"),
          row.names = FALSE)
message(">> Saved: clinical_scored_lite.csv")

message("\n>> STEP 2 COMPLETE: ", length(score_summary), " clocks scored on ",
        format(nrow(results), big.mark = ","), " patients")


# ══════════════════════════════════════════════════════════════════════════
# STEP 3. POPULATION RECENTERING
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 3: POPULATION RECENTERING                        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# Work on a copy of results for recentering
rc <- results

# Identify advancement columns with enough data
adv_cols <- grep("_advance$", names(rc), value = TRUE)
viable_adv <- character(0)
for (ac in adv_cols) {
  n_ok <- sum(!is.na(rc[[ac]]) & is.finite(rc[[ac]]))
  if (n_ok >= 100) viable_adv <- c(viable_adv, ac)
}
message(">> Viable clocks for recentering: ", paste(sub("_advance$", "", viable_adv),
                                                     collapse = ", "))

# ── 3a. Method A: Within-cohort residualization (BA ~ age) ──

cat("\n--- Method A: Within-Cohort Residualization (BA ~ age) ---\n\n")

for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  ba_vals <- rc[[ba_col]]
  ok <- !is.na(ba_vals) & is.finite(ba_vals) & !is.na(rc$age)
  if (sum(ok) < 30) next

  fit_local <- lm(ba_vals[ok] ~ rc$age[ok])
  coefs <- coef(fit_local)

  resid_col <- paste0(ba_col, "_local_resid")
  z_col     <- paste0(ba_col, "_local_z")

  rc[[resid_col]] <- NA_real_
  rc[[resid_col]][ok] <- residuals(fit_local)
  rc[[z_col]] <- NA_real_
  valid_resid <- !is.na(rc[[resid_col]])
  if (sum(valid_resid) > 1)
    rc[[z_col]][valid_resid] <- as.numeric(scale(rc[[resid_col]][valid_resid]))

  cat(sprintf("  %-40s intercept = %+.3f  slope = %.4f  n = %d\n",
              ba_col, coefs[1], coefs[2], sum(ok)))
}

# ── 3b. Method B: Sex-adjusted residualization (BA ~ age + sex) ──

cat("\n--- Method B: Sex-Adjusted Residualization (BA ~ age + sex) ---\n\n")

for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  ba_vals <- rc[[ba_col]]
  ok <- !is.na(ba_vals) & is.finite(ba_vals) &
        !is.na(rc$age) & rc$sex %in% c("F", "M")
  if (sum(ok) < 30) next

  fit_sex <- lm(ba_vals[ok] ~ rc$age[ok] + rc$sex[ok])
  coefs <- coef(fit_sex)

  resid_col <- paste0(ba_col, "_sexadj_resid")
  z_col     <- paste0(ba_col, "_sexadj_z")

  rc[[resid_col]] <- NA_real_
  rc[[resid_col]][ok] <- residuals(fit_sex)
  rc[[z_col]] <- NA_real_
  valid_resid <- !is.na(rc[[resid_col]])
  if (sum(valid_resid) > 1)
    rc[[z_col]][valid_resid] <- as.numeric(scale(rc[[resid_col]][valid_resid]))

  cat(sprintf("  %-40s slope_age = %.4f  slope_sexM = %+.3f  n = %d\n",
              ba_col, coefs[2], ifelse(length(coefs) >= 3, coefs[3], NA), sum(ok)))
}

# ── 3c. Method C: Age-stratified Z-scoring ──

cat("\n--- Method C: Age-Stratified Z-Scoring ---\n\n")

rc$age_stratum <- cut(rc$age,
  breaks = c(0, 30, 40, 50, 60, 70, 80, 120),
  labels = c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
  right = FALSE
)

for (ac in viable_adv) {
  ba_col  <- sub("_advance$", "", ac)
  adv_vals <- rc[[ac]]
  strat_z_col <- paste0(ba_col, "_strat_z")
  rc[[strat_z_col]] <- NA_real_

  for (stratum in levels(rc$age_stratum)) {
    idx <- rc$age_stratum == stratum & !is.na(adv_vals) & is.finite(adv_vals)
    n_s <- sum(idx, na.rm = TRUE)
    if (n_s < 10) next
    x_s <- adv_vals[idx]
    mu_s <- mean(x_s); sd_s <- sd(x_s)
    if (sd_s > 0) rc[[strat_z_col]][idx] <- (x_s - mu_s) / sd_s
  }
}

# ── 3d. Post-recentering QC ──

cat("\n--- Post-Recentering QC (Method B — sex-adjusted, RECOMMENDED) ---\n\n")

sexadj_z_cols <- paste0(sub("_advance$", "", viable_adv), "_sexadj_z")
sexadj_z_cols <- intersect(sexadj_z_cols, names(rc))

if (length(sexadj_z_cols) > 0) {
  post_qc_stats <- lapply(sexadj_z_cols, function(zc) {
    x <- rc[[zc]]
    ok <- !is.na(x) & is.finite(x)
    clock <- sub("_sexadj_z$", "", zc)
    data.frame(
      clock      = clock,
      n          = sum(ok),
      mean       = round(mean(x[ok]), 4),
      sd         = round(sd(x[ok]), 4),
      median     = round(median(x[ok]), 4),
      pct_above0 = round(100 * mean(x[ok] > 0), 1),
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  post_qc <- do.call(rbind, post_qc_stats)
  cat("  (Expected: mean ~ 0, SD ~ 1, pct_above0 ~ 50%)\n\n")
  print(post_qc, row.names = FALSE)
  write.csv(post_qc, file.path(OUTPUT_DIR, "step3_recenter_qc.csv"), row.names = FALSE)
}

# ── 3e. Density plots: before vs after recentering ──

for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  local_col <- paste0(ba_col, "_local_resid")
  if (!local_col %in% names(rc)) next

  p <- plot_recenter_density(rc, ac, local_col)
  label <- gsub("pheno_|phenoage_", "", ba_col)
  ggsave(file.path(OUTPUT_DIR, paste0("step3_recenter_density_", label, ".png")),
         p, width = 9, height = 5, dpi = 150)
  print(p)
}

# ── 3f. Inter-clock correlation (recentered) ──

if (length(sexadj_z_cols) >= 2) {
  cat("\n--- Inter-Clock Spearman Correlation (Recentered) ---\n\n")

  cor_data <- rc[, sexadj_z_cols, drop = FALSE]
  for (v in sexadj_z_cols) cor_data[[v]] <- as.numeric(cor_data[[v]])
  cor_mat <- cor(cor_data, use = "pairwise.complete.obs", method = "spearman")
  colnames(cor_mat) <- gsub("_sexadj_z$", "", colnames(cor_mat))
  rownames(cor_mat) <- gsub("_sexadj_z$", "", rownames(cor_mat))

  print(round(cor_mat, 3))
  write.csv(as.data.frame(cor_mat),
            file.path(OUTPUT_DIR, "step3_interclock_correlation.csv"))

  if (has_corrplot) {
    png(file.path(OUTPUT_DIR, "step3_interclock_heatmap.png"),
        width = 800, height = 700, res = 120)
    corrplot::corrplot(cor_mat, method = "color", type = "upper",
                       tl.col = "black", tl.cex = 0.9,
                       addCoef.col = "black", number.cex = 0.8,
                       title = "Inter-Clock Correlation (Recentered, Sex-Adjusted)",
                       mar = c(0, 0, 2, 0))
    dev.off()
    message(">> Saved: step3_interclock_heatmap.png")
  }
}

# Save recentered dataset
out_cols <- c("patient_id", "age", "sex")
for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  out_cols <- c(out_cols, ba_col, ac,
    paste0(ba_col, "_local_resid"), paste0(ba_col, "_local_z"),
    paste0(ba_col, "_sexadj_resid"), paste0(ba_col, "_sexadj_z"),
    paste0(ba_col, "_strat_z"))
}
out_cols <- intersect(out_cols, names(rc))
write.csv(rc[, out_cols], file.path(OUTPUT_DIR, "clinical_recentered.csv"),
          row.names = FALSE)
# Also save to standard location for run_subclock_discovery compatibility
recenter_out <- file.path(PROJECT_ROOT, "recenter_output")
if (!dir.exists(recenter_out)) dir.create(recenter_out, recursive = TRUE)
write.csv(rc[, out_cols], file.path(recenter_out, "clinical_recentered.csv"),
          row.names = FALSE)

message("\n>> STEP 3 COMPLETE: ", length(viable_adv), " clocks recentered")


# ══════════════════════════════════════════════════════════════════════════
# STEP 4. FAST/SLOW AGER PROFILING & BIOMARKER SIGNATURES
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 4: FAST vs SLOW AGER PROFILING                   ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# Use sex-adjusted Z-scores for profiling (best Z, adjusted for age + sex)
z_cols <- grep("_sexadj_z$", names(rc), value = TRUE)
if (length(z_cols) == 0) z_cols <- grep("_local_z$", names(rc), value = TRUE)

# Pick highest-coverage Z-score as primary
z_coverage <- sapply(z_cols, function(zc) sum(!is.na(rc[[zc]])))
primary_z <- z_cols[which.max(z_coverage)]
message(">> Primary Z-score for profiling: ", primary_z,
        " (n = ", max(z_coverage), ")")

# Load raw clinical data for ALL biomarkers
raw <- cohort
raw$patient_id <- raw$Protocolo
raw$age <- as.numeric(raw$age)
raw <- raw[!is.na(raw$age) & raw$age >= 18 & raw$age <= 120, ]

# Merge recentered scores with raw biomarkers
prof_df <- merge(rc, raw[, c("patient_id", setdiff(names(raw),
            c(names(rc), "Protocolo", "source_file",
              "date_admission", "date_result")))],
            by = "patient_id", all.x = TRUE)

# Convert biomarkers to numeric
all_biomarkers <- c(
  "hemoglobin", "hematocrit", "platelets", "mch", "mchc",
  "neutrophils", "eosinophils", "basophils", "monocytes",
  "lymphocyte_abs", "band_neutrophils",
  "cholesterol", "hdl", "ldl", "non_hdl", "triglycerides", "lp_a",
  "glucose", "hba1c", "insulin",
  "ast", "alt", "ggt", "ldh", "albumin", "alp",
  "ferritin", "transferrin", "transferrin_saturation",
  "creatinine", "bun", "uric_acid", "eGFR",
  "crp", "fibrinogen", "homocysteine",
  "TSH", "free_t4",
  "vitamin_d", "vitamin_b12",
  "wbc", "rbc", "mcv", "rdw", "lymphocyte"
)
all_biomarkers <- intersect(all_biomarkers, names(prof_df))
for (bm in all_biomarkers) prof_df[[bm]] <- suppressWarnings(as.numeric(prof_df[[bm]]))

# ── 4a. Classify fast vs slow agers (quintile-based) ──

ok <- !is.na(prof_df[[primary_z]]) & is.finite(prof_df[[primary_z]])
df_valid <- prof_df[ok, ]
n_valid <- nrow(df_valid)
message(">> Valid subjects for profiling: ", n_valid)

quants <- quantile(df_valid[[primary_z]], probs = c(0.10, 0.25, 0.75, 0.90),
                   na.rm = TRUE)
df_valid$ager_group <- cut(df_valid[[primary_z]],
  breaks = c(-Inf, quants[1], quants[2], quants[3], quants[4], Inf),
  labels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")
)

cat("\n--- Ager Group Distribution ---\n\n")
tbl <- table(df_valid$ager_group)
for (g in names(tbl))
  cat(sprintf("   %-15s n = %6d (%5.1f%%)\n", g, tbl[g], 100 * tbl[g] / n_valid))

# ── 4b. Biomarker profiling across ager groups ──

cat("\n--- Biomarker Profiles: Fast vs Slow Agers ---\n\n")

usable_bm <- character(0)
for (bm in all_biomarkers) {
  n_ok <- sum(!is.na(df_valid[[bm]]) & is.finite(df_valid[[bm]]))
  if (n_ok >= 50) usable_bm <- c(usable_bm, bm)
}
message(">> Usable biomarkers: ", length(usable_bm))

profile_results <- list()

for (bm in usable_bm) {
  x <- df_valid[[bm]]
  g <- df_valid$ager_group
  ok_bm <- !is.na(x) & is.finite(x) & !is.na(g)
  if (sum(ok_bm) < 30) next

  kw <- tryCatch(kruskal.test(x[ok_bm] ~ g[ok_bm]), error = function(e) NULL)
  ct <- cor.test(df_valid[[primary_z]][ok_bm], x[ok_bm], method = "spearman")

  p1 <- x[ok_bm & g == "P1_very_slow"]
  p5 <- x[ok_bm & g == "P5_very_fast"]
  cohens_d <- NA; d_ci_lo <- NA; d_ci_hi <- NA
  n1 <- length(p1); n5 <- length(p5)
  if (n1 >= 5 && n5 >= 5) {
    pooled_sd <- sqrt((var(p1) * (n1 - 1) + var(p5) * (n5 - 1)) / (n1 + n5 - 2))
    if (pooled_sd > 0) {
      cohens_d <- (mean(p5) - mean(p1)) / pooled_sd
      # Analytical SE: Hedges & Olkin (1985)
      se_d <- sqrt((n1 + n5) / (n1 * n5) + cohens_d^2 / (2 * (n1 + n5)))
      d_ci_lo <- cohens_d - 1.96 * se_d
      d_ci_hi <- cohens_d + 1.96 * se_d
    }
  }

  profile_results[[bm]] <- data.frame(
    biomarker  = bm,
    n_total    = sum(ok_bm),
    n_P1       = n1,
    n_P5       = n5,
    rho_z      = round(ct$estimate, 4),
    rho_p      = ct$p.value,
    kw_chi2    = if (!is.null(kw)) round(kw$statistic, 2) else NA,
    kw_p       = if (!is.null(kw)) kw$p.value else NA,
    d_P5_vs_P1 = round(cohens_d, 3),
    d_ci_lo    = round(d_ci_lo, 3),
    d_ci_hi    = round(d_ci_hi, 3),
    mean_P1    = if (n1 >= 5) round(mean(p1), 3) else NA,
    mean_P5    = if (n5 >= 5) round(mean(p5), 3) else NA,
    direction  = ifelse(!is.na(cohens_d),
                        ifelse(cohens_d > 0.1, "HIGHER_in_fast",
                        ifelse(cohens_d < -0.1, "LOWER_in_fast", "~equal")),
                        NA),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

profile_table <- NULL
if (length(profile_results) > 0) {
  profile_table <- do.call(rbind, profile_results)
  profile_table$rho_padj <- p.adjust(profile_table$rho_p, method = "BH")
  profile_table$kw_padj  <- p.adjust(profile_table$kw_p, method = "BH")
  profile_table <- profile_table[order(-abs(profile_table$d_P5_vs_P1)), ]
  rownames(profile_table) <- NULL

  sig <- profile_table[!is.na(profile_table$d_P5_vs_P1) &
                       abs(profile_table$d_P5_vs_P1) > 0.1, ]
  if (nrow(sig) > 0) {
    cat("  Significantly differentiating biomarkers (|d| > 0.1):\n\n")
    print(sig[, c("biomarker", "n_total", "n_P1", "n_P5", "rho_z", "d_P5_vs_P1",
                   "d_ci_lo", "d_ci_hi", "mean_P1", "mean_P5", "direction", "kw_padj")],
          row.names = FALSE)
  }

  write.csv(profile_table, file.path(OUTPUT_DIR, "step4_ager_biomarker_profiles.csv"),
            row.names = FALSE)
  message(">> Saved: step4_ager_biomarker_profiles.csv")
}

# ── 4c. Aging signature heatmap ──

if (!is.null(profile_table) && nrow(profile_table) > 0) {
  sig_bm <- profile_table$biomarker[
    !is.na(profile_table$d_P5_vs_P1) &
    abs(profile_table$d_P5_vs_P1) > 0.15 &
    !is.na(profile_table$kw_padj) &
    profile_table$kw_padj < 0.05
  ]

  if (length(sig_bm) >= 3) {
    cat(sprintf("\n>> Building aging signature from %d significant biomarkers\n\n",
                length(sig_bm)))

    radar_data <- list()
    for (bm in sig_bm) {
      if (!bm %in% names(df_valid)) next
      x <- df_valid[[bm]]
      ok_bm <- !is.na(x) & is.finite(x) & !is.na(df_valid$ager_group)
      if (sum(ok_bm) < 30) next

      x_z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      group_means <- tapply(x_z[ok_bm], df_valid$ager_group[ok_bm], mean)

      for (g in names(group_means)) {
        radar_data[[length(radar_data) + 1]] <- data.frame(
          biomarker = bm, group = g, mean_z = round(group_means[g], 3),
          stringsAsFactors = FALSE, row.names = NULL
        )
      }
    }

    if (length(radar_data) > 0) {
      radar_table <- do.call(rbind, radar_data)
      radar_wide <- pivot_wider(radar_table, names_from = group, values_from = mean_z)
      write.csv(radar_wide, file.path(OUTPUT_DIR, "step4_aging_signature_data.csv"),
                row.names = FALSE)

      radar_table$group <- factor(radar_table$group,
        levels = c("P1_very_slow", "P2_slow", "P3_average",
                   "P4_fast", "P5_very_fast"))

      bm_order <- radar_table %>%
        filter(group %in% c("P1_very_slow", "P5_very_fast")) %>%
        pivot_wider(names_from = group, values_from = mean_z) %>%
        mutate(diff = P5_very_fast - P1_very_slow) %>%
        arrange(diff)
      radar_table$biomarker <- factor(radar_table$biomarker,
                                       levels = bm_order$biomarker)

      p_heat <- ggplot(radar_table, aes(x = group, y = biomarker, fill = mean_z)) +
        geom_tile(colour = "white", linewidth = 0.5) +
        geom_text(aes(label = sprintf("%.2f", mean_z)), size = 3) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                             midpoint = 0, name = "Mean Z") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Accelerated Aging Biomarker Signature",
             subtitle = paste0("Based on ", primary_z, " quintiles"),
             x = "Ager Group", y = "Biomarker")

      ggsave(file.path(OUTPUT_DIR, "step4_aging_signature_heatmap.png"),
             p_heat, width = 10, height = max(6, length(sig_bm) * 0.35 + 2),
             dpi = 150)
      print(p_heat)
      message(">> Saved: step4_aging_signature_heatmap.png")

      cat("\n--- Accelerated Aging Signature ---\n\n")
      print(as.data.frame(radar_wide), row.names = FALSE)
    }
  }
}

# ── 4d. Per-clock fast/slow ager profiling ──

cat("\n--- Per-Clock Ager Classification ---\n\n")
cat("  Comparing ager groups defined by EACH sub-clock independently\n\n")

for (zc in z_cols) {
  clock_label <- sub("_sexadj_z$|_local_z$", "", zc)
  z_vals <- prof_df[[zc]]
  ok_z <- !is.na(z_vals) & is.finite(z_vals)
  if (sum(ok_z) < 100) next

  # Classify using this clock's Z-score
  qts <- quantile(z_vals[ok_z], probs = c(0.25, 0.75))
  slow <- ok_z & z_vals <= qts[1]
  fast <- ok_z & z_vals >= qts[2]

  cat(sprintf("  %s: slow (Q1) = %d, fast (Q4) = %d\n",
              clock_label, sum(slow), sum(fast)))

  # For each usable biomarker, compare fast vs slow
  clock_diffs <- list()
  for (bm in usable_bm) {
    x <- prof_df[[bm]]
    x_slow <- x[slow & !is.na(x) & is.finite(x)]
    x_fast <- x[fast & !is.na(x) & is.finite(x)]
    if (length(x_slow) < 20 || length(x_fast) < 20) next

    ps <- sqrt((var(x_slow) * (length(x_slow) - 1) + var(x_fast) * (length(x_fast) - 1)) /
               (length(x_slow) + length(x_fast) - 2))
    d <- if (ps > 0) (mean(x_fast) - mean(x_slow)) / ps else NA

    clock_diffs[[bm]] <- data.frame(
      clock = clock_label, biomarker = bm,
      d = round(d, 3),
      mean_slow = round(mean(x_slow), 3),
      mean_fast = round(mean(x_fast), 3),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }

  if (length(clock_diffs) > 0) {
    cd <- do.call(rbind, clock_diffs)
    cd <- cd[order(-abs(cd$d)), ]
    # Show top 5
    top <- head(cd, 5)
    for (i in seq_len(nrow(top))) {
      cat(sprintf("    %-25s d = %+.3f  (slow: %.1f → fast: %.1f)\n",
                  top$biomarker[i], top$d[i], top$mean_slow[i], top$mean_fast[i]))
    }
    cat("\n")
  }
}

# ── 4d-bis. FULL biomarker profile using pheno_max_v2 as classifier ──

cat("\n##########################################################\n")
cat("#  FULL BIOMARKER PROFILE — pheno_max_v2 classifier      #\n")
cat("##########################################################\n\n")

maxv2_z <- "pheno_pheno_max_v2_sexadj_z"
if (!maxv2_z %in% names(prof_df)) maxv2_z <- "pheno_pheno_max_v2_local_z"

if (maxv2_z %in% names(prof_df)) {
  z_v2 <- prof_df[[maxv2_z]]
  ok_v2 <- !is.na(z_v2) & is.finite(z_v2)
  n_v2 <- sum(ok_v2)
  message(">> pheno_max_v2 classifier: ", n_v2, " subjects")

  # Quintile classification
  qts_v2 <- quantile(z_v2[ok_v2], probs = c(0.10, 0.25, 0.75, 0.90), na.rm = TRUE)
  prof_df$ager_maxv2 <- NA_character_
  prof_df$ager_maxv2[ok_v2] <- as.character(cut(z_v2[ok_v2],
    breaks = c(-Inf, qts_v2[1], qts_v2[2], qts_v2[3], qts_v2[4], Inf),
    labels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")
  ))

  cat("--- Ager Groups (pheno_max_v2) ---\n\n")
  tbl_v2 <- table(prof_df$ager_maxv2[ok_v2])
  for (g in names(tbl_v2))
    cat(sprintf("   %-15s n = %6d (%5.1f%%)\n", g, tbl_v2[g], 100 * tbl_v2[g] / n_v2))

  # Full profiling
  v2_valid <- prof_df[ok_v2, ]
  v2_results <- list()

  for (bm in usable_bm) {
    x <- v2_valid[[bm]]
    g <- v2_valid$ager_maxv2
    ok_bm <- !is.na(x) & is.finite(x) & !is.na(g)
    if (sum(ok_bm) < 30) next

    kw <- tryCatch(kruskal.test(x[ok_bm] ~ factor(g[ok_bm])), error = function(e) NULL)
    ct <- cor.test(z_v2[ok_v2][ok_bm], x[ok_bm], method = "spearman")

    p1 <- x[ok_bm & g == "P1_very_slow"]
    p5 <- x[ok_bm & g == "P5_very_fast"]
    cohens_d <- NA; d_ci_lo <- NA; d_ci_hi <- NA
    n1 <- length(p1); n5 <- length(p5)
    if (n1 >= 5 && n5 >= 5) {
      pooled_sd <- sqrt((var(p1) * (n1 - 1) + var(p5) * (n5 - 1)) / (n1 + n5 - 2))
      if (pooled_sd > 0) {
        cohens_d <- (mean(p5) - mean(p1)) / pooled_sd
        se_d <- sqrt((n1 + n5) / (n1 * n5) + cohens_d^2 / (2 * (n1 + n5)))
        d_ci_lo <- cohens_d - 1.96 * se_d
        d_ci_hi <- cohens_d + 1.96 * se_d
      }
    }

    v2_results[[bm]] <- data.frame(
      biomarker  = bm,
      n_total    = sum(ok_bm),
      n_P1       = n1,
      n_P5       = n5,
      rho_z      = round(ct$estimate, 4),
      rho_p      = ct$p.value,
      kw_chi2    = if (!is.null(kw)) round(kw$statistic, 2) else NA,
      kw_p       = if (!is.null(kw)) kw$p.value else NA,
      d_P5_vs_P1 = round(cohens_d, 3),
      d_ci_lo    = round(d_ci_lo, 3),
      d_ci_hi    = round(d_ci_hi, 3),
      mean_P1    = if (n1 >= 5) round(mean(p1), 3) else NA,
      mean_P5    = if (n5 >= 5) round(mean(p5), 3) else NA,
      direction  = ifelse(!is.na(cohens_d),
                          ifelse(cohens_d > 0.1, "HIGHER_in_fast",
                          ifelse(cohens_d < -0.1, "LOWER_in_fast", "~equal")),
                          NA),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }

  if (length(v2_results) > 0) {
    v2_table <- do.call(rbind, v2_results)
    v2_table$rho_padj <- p.adjust(v2_table$rho_p, method = "BH")
    v2_table$kw_padj  <- p.adjust(v2_table$kw_p, method = "BH")
    v2_table <- v2_table[order(-abs(v2_table$d_P5_vs_P1)), ]
    rownames(v2_table) <- NULL

    sig_v2 <- v2_table[!is.na(v2_table$d_P5_vs_P1) &
                        abs(v2_table$d_P5_vs_P1) > 0.1, ]
    if (nrow(sig_v2) > 0) {
      cat("\n  Biomarkers differentiating fast vs slow agers (pheno_max_v2, |d| > 0.1):\n\n")
      print(sig_v2[, c("biomarker", "n_total", "n_P1", "n_P5", "rho_z", "d_P5_vs_P1",
                        "d_ci_lo", "d_ci_hi", "mean_P1", "mean_P5", "direction", "kw_padj")],
            row.names = FALSE)
    }

    write.csv(v2_table, file.path(OUTPUT_DIR, "step4_maxv2_biomarker_profiles.csv"),
              row.names = FALSE)
    message(">> Saved: step4_maxv2_biomarker_profiles.csv")

    # Heatmap for pheno_max_v2
    sig_bm_v2 <- v2_table$biomarker[
      !is.na(v2_table$d_P5_vs_P1) &
      abs(v2_table$d_P5_vs_P1) > 0.15 &
      !is.na(v2_table$kw_padj) &
      v2_table$kw_padj < 0.05
    ]

    if (length(sig_bm_v2) >= 3) {
      radar_v2 <- list()
      for (bm in sig_bm_v2) {
        if (!bm %in% names(v2_valid)) next
        x <- v2_valid[[bm]]
        ok_bm <- !is.na(x) & is.finite(x) & !is.na(v2_valid$ager_maxv2)
        if (sum(ok_bm) < 30) next

        x_z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
        group_means <- tapply(x_z[ok_bm], v2_valid$ager_maxv2[ok_bm], mean)
        for (g in names(group_means)) {
          radar_v2[[length(radar_v2) + 1]] <- data.frame(
            biomarker = bm, group = g, mean_z = round(group_means[g], 3),
            stringsAsFactors = FALSE, row.names = NULL
          )
        }
      }

      if (length(radar_v2) > 0) {
        radar_v2_df <- do.call(rbind, radar_v2)
        radar_v2_wide <- pivot_wider(radar_v2_df, names_from = group, values_from = mean_z)
        write.csv(radar_v2_wide, file.path(OUTPUT_DIR, "step4_maxv2_aging_signature.csv"),
                  row.names = FALSE)

        radar_v2_df$group <- factor(radar_v2_df$group,
          levels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast"))

        bm_ord_v2 <- radar_v2_df %>%
          filter(group %in% c("P1_very_slow", "P5_very_fast")) %>%
          pivot_wider(names_from = group, values_from = mean_z) %>%
          mutate(diff = P5_very_fast - P1_very_slow) %>%
          arrange(diff)
        radar_v2_df$biomarker <- factor(radar_v2_df$biomarker, levels = bm_ord_v2$biomarker)

        p_v2 <- ggplot(radar_v2_df, aes(x = group, y = biomarker, fill = mean_z)) +
          geom_tile(colour = "white", linewidth = 0.5) +
          geom_text(aes(label = sprintf("%.2f", mean_z)), size = 3) +
          scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                               midpoint = 0, name = "Mean Z") +
          theme_minimal(base_size = 12) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "Accelerated Aging Signature (pheno_max_v2 classifier)",
               subtitle = sprintf("6/9 Levine biomarkers — n = %s",
                                   format(n_v2, big.mark = ",")),
               x = "Ager Group", y = "Biomarker")

        ggsave(file.path(OUTPUT_DIR, "step4_maxv2_aging_signature_heatmap.png"),
               p_v2, width = 10, height = max(6, length(sig_bm_v2) * 0.35 + 2),
               dpi = 150)
        print(p_v2)
        message(">> Saved: step4_maxv2_aging_signature_heatmap.png")

        cat("\n--- Accelerated Aging Signature (pheno_max_v2) ---\n\n")
        print(as.data.frame(radar_v2_wide), row.names = FALSE)
      }
    }
  }
} else {
  message(">> pheno_max_v2 Z-score not found — skipping")
}

# ── 4e. Biomarker profiles by age decade and sex ──

cat("\n--- Biomarker Trends by Age Decade & Sex ---\n\n")

prof_df$age_decade <- cut(prof_df$age,
  breaks = c(0, 30, 40, 50, 60, 70, 80, 120),
  labels = c("<30", "30s", "40s", "50s", "60s", "70s", "80+"),
  right = FALSE
)

# Select key biomarkers that differentiate fast/slow agers
key_bm <- character(0)
if (!is.null(profile_table) && nrow(profile_table) > 0) {
  key_bm <- profile_table$biomarker[
    !is.na(profile_table$d_P5_vs_P1) &
    abs(profile_table$d_P5_vs_P1) > 0.15 &
    !is.na(profile_table$kw_padj) &
    profile_table$kw_padj < 0.05
  ]
  key_bm <- head(key_bm, 12)  # Top 12 for readability
}

if (length(key_bm) >= 3) {
  for (bm in key_bm) {
    if (!bm %in% names(prof_df)) next
    x <- suppressWarnings(as.numeric(prof_df[[bm]]))
    ok <- !is.na(x) & is.finite(x) & !is.na(prof_df$age_decade) &
          prof_df$sex %in% c("F", "M")
    if (sum(ok) < 100) next

    plot_df <- data.frame(
      value = x[ok],
      age_decade = prof_df$age_decade[ok],
      sex = prof_df$sex[ok]
    )

    p_bm <- ggplot(plot_df, aes(x = age_decade, y = value, fill = sex)) +
      geom_boxplot(alpha = 0.6, outlier.size = 0.3, outlier.alpha = 0.1) +
      scale_fill_manual(values = c("F" = "firebrick", "M" = "steelblue")) +
      theme_minimal(base_size = 12) +
      labs(title = paste0("Biomarker by Age & Sex: ", bm),
           x = "Age Decade", y = bm, fill = "Sex")
    print(p_bm)

    ggsave(file.path(OUTPUT_DIR, paste0("step4_biomarker_age_sex_", bm, ".png")),
           p_bm, width = 9, height = 5, dpi = 150)
  }
  message(">> Saved biomarker × age × sex plots for ", length(key_bm), " biomarkers")
}

# ── 4f. Advancement by age decade (recentered) for all clocks ──

cat("\n--- Advancement Distribution by Age Decade (Recentered) ---\n\n")

for (zc in sexadj_z_cols) {
  clock_label <- sub("_sexadj_z$", "", zc)
  short_label <- gsub("pheno_", "", clock_label)

  ok <- !is.na(rc[[zc]]) & is.finite(rc[[zc]]) & !is.na(rc$age_stratum)
  if (sum(ok) < 100) next

  p_box <- ggplot(rc[ok, ], aes(x = age_stratum, y = .data[[zc]])) +
    geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Recentered Z-Score by Age: ", short_label),
         subtitle = "After sex-adjusted residualization",
         x = "Age Group", y = "Z-score (0 = population average)")
  print(p_box)

  ggsave(file.path(OUTPUT_DIR, paste0("step4_zscore_by_age_", short_label, ".png")),
         p_box, width = 8, height = 5, dpi = 150)
}


# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  PIPELINE COMPLETE                                      ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

cat("  STEP 1 — NHANES Training\n")
if (exists("nhanes_table") && !is.null(nhanes_table)) {
  for (i in seq_len(nrow(nhanes_table))) {
    ct <- nhanes_table[i, ]
    cat(sprintf("    %-35s r = %.4f  RMSE = %.2f  (n = %s)\n",
                ct$clock, ct$r, ct$RMSE, format(ct$n, big.mark = ",")))
  }
}

cat("\n  STEP 2 — Argentine Scoring\n")
cat(sprintf("    Cohort: %s patients (%.0f%% F, age %.0f ± %.0f)\n",
            format(nrow(results), big.mark = ","),
            100 * mean(results$sex == "F", na.rm = TRUE),
            mean(results$age, na.rm = TRUE),
            sd(results$age, na.rm = TRUE)))
if (exists("arg_table") && !is.null(arg_table)) {
  for (i in seq_len(nrow(arg_table))) {
    ct <- arg_table[i, ]
    cat(sprintf("    %-35s r = %.4f  RMSE = %.2f  n = %s\n",
                ct$clock, ct$r, ct$RMSE, format(ct$n, big.mark = ",")))
  }
}

cat("\n  STEP 3 — Recentering\n")
cat(sprintf("    Clocks recentered: %d\n", length(viable_adv)))
cat("    Methods: A (age-residual), B (age+sex), C (age-stratified)\n")

cat("\n  STEP 4 — Fast/Slow Ager Profiling\n")
cat(sprintf("    Subjects profiled: %s\n", format(n_valid, big.mark = ",")))
cat(sprintf("    Primary Z-score: %s\n", primary_z))
if (!is.null(profile_table)) {
  n_sig <- sum(!is.na(profile_table$d_P5_vs_P1) &
               abs(profile_table$d_P5_vs_P1) > 0.15 &
               !is.na(profile_table$kw_padj) &
               profile_table$kw_padj < 0.05)
  cat(sprintf("    Significant biomarkers (|d| > 0.15): %d / %d\n",
              n_sig, nrow(profile_table)))
}

cat("\n  Output directory: ", OUTPUT_DIR, "\n")
cat("  Files saved:\n")
saved <- list.files(OUTPUT_DIR, full.names = FALSE)
for (f in saved) cat("    - ", f, "\n")

cat("\n  Key data files:\n")
cat("    - bioage_deployment_bundle.rds    (trained models)\n")
cat("    - clinical_scored_lite.csv        (raw scores)\n")
cat("    - final_output/clinical_recentered.csv (recentered)\n")
cat("    - final_output/step4_aging_signature_heatmap.png\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  Done. Results ready for manuscript.                    ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")

# Store everything
final_results <- list(
  nhanes_perf     = if (exists("nhanes_table")) nhanes_table else NULL,
  argentina_perf  = if (exists("arg_table")) arg_table else NULL,
  scored_data     = results,
  recentered_data = rc,
  profile_table   = profile_table,
  bundle          = bundle,
  primary_z       = primary_z
)
message(">> All results in 'final_results'")
