###############################################################################
# run_clinical_lite.R — Score clinical cohort using BioAge::phenoage_calc()
#
# Uses the BioAge package function directly instead of reimplementing the
# Gompertz formula.  This guarantees the scoring matches the training exactly.
#
# Requirements:
#   - BioAge package installed
#   - bioage_deployment_bundle.rds (from run_pipeline.R)
#   - phenoage_master_PARTE1.csv (+ optional PARTE2.csv)
#
# Usage:  source("run_clinical_lite.R")
###############################################################################

cat("==========================================================\n")
cat("  BIOLOGICAL AGE — CLINICAL SCORING (LITE)\n")
cat("  Using BioAge::phenoage_calc() directly\n")
cat("==========================================================\n\n")

# ── 0. Setup ──────────────────────────────────────────────────────────────

PROJECT_ROOT <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) getwd())
PROJECT_ROOT <- normalizePath(PROJECT_ROOT)
setwd(PROJECT_ROOT)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(BioAge)
  library(flexsurv)
})

source(file.path(PROJECT_ROOT, "R", "config.R"), local = FALSE)
source(file.path(PROJECT_ROOT, "R", "export.R"), local = FALSE)

# ── 1. Load deployment bundle ─────────────────────────────────────────────

BUNDLE_PATH <- file.path(PROJECT_ROOT, "bioage_deployment_bundle.rds")
if (!file.exists(BUNDLE_PATH)) {
  stop("Bundle not found. Run source('run_pipeline.R') first.")
}

message(">> Loading bundle...")
bundle <- load_bundle(BUNDLE_PATH)
message(">> Clocks: ", paste(bundle$meta$clocks_trained, collapse = ", "))

# ── 2. Load clinical data ────────────────────────────────────────────────

message("\n>> Loading clinical data...")
f1 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE2.csv")

parte1 <- read.csv(f1, stringsAsFactors = FALSE)
cohort <- parte1
if (file.exists(f2)) {
  parte2 <- read.csv(f2, stringsAsFactors = FALSE)
  cohort <- rbind(parte1, parte2)
}
message(">> Rows: ", nrow(cohort))

# ── 3. Transform to NHANES column names ──────────────────────────────────

message(">> Transforming biomarkers...")

df <- cohort %>%
  mutate(
    patient_id = Protocolo,
    age = as.numeric(age),

    # Levine 9 biomarkers (NHANES naming)
    albumin = as.numeric(albumin),
    alp     = as.numeric(alp),
    glucose = as.numeric(glucose),
    crp     = as.numeric(crp),
    lncrp   = ifelse(is.na(crp) | crp <= 0, NA, log(crp)),
    creat   = as.numeric(creatinine),
    lncreat = ifelse(is.na(creat) | creat <= 0, NA, log(creat)),
    lymph   = as.numeric(lymphocyte),
    mcv     = as.numeric(mcv),
    rdw     = as.numeric(rdw),
    wbc     = ifelse(as.numeric(wbc) > 300, as.numeric(wbc) / 1000, as.numeric(wbc)),

    # Sub-clock biomarkers
    rbc        = as.numeric(rbc),
    rbc        = ifelse(rbc > 100, rbc / 1e6, rbc),
    ggt        = as.numeric(ggt),
    insulin    = as.numeric(insulin),
    trig       = as.numeric(triglycerides),
    totchol    = as.numeric(cholesterol),
    hba1c      = as.numeric(hba1c),
    vitaminB12 = as.numeric(vitamin_b12),
    bun        = as.numeric(bun),
    uap        = as.numeric(uric_acid),
    gender     = ifelse(sex == "M", 1, 2)
  )

# Filter invalid ages
df <- df[!is.na(df$age) & df$age >= 1 & df$age <= 120, ]
message(">> After age filter: ", nrow(df), " rows")
message(">> Age: ", round(mean(df$age, na.rm=TRUE), 1), " +/- ",
        round(sd(df$age, na.rm=TRUE), 1), " years")

# ══════════════════════════════════════════════════════════════════════════
# 4. BIOMARKER AVAILABILITY
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("BIOMARKER AVAILABILITY")
message(strrep("=", 60))

all_bm <- c("albumin", "alp", "glucose", "lncrp", "lncreat",
            "lymph", "mcv", "rdw", "wbc",
            "rbc", "ggt", "insulin", "trig", "totchol",
            "vitaminB12", "hba1c", "bun")
for (v in all_bm) {
  if (v %in% names(df)) {
    n_ok <- sum(!is.na(df[[v]]))
    cat(sprintf("  %-14s %6d / %d  (%5.1f%%)\n",
                v, n_ok, nrow(df), 100 * n_ok / nrow(df)))
  }
}

# ══════════════════════════════════════════════════════════════════════════
# 5. SCORE CLOCKS USING BioAge::phenoage_calc() DIRECTLY
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("SCORING CLOCKS VIA BioAge::phenoage_calc()")
message(strrep("=", 60))

# Define clocks to score: name → list(biomarkers, fit)
clocks_to_score <- list()

# Original PhenoAge (NHANES-trained)
if (!is.null(bundle$models$pheno_orig)) {
  clocks_to_score$phenoage_orig <- list(
    biomarkers = BIOMARKERS_PHENO_LEVINE,
    fit        = bundle$models$pheno_orig
  )
}

# Sub-clocks
for (nm in names(bundle$subclock_models)) {
  panel_name <- nm
  if (panel_name %in% names(bundle$panels$subclocks)) {
    clocks_to_score[[paste0("pheno_", nm)]] <- list(
      biomarkers = bundle$panels$subclocks[[panel_name]],
      fit        = bundle$subclock_models[[nm]]
    )
  }
}

# Master results data frame
results <- data.frame(
  patient_id = df$patient_id,
  age        = df$age,
  sex        = df$sex,
  stringsAsFactors = FALSE
)

# Score each clock
score_summary <- list()

for (clock_name in names(clocks_to_score)) {
  clock <- clocks_to_score[[clock_name]]
  bm    <- clock$biomarkers
  fit   <- clock$fit

  # Check which biomarker columns exist and have data
  bm_plus_age <- c(bm, "age")
  missing_cols <- setdiff(bm_plus_age, names(df))

  if (length(missing_cols) > 0) {
    message("\n>> ", clock_name, ": SKIP — missing columns: ",
            paste(missing_cols, collapse = ", "))
    next
  }

  # Find complete cases (patients with ALL biomarkers for this clock)
  complete_mask <- complete.cases(df[, bm_plus_age, drop = FALSE])
  n_complete <- sum(complete_mask)

  if (n_complete == 0) {
    message("\n>> ", clock_name, ": SKIP — 0 complete cases")
    next
  }

  message("\n>> ", clock_name, ": ", n_complete, " complete cases (",
          round(100 * n_complete / nrow(df), 1), "%)")
  message("   Biomarkers: ", paste(bm, collapse = ", "))

  # Subset to complete cases only
  dat_subset <- df[complete_mask, c(bm_plus_age), drop = FALSE]

  # Call BioAge::phenoage_calc() directly — it handles the entire formula
  calc_result <- tryCatch({
    phenoage_calc(
      data       = dat_subset,
      biomarkers = bm,
      fit        = fit
    )
  }, error = function(e) {
    message("   ERROR: ", e$message)
    NULL
  })

  if (is.null(calc_result)) next

  # Extract results
  pa_vals  <- calc_result$data$phenoage
  pa_adv   <- calc_result$data$phenoage_advance
  n_valid  <- sum(!is.na(pa_vals) & is.finite(pa_vals))

  message("   -> ", n_valid, " / ", n_complete, " valid PhenoAge values")

  if (n_valid > 0) {
    # Map back to full dataset
    ba_col  <- clock_name
    adv_col <- paste0(clock_name, "_advance")
    results[[ba_col]]  <- NA_real_
    results[[adv_col]] <- NA_real_
    results[[ba_col]][complete_mask]  <- pa_vals
    results[[adv_col]][complete_mask] <- pa_adv

    score_summary[[clock_name]] <- data.frame(
      clock   = clock_name,
      n_scored = n_valid,
      pct     = round(100 * n_valid / nrow(df), 1),
      stringsAsFactors = FALSE
    )
  }
}

# ══════════════════════════════════════════════════════════════════════════
# 6. STATISTICAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("RESULTS SUMMARY")
message(strrep("=", 60))

# Cohort demographics
cat("\n--- Cohort Demographics ---\n")
cat(sprintf("  Total:    %d subjects\n", nrow(results)))
cat(sprintf("  Female:   %d (%.1f%%)\n",
            sum(results$sex == "F", na.rm=TRUE),
            100 * mean(results$sex == "F", na.rm=TRUE)))
cat(sprintf("  Male:     %d (%.1f%%)\n",
            sum(results$sex == "M", na.rm=TRUE),
            100 * mean(results$sex == "M", na.rm=TRUE)))
cat(sprintf("  Age:      %.1f +/- %.1f  [%d - %d]\n",
            mean(results$age, na.rm=TRUE), sd(results$age, na.rm=TRUE),
            min(results$age, na.rm=TRUE), max(results$age, na.rm=TRUE)))

# Scoring summary
if (length(score_summary) > 0) {
  cat("\n--- Clocks Scored ---\n")
  ss <- do.call(rbind, score_summary)
  print(ss, row.names = FALSE)
} else {
  cat("\n>> WARNING: No clocks scored any patients.\n")
}

# Correlation with age
cat("\n--- Correlation of Biological Age with Chronological Age ---\n")
cat("  (A good clock should have r > 0.8)\n\n")

ba_cols <- names(score_summary)
cor_stats <- list()

for (bc in ba_cols) {
  vals <- results[[bc]]
  ages <- results$age
  ok   <- !is.na(vals) & !is.na(ages) & is.finite(vals)
  n_ok <- sum(ok)

  if (n_ok >= 10) {
    ct     <- cor.test(ages[ok], vals[ok])
    lm_fit <- lm(vals[ok] ~ ages[ok])
    s      <- summary(lm_fit)

    cor_stats[[bc]] <- data.frame(
      clock     = bc,
      n         = n_ok,
      r         = round(ct$estimate, 4),
      r_sq      = round(s$r.squared, 4),
      p         = ifelse(ct$p.value < 2.2e-16, "< 2.2e-16",
                         sprintf("%.2e", ct$p.value)),
      slope     = round(coef(lm_fit)[2], 4),
      intercept = round(coef(lm_fit)[1], 2),
      RMSE      = round(sqrt(mean(s$residuals^2)), 2),
      stringsAsFactors = FALSE
    )

    cat(sprintf("  %-40s r = %.4f  (n = %d)\n", bc, ct$estimate, n_ok))
  }
}

cor_table <- NULL
if (length(cor_stats) > 0) {
  cor_table <- do.call(rbind, cor_stats)
  rownames(cor_table) <- NULL
  cat("\n")
  print(cor_table, row.names = FALSE)
}

# Advancement statistics
adv_cols <- grep("_advance$", names(results), value = TRUE)
if (length(adv_cols) > 0) {
  cat("\n--- Advancement (BA - CA) ---\n")
  for (ac in adv_cols) {
    vals <- results[[ac]]
    ok <- !is.na(vals) & is.finite(vals)
    if (sum(ok) >= 10) {
      cat(sprintf("  %-40s mean = %+.2f  SD = %.2f  median = %+.2f  (n = %d)\n",
                  sub("_advance$", "", ac),
                  mean(vals[ok]), sd(vals[ok]), median(vals[ok]), sum(ok)))
    }
  }
}

# ══════════════════════════════════════════════════════════════════════════
# 7. PLOTS
# ══════════════════════════════════════════════════════════════════════════

message("\n>> Generating plots...")

# Plot: BA vs CA for each clock
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

# Plot: Advancement distribution
valid_adv <- character(0)
for (ac in adv_cols) {
  if (sum(!is.na(results[[ac]]) & is.finite(results[[ac]])) > 30)
    valid_adv <- c(valid_adv, ac)
}

if (length(valid_adv) > 0) {
  long_adv <- pivot_longer(
    results[, c("age", valid_adv), drop = FALSE],
    cols = all_of(valid_adv),
    names_to = "clock", values_to = "advancement"
  )
  long_adv$clock <- gsub("pheno_|phenoage_|_advance", "", long_adv$clock)

  p_adv <- ggplot(long_adv[!is.na(long_adv$advancement) &
                            is.finite(long_adv$advancement), ],
                  aes(x = advancement)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                   fill = "steelblue", alpha = 0.6, colour = "white") +
    geom_density(colour = "navy", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
    facet_wrap(~ clock, scales = "free") +
    theme_minimal(base_size = 12) +
    labs(title = "Distribution of Biological Age Advancement (BA - CA)",
         x = "Advancement (years)", y = "Density")
  print(p_adv)
}

# Plot: Age distribution
p_age <- ggplot(results, aes(x = age)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7, colour = "white") +
  theme_minimal(base_size = 13) +
  labs(title = "Age Distribution",
       subtitle = sprintf("N = %s, Mean = %.1f",
                          format(nrow(results), big.mark=","),
                          mean(results$age, na.rm=TRUE)),
       x = "Age (years)", y = "Count")
print(p_age)

# Plot: BA vs CA by sex
for (bc in ba_cols) {
  vals <- results[[bc]]
  ok <- !is.na(vals) & is.finite(vals)
  if (sum(ok) < 50) next
  label <- gsub("pheno_|phenoage_", "", bc)

  p_sex <- ggplot(results[ok, ], aes(x = age, y = .data[[bc]], colour = sex)) +
    geom_point(alpha = 0.05, size = 0.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray40") +
    scale_colour_manual(values = c("F" = "firebrick", "M" = "steelblue")) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("BA vs CA by Sex: ", label),
         x = "Chronological Age", y = "Biological Age", colour = "Sex")
  print(p_sex)
}

# Plot: Advancement by age decade
if (length(valid_adv) > 0) {
  results$age_decade <- cut(results$age,
    breaks = c(0, 30, 40, 50, 60, 70, 80, 120),
    labels = c("<30", "30s", "40s", "50s", "60s", "70s", "80+"),
    right = FALSE)

  for (ac in valid_adv) {
    ok <- !is.na(results[[ac]]) & is.finite(results[[ac]])
    if (sum(ok) < 50) next
    label <- gsub("pheno_|phenoage_|_advance", "", ac)

    p_box <- ggplot(results[ok, ], aes(x = age_decade, y = .data[[ac]])) +
      geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.size = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      theme_minimal(base_size = 13) +
      labs(title = paste0("Advancement by Decade: ", label),
           x = "Age Decade", y = "Advancement (years)")
    print(p_box)
  }
}

# ══════════════════════════════════════════════════════════════════════════
# 8. SAVE
# ══════════════════════════════════════════════════════════════════════════

message("\n>> Saving outputs...")

write.csv(results, file.path(PROJECT_ROOT, "clinical_scored_lite.csv"),
          row.names = FALSE)
message(">> Saved: clinical_scored_lite.csv")

if (!is.null(cor_table)) {
  write.csv(cor_table, file.path(PROJECT_ROOT, "clinical_stats_lite.csv"),
            row.names = FALSE)
  message(">> Saved: clinical_stats_lite.csv")
}

# ══════════════════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════════════════

n_scored <- length(score_summary)
cat("\n")
cat("==========================================================\n")
cat("  SCORING COMPLETE\n")
cat("==========================================================\n")
cat("  Cohort:        ", format(nrow(results), big.mark=","), " subjects\n")
cat("  Clocks scored: ", n_scored, "\n")
if (!is.null(cor_table) && nrow(cor_table) > 0) {
  best <- cor_table[which.max(abs(cor_table$r)), ]
  cat("  Best r(BA,CA): ", best$r, " (", best$clock, ")\n")
}
cat("\n  Output: clinical_scored_lite.csv\n")
cat("          clinical_stats_lite.csv\n")
cat("==========================================================\n")

# Store results in environment
clinical_results_lite <- list(
  data      = results,
  cor_table = cor_table,
  bundle    = bundle
)
message(">> Results in 'clinical_results_lite'")
