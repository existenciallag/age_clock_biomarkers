###############################################################################
# run_clinical.R — Apply biological age clocks to preclinical cohort data
#
# Reads master_PARTE1.csv + master_PARTE2.csv, transforms biomarkers to
# NHANES conventions, computes PhenoAge (priority) + all available clocks,
# and produces statistical summaries with correlation vs age.
#
# Usage:
#   source("run_pipeline.R")   # Train models on NHANES first (once)
#   source("run_clinical.R")   # Then score your cohort
#
# Output:
#   - clinical_scored.csv       : wide-format scored dataset
#   - clinical_report.csv       : per-patient per-clock report
#   - clinical_stats_summary.csv: statistical summary (R vs age, etc.)
#   - Plots displayed in viewer
###############################################################################

cat("==========================================================\n")
cat("  BIOLOGICAL AGE — PRECLINICAL COHORT SCORING\n")
cat("  Applying NHANES-trained models to clinical lab data\n")
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
})

# Source modules
for (f in c("config.R", "export.R", "clinical_score.R", "residualize.R")) {
  source(file.path(PROJECT_ROOT, "R", f), local = FALSE)
}

# ── 1. Load deployment bundle ─────────────────────────────────────────────

BUNDLE_PATH <- file.path(PROJECT_ROOT, "bioage_deployment_bundle.rds")

if (!file.exists(BUNDLE_PATH)) {
  stop("Deployment bundle not found. Run source('run_pipeline.R') first.")
}

message(">> Loading deployment bundle...")
bundle <- load_bundle(BUNDLE_PATH)
message(">> Clocks: ", paste(bundle$meta$clocks_trained, collapse = ", "))

# ── 2. Load and combine cohort data ──────────────────────────────────────

message("\n>> Loading preclinical cohort data...")

f1 <- file.path(PROJECT_ROOT, "master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "master_PARTE2.csv")

if (!file.exists(f1)) stop("master_PARTE1.csv not found")

parte1 <- read.csv(f1, stringsAsFactors = FALSE)
message(">> PARTE 1: ", nrow(parte1), " rows")

if (file.exists(f2)) {
  parte2 <- read.csv(f2, stringsAsFactors = FALSE)
  message(">> PARTE 2: ", nrow(parte2), " rows")
  cohort <- rbind(parte1, parte2)
} else {
  cohort <- parte1
}

message(">> Combined cohort: ", nrow(cohort), " rows, ", ncol(cohort), " columns")

# ── 3. Transform biomarkers to NHANES conventions ────────────────────────

message("\n>> Transforming biomarkers to NHANES conventions...")

# Your lab data -> NHANES variable names
# PhenoAge original requires: albumin_gL, alp, lncrp, totchol,
#   lncreat_umol, hba1c, sbp, bun, uap, lymph, mcv, wbc

df <- cohort %>%
  mutate(
    patient_id = Protocolo,

    # albumin: your data is in g/dL, NHANES uses g/L (multiply by 10)
    # and also g/dL for some panels
    albumin    = as.numeric(albumin),
    albumin_gL = albumin * 10,

    # ALP: direct (U/L)
    alp   = as.numeric(ALP),
    lnalp = ifelse(is.na(alp) | alp <= 0, NA, log(alp)),

    # CRP: your data is raw mg/L, NHANES uses log(CRP)
    crp   = as.numeric(CRP),
    lncrp = ifelse(is.na(crp) | crp <= 0, NA, log(crp)),

    # Total cholesterol (mg/dL)
    totchol = as.numeric(total_cholesterol),

    # Creatinine: your data is mg/dL, NHANES uses log(creat in umol/L)
    # umol/L = mg/dL * 88.42
    creat       = as.numeric(creatinine),
    creat_umol  = creat * 88.42,
    lncreat     = ifelse(is.na(creat) | creat <= 0, NA, log(creat)),
    lncreat_umol = ifelse(is.na(creat_umol) | creat_umol <= 0, NA, log(creat_umol)),

    # HbA1c (%)
    hba1c   = as.numeric(HbA1c),
    lnhba1c = ifelse(is.na(hba1c) | hba1c <= 0, NA, log(hba1c)),

    # SBP: not available in your data — will limit some clocks
    sbp = NA_real_,

    # BUN (mg/dL)
    bun   = as.numeric(BUN),
    lnbun = ifelse(is.na(bun) | bun <= 0, NA, log(bun)),

    # Uric acid (mg/dL)
    uap   = as.numeric(uric_acid),
    lnuap = ifelse(is.na(uap) | uap <= 0, NA, log(uap)),

    # Lymphocytes (%) - your data is already percentage
    lymph = as.numeric(lymphocytes),

    # MCV (fL) - direct
    mcv = as.numeric(MCV),

    # WBC (K/uL) - your data is absolute count, NHANES uses K/uL
    # If your WBC is in cells/uL (e.g., 6850), divide by 1000
    wbc = ifelse(as.numeric(WBC) > 300, as.numeric(WBC) / 1000, as.numeric(WBC)),

    # ── Additional biomarkers for sub-clocks ──
    ggt       = as.numeric(GGT),
    insulin_raw = as.numeric(insulin),
    # BioAge uses "insulin" but it's fasting insulin (uU/mL)
    insulin   = as.numeric(insulin),

    trig      = as.numeric(triglycerides),

    rdw       = as.numeric(RDW),
    rbc       = as.numeric(RBC),
    # If RBC is in cells/uL (e.g., 5940000), convert to M/uL
    rbc       = ifelse(rbc > 100, rbc / 1e6, rbc),

    vitaminB12 = as.numeric(vitamin_B12),

    # Gender (NHANES: 1=male, 2=female)
    gender = ifelse(sex == "M", 1, 2),

    # Age
    age = as.numeric(age)
  )

message(">> Transformation complete")

# Show what we have
cat("\n--- Biomarker Availability ---\n")
pheno_vars <- c("albumin_gL", "alp", "lncrp", "totchol", "lncreat_umol",
                "hba1c", "sbp", "bun", "uap", "lymph", "mcv", "wbc")
for (v in pheno_vars) {
  n_ok <- sum(!is.na(df[[v]]))
  cat(sprintf("  %-16s %6d / %d  (%.1f%%)\n",
              v, n_ok, nrow(df), 100 * n_ok / nrow(df)))
}

# ── 4. Filter eligible subjects ──────────────────────────────────────────

message("\n>> Filtering eligible subjects...")

# Use the pre-computed eligibility flags
df$eligible_MAIN_PHENOAGE <- as.integer(df$eligible_MAIN_PHENOAGE)
df$eligible_HEMATOLOGIC   <- as.integer(df$eligible_HEMATOLOGIC)
df$eligible_B12_GLYCATION <- as.integer(df$eligible_B12_GLYCATION)
df$eligible_RENAL         <- as.integer(df$eligible_RENAL)

n_elig_pheno <- sum(df$eligible_MAIN_PHENOAGE == 1, na.rm = TRUE)
n_elig_hema  <- sum(df$eligible_HEMATOLOGIC == 1, na.rm = TRUE)
n_elig_b12   <- sum(df$eligible_B12_GLYCATION == 1, na.rm = TRUE)
n_elig_renal <- sum(df$eligible_RENAL == 1, na.rm = TRUE)

cat(sprintf("  Main PhenoAge eligible:    %d / %d  (%.1f%%)\n",
            n_elig_pheno, nrow(df), 100 * n_elig_pheno / nrow(df)))
cat(sprintf("  Hematologic eligible:      %d / %d  (%.1f%%)\n",
            n_elig_hema, nrow(df), 100 * n_elig_hema / nrow(df)))
cat(sprintf("  B12/Glycation eligible:    %d / %d  (%.1f%%)\n",
            n_elig_b12, nrow(df), 100 * n_elig_b12 / nrow(df)))
cat(sprintf("  Renal eligible:            %d / %d  (%.1f%%)\n",
            n_elig_renal, nrow(df), 100 * n_elig_renal / nrow(df)))

# ══════════════════════════════════════════════════════════════════════════
# 5. SCORE PHENOAGE (PRIORITY)
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 1: PhenoAge scoring (priority)")
message(strrep("=", 60))

# PhenoAge original uses these exact 12 biomarkers INCLUDING SBP.
# SBP (systolic blood pressure) is NOT available in our clinical lab data.
#
# Solution: We added a 'full_no_sbp' sub-clock to config.R that trains
# PhenoAge on the same 11 biomarkers minus SBP. This is the PRIMARY
# PhenoAge-like clock for clinical deployment.
#
# The canonical PhenoAge, KDM, and HD will return NA (expected behavior).

# Score using the bundle
scored <- score_patients(
  patient_data = df,
  bundle       = bundle,
  name_map     = NULL  # Already mapped above
)

n_phenoage <- sum(!is.na(scored$phenoage_orig))
message(">> PhenoAge original scored: ", n_phenoage, " / ", nrow(scored))

if (n_phenoage == 0) {
  message(">> PhenoAge original: all NA (expected — SBP not available)")
  message(">> The 'pheno_full_no_sbp' clock is the primary PhenoAge")
  message("   approximation for this cohort (11/12 biomarkers, no SBP).")
}

# Check the full_no_sbp clock specifically
n_full_no_sbp <- sum(!is.na(scored$pheno_full_no_sbp), na.rm = TRUE)
if (n_full_no_sbp > 0) {
  message(">> pheno_full_no_sbp scored: ", n_full_no_sbp, " / ", nrow(scored),
          " — THIS IS YOUR PRIMARY CLOCK")
} else {
  message(">> pheno_full_no_sbp: 0 scored. Re-run run_pipeline.R to train it.")
  message("   (The 'full_no_sbp' panel was just added to config.R)")
}

# ══════════════════════════════════════════════════════════════════════════
# 6. SCORE ALL AVAILABLE SUB-CLOCKS
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 2: Scoring all available clocks")
message(strrep("=", 60))

# Check which clocks produced non-NA values
clock_cols <- grep("_advance$", names(scored), value = TRUE)
cat("\n--- Clock Scoring Results ---\n")
for (cc in clock_cols) {
  ba_col <- sub("_advance$", "", cc)
  n_ok <- sum(!is.na(scored[[cc]]))
  cat(sprintf("  %-40s %6d / %d scored\n", ba_col, n_ok, nrow(scored)))
}

# ══════════════════════════════════════════════════════════════════════════
# 7. STATISTICAL SUMMARY — CORRELATION WITH AGE (PRIORITY)
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 3: Statistical summary — correlation with age")
message(strrep("=", 60))

# Merge scored data back with original cohort info
results <- cbind(
  df[, c("patient_id", "age", "sex", "source_file",
         "eligible_MAIN_PHENOAGE", "eligible_HEMATOLOGIC",
         "eligible_B12_GLYCATION", "eligible_RENAL"), drop = FALSE],
  scored[, setdiff(names(scored), c("patient_id", "age")), drop = FALSE]
)

# ── 7a. Cohort demographics ──
cat("\n--- Cohort Demographics ---\n")
cat(sprintf("  Total subjects:    %d\n", nrow(results)))
cat(sprintf("  PARTE 1:           %d\n", sum(results$source_file == "PARTE1", na.rm = TRUE)))
cat(sprintf("  PARTE 2:           %d\n", sum(results$source_file == "PARTE2", na.rm = TRUE)))
cat(sprintf("  Female:            %d (%.1f%%)\n",
            sum(results$sex == "F", na.rm = TRUE),
            100 * mean(results$sex == "F", na.rm = TRUE)))
cat(sprintf("  Male:              %d (%.1f%%)\n",
            sum(results$sex == "M", na.rm = TRUE),
            100 * mean(results$sex == "M", na.rm = TRUE)))
cat(sprintf("  Age range:         %.0f - %.0f years\n",
            min(results$age, na.rm = TRUE), max(results$age, na.rm = TRUE)))
cat(sprintf("  Age mean (SD):     %.1f (%.1f)\n",
            mean(results$age, na.rm = TRUE), sd(results$age, na.rm = TRUE)))
cat(sprintf("  Age median [IQR]:  %.0f [%.0f - %.0f]\n",
            median(results$age, na.rm = TRUE),
            quantile(results$age, 0.25, na.rm = TRUE),
            quantile(results$age, 0.75, na.rm = TRUE)))

# ── 7b. Correlation table: each BA clock vs chronological age ──
cat("\n--- Correlation of Biological Age with Chronological Age ---\n")
cat("(KEY validation: a good clock should correlate r > 0.8 with age)\n")
cat("(PRIMARY clock for this cohort: pheno_full_no_sbp)\n\n")

ba_cols <- grep("^phenoage_|^pheno_", names(results), value = TRUE)
ba_cols <- ba_cols[!grepl("_advance|_resid|_z$", ba_cols)]

cor_stats <- list()
for (bc in ba_cols) {
  vals <- results[[bc]]
  ages <- results$age
  complete <- !is.na(vals) & !is.na(ages)
  n_complete <- sum(complete)

  if (n_complete >= 10) {
    ct <- cor.test(ages[complete], vals[complete])
    lm_fit <- lm(vals[complete] ~ ages[complete])
    s <- summary(lm_fit)

    cor_stats[[bc]] <- data.frame(
      clock      = bc,
      n          = n_complete,
      r          = ct$estimate,
      r_squared  = s$r.squared,
      p_value    = ct$p.value,
      slope      = coef(lm_fit)[2],
      intercept  = coef(lm_fit)[1],
      RMSE       = sqrt(mean(s$residuals^2)),
      stringsAsFactors = FALSE
    )
  }
}

if (length(cor_stats) > 0) {
  cor_table <- do.call(rbind, cor_stats)
  rownames(cor_table) <- NULL
  cor_table <- cor_table[order(-abs(cor_table$r)), ]

  # Format for display
  disp <- cor_table
  disp$r <- sprintf("%.4f", disp$r)
  disp$r_squared <- sprintf("%.4f", disp$r_squared)
  disp$p_value <- ifelse(cor_table$p_value < 2.2e-16, "< 2.2e-16",
                         sprintf("%.2e", cor_table$p_value))
  disp$slope <- sprintf("%.4f", cor_table$slope)
  disp$intercept <- sprintf("%.2f", cor_table$intercept)
  disp$RMSE <- sprintf("%.2f", cor_table$RMSE)
  print(disp, row.names = FALSE)
} else {
  message(">> No clocks had enough data to compute correlations")
  cor_table <- NULL
}

# ── 7c. Advancement correlations ──
cat("\n--- Advancement Statistics (BA - CA) ---\n")
adv_cols <- grep("_advance$", names(results), value = TRUE)
adv_stats <- list()
for (ac in adv_cols) {
  vals <- results[[ac]]
  complete <- !is.na(vals)
  if (sum(complete) >= 10) {
    adv_stats[[ac]] <- data.frame(
      clock       = sub("_advance$", "", ac),
      n           = sum(complete),
      mean_adv    = mean(vals[complete]),
      sd_adv      = sd(vals[complete]),
      median_adv  = median(vals[complete]),
      pct_older   = 100 * mean(vals[complete] > 0),
      r_with_age  = cor(results$age[complete], vals[complete]),
      stringsAsFactors = FALSE
    )
  }
}
if (length(adv_stats) > 0) {
  adv_table <- do.call(rbind, adv_stats)
  rownames(adv_table) <- NULL
  print(adv_table, row.names = FALSE, digits = 3)
}

# ══════════════════════════════════════════════════════════════════════════
# 8. PLOTS (displayed in viewer)
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 4: Generating plots")
message(strrep("=", 60))

# ── Plot 1: BA vs Age scatter for each scored clock ──
for (bc in ba_cols) {
  vals <- results[[bc]]
  n_ok <- sum(!is.na(vals))
  if (n_ok < 30) next

  r_val <- cor(results$age, vals, use = "complete.obs")
  clean_name <- gsub("pheno_|phenoage_|_orig", "", bc)

  p <- ggplot(results[!is.na(vals), ], aes(x = age, y = .data[[bc]])) +
    geom_point(alpha = 0.08, size = 0.5, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "red", linewidth = 0.7) +
    geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Biological Age vs Chronological Age: ", clean_name),
         subtitle = sprintf("r = %.3f, n = %d", r_val, n_ok),
         x = "Chronological Age (years)",
         y = paste0(clean_name, " (Biological Age)"))
  print(p)
}

# ── Pre-compute which advancement columns have data (used by Plots 2 & 4) ──
valid_adv <- character(0)
if (length(adv_cols) > 0) {
  has_data <- vapply(adv_cols, function(x) sum(!is.na(results[[x]])) > 30, logical(1))
  valid_adv <- adv_cols[has_data]
}

# ── Plot 2: Advancement distribution ──
if (length(valid_adv) > 0) {
  scored_adv <- results[, c("age", valid_adv), drop = FALSE]

  if (length(valid_adv) > 0) {
    long_adv <- pivot_longer(scored_adv, cols = all_of(valid_adv),
                             names_to = "clock", values_to = "advancement")
    long_adv$clock <- gsub("pheno_|phenoage_|_advance|_orig", "", long_adv$clock)

    p_adv <- ggplot(long_adv[!is.na(long_adv$advancement), ],
                    aes(x = advancement)) +
      geom_histogram(aes(y = after_stat(density)), bins = 50,
                     fill = "steelblue", alpha = 0.6, colour = "white") +
      geom_density(colour = "navy", linewidth = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
      facet_wrap(~ clock, scales = "free") +
      theme_minimal(base_size = 12) +
      labs(title = "Distribution of Biological Age Advancement (BA - CA)",
           subtitle = "Values > 0 = biologically older than chronological age",
           x = "Advancement (years)", y = "Density")
    print(p_adv)
  }
}

# ── Plot 3: Age distribution of cohort ──
p_age <- ggplot(results, aes(x = age)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7, colour = "white") +
  theme_minimal(base_size = 13) +
  labs(title = "Age Distribution of Clinical Cohort",
       subtitle = sprintf("N = %d, Mean = %.1f, Median = %.0f",
                          nrow(results), mean(results$age, na.rm = TRUE),
                          median(results$age, na.rm = TRUE)),
       x = "Age (years)", y = "Count")
print(p_age)

# ── Plot 4: Advancement by age decade ──
if (length(valid_adv) > 0) {
  results$age_decade <- cut(results$age,
    breaks = c(0, 30, 40, 50, 60, 70, 80, 120),
    labels = c("<30", "30s", "40s", "50s", "60s", "70s", "80+"),
    right = FALSE)

  for (ac in valid_adv) {
    n_ok <- sum(!is.na(results[[ac]]))
    if (n_ok < 50) next
    clean <- gsub("pheno_|phenoage_|_advance|_orig", "", ac)

    p_box <- ggplot(results[!is.na(results[[ac]]), ],
                    aes(x = age_decade, y = .data[[ac]])) +
      geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.size = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      theme_minimal(base_size = 13) +
      labs(title = paste0("Advancement by Age Decade: ", clean),
           x = "Age Decade", y = "Advancement (years)")
    print(p_box)
  }
}

# ── Plot 5: Sex differences ──
for (bc in ba_cols) {
  n_ok <- sum(!is.na(results[[bc]]))
  if (n_ok < 50) next
  clean <- gsub("pheno_|phenoage_|_orig", "", bc)

  p_sex <- ggplot(results[!is.na(results[[bc]]), ],
                  aes(x = age, y = .data[[bc]], colour = sex)) +
    geom_point(alpha = 0.05, size = 0.4) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray40") +
    scale_colour_manual(values = c("F" = "firebrick", "M" = "steelblue")) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("BA vs Age by Sex: ", clean),
         x = "Chronological Age", y = "Biological Age",
         colour = "Sex")
  print(p_sex)
}

# ══════════════════════════════════════════════════════════════════════════
# 9. SAVE OUTPUTS
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 5: Saving outputs")
message(strrep("=", 60))

# Save scored data
write.csv(results, file.path(PROJECT_ROOT, "clinical_scored.csv"),
          row.names = FALSE)
message(">> Saved: clinical_scored.csv")

# Save correlation table
if (!is.null(cor_table)) {
  write.csv(cor_table, file.path(PROJECT_ROOT, "clinical_stats_summary.csv"),
            row.names = FALSE)
  message(">> Saved: clinical_stats_summary.csv")
}

# Save clinical report
tryCatch({
  report <- clinical_report_table(scored, bundle)
  write.csv(report, file.path(PROJECT_ROOT, "clinical_report.csv"),
            row.names = FALSE)
  message(">> Saved: clinical_report.csv")
}, error = function(e) {
  message("  ! Report table: ", e$message)
})

# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════

n_clocks <- sum(vapply(ba_cols, function(x) sum(!is.na(results[[x]])) > 0, logical(1)))

cat("\n")
cat("==========================================================\n")
cat("  CLINICAL SCORING COMPLETE\n")
cat("==========================================================\n")
cat("  Cohort:           ", nrow(results), " subjects\n")
cat("  PARTE 1:          ", sum(results$source_file == "PARTE1", na.rm = TRUE), "\n")
cat("  PARTE 2:          ", sum(results$source_file == "PARTE2", na.rm = TRUE), "\n")
cat("  Clocks scored:    ", n_clocks, "\n")
if (!is.null(cor_table) && nrow(cor_table) > 0) {
  best <- cor_table[1, ]
  cat("  Best r(BA, age):  ", sprintf("%.4f", best$r), " (", best$clock, ")\n")
}
cat("\n  IMPORTANT: SBP (blood pressure) is NOT available in this\n")
cat("  cohort. PhenoAge original, KDM, and HD require SBP and\n")
cat("  will show NA. Use 'pheno_full_no_sbp' as the primary\n")
cat("  PhenoAge-like clock (11/12 original biomarkers, no SBP).\n")
cat("\n  Output files:\n")
cat("    - clinical_scored.csv\n")
cat("    - clinical_stats_summary.csv\n")
cat("    - clinical_report.csv\n")
cat("==========================================================\n")

# Store for interactive use
clinical_results <- list(
  data      = results,
  scored    = scored,
  cor_table = cor_table,
  adv_table = if (exists("adv_table")) adv_table else NULL,
  bundle    = bundle
)

message("\n>> Results in 'clinical_results'")
message(">> Access: clinical_results$data, clinical_results$cor_table")
