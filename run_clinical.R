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

f1 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE2.csv")

if (!file.exists(f1)) stop("phenoage_master_PARTE1.csv not found")

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

# Levine 2018 PhenoAge uses these 9 biomarkers (+ age):
#   albumin (g/dL), alp (U/L), glucose (mg/dL), lncrp (log CRP),
#   lncreat (log creatinine mg/dL), lymph (%), mcv (fL), rdw (%),
#   wbc (1000 cells/uL)

df <- cohort %>%
  mutate(
    patient_id = Protocolo,

    # Age — filter invalid values (some may be year values like 2025)
    age = as.numeric(age),

    # ── Levine 2018 PhenoAge: 9 biomarkers ──

    # 1. Albumin (g/dL) — direct
    albumin = as.numeric(albumin),

    # 2. ALP (U/L) — direct (lowercase in new data)
    alp = as.numeric(alp),

    # 3. Glucose (mg/dL) — NOW AVAILABLE
    glucose = as.numeric(glucose),

    # 4. CRP: raw mg/L → log(CRP)
    #    New data has both 'crp' and 'log_crp'; use raw and compute log
    crp   = as.numeric(crp),
    lncrp = ifelse(is.na(crp) | crp <= 0, NA, log(crp)),

    # 5. Creatinine: mg/dL → log(creatinine mg/dL)
    creat   = as.numeric(creatinine),
    lncreat = ifelse(is.na(creat) | creat <= 0, NA, log(creat)),

    # 6. Lymphocyte % — column is 'lymphocyte' (not 'lymphocytes')
    lymph = as.numeric(lymphocyte),

    # 7. MCV (fL) — direct (lowercase in new data)
    mcv = as.numeric(mcv),

    # 8. RDW (%) — direct (lowercase in new data)
    rdw = as.numeric(rdw),

    # 9. WBC: if in cells/uL (>300), convert to 1000 cells/uL
    wbc = ifelse(as.numeric(wbc) > 300, as.numeric(wbc) / 1000, as.numeric(wbc)),

    # ── Additional biomarkers for sub-clocks ──

    # RBC: if in cells/uL (>100), convert to M/uL
    rbc = as.numeric(rbc),
    rbc = ifelse(rbc > 100, rbc / 1e6, rbc),

    # GGT (U/L) — direct (lowercase in new data)
    ggt = as.numeric(ggt),

    # Insulin (uU/mL)
    insulin = as.numeric(insulin),

    # Triglycerides (mg/dL)
    trig = as.numeric(triglycerides),

    # Total cholesterol (mg/dL) — column is 'cholesterol' in new data
    totchol = as.numeric(cholesterol),

    # HbA1c (%)
    hba1c = as.numeric(hba1c),

    # Vitamin B12 (pg/mL) — column is 'vitamin_b12' in new data
    vitaminB12 = as.numeric(vitamin_b12),

    # BUN (mg/dL)
    bun = as.numeric(bun),

    # Uric acid (mg/dL)
    uap = as.numeric(uric_acid),

    # Gender (NHANES: 1=male, 2=female)
    gender = ifelse(sex == "M", 1, 2)
  )

message(">> Transformation complete")

# ── 3b. Data quality: filter invalid ages ────────────────────────────────

n_before <- nrow(df)
n_bad_age <- sum(is.na(df$age) | df$age < 1 | df$age > 120, na.rm = TRUE)

if (n_bad_age > 0) {
  message(">> WARNING: ", n_bad_age, " rows with invalid age (NA, <1, or >120)")
  message("   Age range before filter: ", min(df$age, na.rm = TRUE), " - ",
          max(df$age, na.rm = TRUE))
  df <- df[!is.na(df$age) & df$age >= 1 & df$age <= 120, ]
  message("   After filter: ", nrow(df), " rows (removed ", n_before - nrow(df), ")")
}

cat(sprintf("  Age range: %.0f - %.0f years\n",
            min(df$age, na.rm = TRUE), max(df$age, na.rm = TRUE)))
cat(sprintf("  Age mean (SD): %.1f (%.1f)\n",
            mean(df$age, na.rm = TRUE), sd(df$age, na.rm = TRUE)))

# ══════════════════════════════════════════════════════════════════════════
# 4. BIOMARKER AVAILABILITY REPORT
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("BIOMARKER AVAILABILITY REPORT")
message(strrep("=", 60))

cat("\n--- Levine 2018 PhenoAge (9 biomarkers + age) ---\n")
levine_vars <- c("albumin", "alp", "glucose", "lncrp", "lncreat",
                  "lymph", "mcv", "rdw", "wbc")
for (v in levine_vars) {
  if (v %in% names(df)) {
    n_ok <- sum(!is.na(df[[v]]))
    cat(sprintf("  %-12s %6d / %d  (%5.1f%%)  %s\n",
                v, n_ok, nrow(df), 100 * n_ok / nrow(df),
                ifelse(n_ok == 0, "<< MISSING", "")))
  } else {
    cat(sprintf("  %-12s      0 / %d  (  0.0%%)  << COLUMN NOT IN DATA\n",
                v, nrow(df)))
  }
}

cat("\n--- Sub-clock biomarkers ---\n")
sub_vars <- c("rdw", "mcv", "rbc", "wbc", "lymph",  # hema
              "ggt", "insulin", "alp",                 # hepatic
              "trig", "totchol",                        # lipid
              "vitaminB12", "hba1c",                    # micronutrient
              "lncreat", "bun")                         # renal
sub_vars <- unique(sub_vars)
for (v in sub_vars) {
  if (v %in% names(df)) {
    n_ok <- sum(!is.na(df[[v]]))
    cat(sprintf("  %-14s %6d / %d  (%5.1f%%)\n",
                v, n_ok, nrow(df), 100 * n_ok / nrow(df)))
  }
}

# Which sub-clocks can score?
cat("\n--- Expected sub-clock coverage ---\n")
for (nm in names(PANELS_SUB)) {
  panel <- PANELS_SUB[[nm]]
  avail <- vapply(panel, function(v) {
    if (v %in% names(df)) sum(!is.na(df[[v]])) else 0L
  }, integer(1))
  # A row can be scored only if ALL biomarkers are non-NA
  complete_rows <- sum(complete.cases(df[, intersect(panel, names(df)), drop = FALSE]))
  cat(sprintf("  %-30s  %6d complete rows  (%5.1f%%)\n",
              nm, complete_rows, 100 * complete_rows / nrow(df)))
}

# ── 4b. Eligibility flags ──
cat("\n--- Pre-computed eligibility flags ---\n")
for (flag in c("eligible_MAIN_PHENOAGE", "eligible_HEMATOLOGIC",
               "eligible_B12_GLYCATION", "eligible_RENAL")) {
  if (flag %in% names(df)) {
    n_elig <- sum(as.integer(df[[flag]]) == 1, na.rm = TRUE)
    cat(sprintf("  %-30s %6d / %d  (%5.1f%%)\n",
                flag, n_elig, nrow(df), 100 * n_elig / nrow(df)))
  }
}

# ══════════════════════════════════════════════════════════════════════════
# 5. SCORE PHENOAGE (PRIORITY)
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 1: Scoring all clocks")
message(strrep("=", 60))

# All 9 Levine 2018 PhenoAge biomarkers are now available (including glucose).
# The published PhenoAge scoring will use hardcoded Levine 2018 coefficients.
# Bundle-based sub-clocks also scored for organ-system breakdown.

scored <- score_patients(
  patient_data = df,
  bundle       = bundle,
  name_map     = NULL
)

# ══════════════════════════════════════════════════════════════════════════
# 6. CHECK SCORING RESULTS
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 2: Scoring results summary")
message(strrep("=", 60))

clock_cols <- grep("_advance$", names(scored), value = TRUE)
cat("\n--- Clock Scoring Results ---\n")
for (cc in clock_cols) {
  ba_col <- sub("_advance$", "", cc)
  n_ok <- sum(!is.na(scored[[cc]]))
  cat(sprintf("  %-40s %6d / %d scored  (%5.1f%%)\n",
              ba_col, n_ok, nrow(scored), 100 * n_ok / nrow(scored)))
}

if (length(clock_cols) == 0) {
  message("\n>> WARNING: No clocks scored any patients!")
  message(">> This likely means the bundle was trained with different")
  message("   biomarker names than what's in the clinical data.")
  message(">> ACTION: Re-run source('run_pipeline.R') to retrain with")
  message("   the corrected config.R, then re-run this script.")
}

# ══════════════════════════════════════════════════════════════════════════
# 7. STATISTICAL SUMMARY — CORRELATION WITH AGE
# ══════════════════════════════════════════════════════════════════════════

message("\n", strrep("=", 60))
message("STEP 3: Statistical summary — correlation with age")
message(strrep("=", 60))

# Merge scored data back with original cohort info
results <- cbind(
  df[, c("patient_id", "age", "sex", "source_file"), drop = FALSE],
  scored[, setdiff(names(scored), c("patient_id", "age")), drop = FALSE]
)

# ── 7a. Cohort demographics ──
cat("\n--- Cohort Demographics ---\n")
cat(sprintf("  Total subjects:    %d\n", nrow(results)))
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
cat("(KEY validation: a good clock should correlate r > 0.8 with age)\n\n")

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

cor_table <- NULL
if (length(cor_stats) > 0) {
  cor_table <- do.call(rbind, cor_stats)
  rownames(cor_table) <- NULL
  cor_table <- cor_table[order(-abs(cor_table$r)), ]

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
}

# ── 7c. Advancement statistics ──
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
# 8. PLOTS
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

# ── Pre-compute valid advancement columns ──
valid_adv <- character(0)
if (length(adv_cols) > 0) {
  has_data <- vapply(adv_cols, function(x) sum(!is.na(results[[x]])) > 30, logical(1))
  valid_adv <- adv_cols[has_data]
}

# ── Plot 2: Advancement distribution ──
if (length(valid_adv) > 0) {
  long_adv <- pivot_longer(
    results[, c("age", valid_adv), drop = FALSE],
    cols = all_of(valid_adv),
    names_to = "clock", values_to = "advancement"
  )
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

# ── Plot 3: Age distribution ──
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

write.csv(results, file.path(PROJECT_ROOT, "clinical_scored.csv"),
          row.names = FALSE)
message(">> Saved: clinical_scored.csv")

if (!is.null(cor_table)) {
  write.csv(cor_table, file.path(PROJECT_ROOT, "clinical_stats_summary.csv"),
            row.names = FALSE)
  message(">> Saved: clinical_stats_summary.csv")
}

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

n_clocks <- length(clock_cols)

cat("\n")
cat("==========================================================\n")
cat("  CLINICAL SCORING COMPLETE\n")
cat("==========================================================\n")
cat("  Cohort:           ", nrow(results), " subjects\n")
cat("  Clocks scored:    ", n_clocks, "\n")
if (!is.null(cor_table) && nrow(cor_table) > 0) {
  best <- cor_table[1, ]
  cat("  Best r(BA, age):  ", sprintf("%.4f", best$r), " (", best$clock, ")\n")
}
cat("\n  NOTES:\n")
cat("  - All 9 Levine PhenoAge biomarkers available (incl. glucose).\n")
cat("  - SBP missing: KDM and HD will show NA (expected).\n")
cat("\n  IMPORTANT: If 0 clocks scored, re-run run_pipeline.R\n")
cat("  to retrain with the corrected config.R biomarker panels.\n")
cat("\n  Output files:\n")
cat("    - clinical_scored.csv\n")
cat("    - clinical_stats_summary.csv\n")
cat("    - clinical_report.csv\n")
cat("==========================================================\n")

clinical_results <- list(
  data      = results,
  scored    = scored,
  cor_table = cor_table,
  adv_table = if (exists("adv_table")) adv_table else NULL,
  bundle    = bundle
)

message("\n>> Results in 'clinical_results'")
message(">> Access: clinical_results$data, clinical_results$cor_table")
