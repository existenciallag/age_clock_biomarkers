###############################################################################
# run_recenter.R — Population recentering & QC comparison pipeline
#
# PURPOSE
# -------
# The clocks were trained on NHANES III and validated on NHANES IV (US pop).
# When applied to the Argentine clinical cohort, systematic offsets appear
# (all advancement values shift negative) due to population, assay, and
# temporal differences.
#
# This pipeline:
#   1. Loads the scored Argentine data (clinical_scored_lite.csv)
#   2. Loads the NHANES IV reference from the deployment bundle
#   3. Compares biomarker distributions (SMD, KS, Welch t)
#   4. Recenters advancement using THREE methods:
#      A) Within-cohort residualization: lm(BA ~ age) within Argentine data
#      B) Within-cohort residualization + sex: lm(BA ~ age + sex)
#      C) Age- and sex-stratified Z-scoring (non-parametric)
#   5. Computes dual-reference Z-scores (NHANES ref + local ref)
#   6. Validates recentered clocks via health-proxy correlations
#   7. Generates QC diagnostics and comparison plots
#
# STATISTICAL JUSTIFICATION
# -------------------------
# Recentering via within-cohort residualization is the standard approach
# used in NAKO (Germany, n=173K), UK Biobank, CHARLS (China), and is
# exactly what Levine defines as "PhenoAgeAccel" — the residual from
# regressing PhenoAge on chronological age within the study population.
#
# This does NOT alter model coefficients or Gompertz parameters.
# It preserves the relative ranking of patients (who ages faster/slower)
# while removing the calibration offset.
#
# REQUIREMENTS
# ------------
#   - clinical_scored_lite.csv (from run_clinical_lite.R)
#   - bioage_deployment_bundle.rds (from run_pipeline.R)
#   - phenoage_master_PARTE1.csv (+ optional PARTE2.csv)
#
# USAGE:  source("run_recenter.R")
###############################################################################

cat("==========================================================\n")
cat("  BIOLOGICAL AGE — POPULATION RECENTERING & QC\n")
cat("  Within-cohort recalibration | NHANES comparison\n")
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
})

source(file.path(PROJECT_ROOT, "R", "config.R"),       local = FALSE)
source(file.path(PROJECT_ROOT, "R", "export.R"),        local = FALSE)
source(file.path(PROJECT_ROOT, "R", "residualize.R"),   local = FALSE)
source(file.path(PROJECT_ROOT, "R", "qc_recenter.R"),   local = FALSE)

OUTPUT_DIR <- file.path(PROJECT_ROOT, "recenter_output")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────

scored_path <- file.path(PROJECT_ROOT, "clinical_scored_lite.csv")
bundle_path <- file.path(PROJECT_ROOT, "bioage_deployment_bundle.rds")

if (!file.exists(scored_path))
  stop("clinical_scored_lite.csv not found. Run source('run_clinical_lite.R') first.")
if (!file.exists(bundle_path))
  stop("bioage_deployment_bundle.rds not found. Run source('run_pipeline.R') first.")

message(">> Loading scored data...")
df <- read.csv(scored_path, stringsAsFactors = FALSE)
message(">> Subjects: ", nrow(df))

message(">> Loading deployment bundle...")
bundle <- load_bundle(bundle_path)

# Load raw clinical data for biomarker-level comparison
message(">> Loading raw clinical data for biomarker comparison...")
f1 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE2.csv")
raw_clin <- read.csv(f1, stringsAsFactors = FALSE)
if (file.exists(f2)) raw_clin <- rbind(raw_clin, read.csv(f2, stringsAsFactors = FALSE))

# Transform raw clinical biomarkers to NHANES naming (same as run_clinical_lite.R)
raw_clin <- raw_clin %>%
  mutate(
    age     = as.numeric(age),
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
    rbc     = as.numeric(rbc),
    rbc     = ifelse(rbc > 100, rbc / 1e6, rbc),
    ggt     = as.numeric(ggt),
    insulin    = as.numeric(insulin),
    trig       = as.numeric(triglycerides),
    totchol    = as.numeric(cholesterol),
    hba1c      = as.numeric(hba1c),
    vitaminB12 = as.numeric(vitamin_b12),
    bun        = as.numeric(bun),
    uap        = as.numeric(uric_acid)
  )
raw_clin <- raw_clin[!is.na(raw_clin$age) & raw_clin$age >= 1 & raw_clin$age <= 120, ]

# ── 2. Extract NHANES IV reference biomarkers from BioAge ─────────────────

message("\n>> Extracting NHANES IV reference data from BioAge package...")

# Reconstruct NHANES IV data by calling phenoage_nhanes on the Levine panel.
# The $data component contains the NHANES IV projection with raw biomarkers.
nhanes_ref <- tryCatch({
  res <- phenoage_nhanes(BIOMARKERS_PHENO_LEVINE)
  res$data
}, error = function(e) {
  message(">> Could not extract NHANES IV data: ", e$message)
  NULL
})

if (!is.null(nhanes_ref)) {
  message(">> NHANES IV reference: ", nrow(nhanes_ref), " subjects")
} else {
  message(">> WARNING: NHANES IV reference not available; skipping biomarker comparison")
}

# Identify advancement columns
adv_cols <- grep("_advance$", names(df), value = TRUE)
ba_cols  <- sub("_advance$", "", adv_cols)
message(">> Clocks to recenter: ", paste(ba_cols, collapse = ", "))

###############################################################################
#   SECTION 1: BIOMARKER DISTRIBUTION COMPARISON
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 1: BIOMARKER DISTRIBUTION COMPARISON          #\n")
cat("##########################################################\n\n")

# Biomarkers shared between Levine panel and sub-clocks
comparison_vars <- c("albumin", "alp", "glucose", "lncrp", "lncreat",
                     "lymph", "mcv", "rdw", "wbc",
                     "rbc", "trig", "totchol", "hba1c", "bun")

# Argentine biomarker summary
arg_bm_summary <- biomarker_summary(raw_clin, comparison_vars, label = "Argentina")

if (!is.null(nhanes_ref)) {
  # NHANES biomarker summary
  nhanes_bm_summary <- biomarker_summary(nhanes_ref, comparison_vars, label = "NHANES_IV")

  # Compare
  bm_comparison <- compare_biomarkers(
    data_a  = nhanes_ref,
    data_b  = raw_clin,
    vars    = comparison_vars,
    label_a = "NHANES_IV",
    label_b = "Argentina"
  )

  if (!is.null(bm_comparison) && nrow(bm_comparison) > 0) {
    cat("--- Biomarker Comparison: Argentina vs NHANES IV ---\n\n")
    print(bm_comparison[, c("variable", "n_ref", "mean_ref", "sd_ref",
                            "n_target", "mean_target", "sd_target",
                            "SMD", "smd_flag")], row.names = FALSE)

    write.csv(bm_comparison,
              file.path(OUTPUT_DIR, "biomarker_comparison_nhanes.csv"),
              row.names = FALSE)

    # SMD bar chart
    p_smd <- plot_smd_chart(bm_comparison)
    ggsave(file.path(OUTPUT_DIR, "smd_biomarker_comparison.png"),
           p_smd, width = 10, height = 7, dpi = 150)
    message(">> Saved: smd_biomarker_comparison.png")

    # Flag variables with large SMD
    large_smd <- bm_comparison$variable[abs(bm_comparison$SMD) > 0.5]
    if (length(large_smd) > 0) {
      cat("\n>> WARNING: Large distributional differences (|SMD| > 0.5):\n")
      cat("   ", paste(large_smd, collapse = ", "), "\n")
      cat("   These biomarkers may contribute to systematic clock offset.\n")
    }
  }

  # Combined summary table
  combined_bm <- rbind(nhanes_bm_summary, arg_bm_summary)
  write.csv(combined_bm,
            file.path(OUTPUT_DIR, "biomarker_summary_combined.csv"),
            row.names = FALSE)
} else {
  cat(">> NHANES IV data unavailable. Reporting Argentine-only summary.\n\n")
  print(arg_bm_summary, row.names = FALSE)
  write.csv(arg_bm_summary,
            file.path(OUTPUT_DIR, "biomarker_summary_argentina.csv"),
            row.names = FALSE)
}

###############################################################################
#   SECTION 2: CLOCK QC — BEFORE RECENTERING
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 2: CLOCK QC — BEFORE RECENTERING              #\n")
cat("##########################################################\n\n")

pre_qc <- clinical_qc_table(df, adv_cols)
cat("--- Pre-Recentering Clock QC ---\n\n")
print(pre_qc, row.names = FALSE)

write.csv(pre_qc, file.path(OUTPUT_DIR, "qc_pre_recenter.csv"), row.names = FALSE)

# Identify which clocks have enough data (>= 100 subjects)
viable_clocks <- pre_qc$clock[pre_qc$n_scored >= 100]
viable_adv    <- paste0(viable_clocks, "_advance")
viable_adv    <- intersect(viable_adv, adv_cols)
message("\n>> Viable clocks (n >= 100): ", paste(viable_clocks, collapse = ", "))

if (length(viable_adv) == 0)
  stop("No clocks have >= 100 scored subjects. Cannot proceed.")

###############################################################################
#   SECTION 3: RECENTERING — THREE METHODS
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 3: RECENTERING (3 METHODS)                    #\n")
cat("##########################################################\n\n")

# Ensure sex is available as factor
df$sex <- as.character(df$sex)

# ── Method A: Within-cohort residualization (lm(BA ~ age)) ────────────────

cat("--- Method A: Within-Cohort Residualization (BA ~ age) ---\n\n")

for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  ba_vals <- df[[ba_col]]
  ok <- !is.na(ba_vals) & is.finite(ba_vals) & !is.na(df$age)

  if (sum(ok) < 30) {
    message(">> ", ba_col, ": too few observations, skipping method A")
    next
  }

  # Fit local regression: BA ~ age
  fit_local <- lm(ba_vals[ok] ~ df$age[ok])
  coefs <- coef(fit_local)

  # Compute recentered advancement as residual
  resid_col  <- paste0(ba_col, "_local_resid")
  z_col      <- paste0(ba_col, "_local_z")

  df[[resid_col]] <- NA_real_
  df[[resid_col]][ok] <- residuals(fit_local)
  df[[z_col]] <- NA_real_
  valid_resid <- !is.na(df[[resid_col]])
  if (sum(valid_resid) > 1) {
    df[[z_col]][valid_resid] <- as.numeric(scale(df[[resid_col]][valid_resid]))
  }

  cat(sprintf("  %-40s intercept = %+.3f  slope = %.4f  n = %d\n",
              ba_col, coefs[1], coefs[2], sum(ok)))
  cat(sprintf("  %-40s mean_resid = %+.4f  SD = %.3f\n",
              "", mean(df[[resid_col]], na.rm = TRUE),
              sd(df[[resid_col]], na.rm = TRUE)))
}

# ── Method B: Within-cohort residualization (BA ~ age + sex) ──────────────

cat("\n--- Method B: Within-Cohort Residualization (BA ~ age + sex) ---\n\n")

for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  ba_vals <- df[[ba_col]]
  ok <- !is.na(ba_vals) & is.finite(ba_vals) &
        !is.na(df$age) & df$sex %in% c("F", "M")

  if (sum(ok) < 30) {
    message(">> ", ba_col, ": too few observations, skipping method B")
    next
  }

  fit_sex <- lm(ba_vals[ok] ~ df$age[ok] + df$sex[ok])
  coefs <- coef(fit_sex)

  resid_col <- paste0(ba_col, "_sexadj_resid")
  z_col     <- paste0(ba_col, "_sexadj_z")

  df[[resid_col]] <- NA_real_
  df[[resid_col]][ok] <- residuals(fit_sex)
  df[[z_col]] <- NA_real_
  valid_resid <- !is.na(df[[resid_col]])
  if (sum(valid_resid) > 1) {
    df[[z_col]][valid_resid] <- as.numeric(scale(df[[resid_col]][valid_resid]))
  }

  cat(sprintf("  %-40s intercept = %+.3f  slope_age = %.4f  slope_sexM = %+.3f  n = %d\n",
              ba_col, coefs[1], coefs[2],
              ifelse(length(coefs) >= 3, coefs[3], NA), sum(ok)))
}

# ── Method C: Age-stratified Z-scoring (non-parametric) ───────────────────

cat("\n--- Method C: Age-Stratified Z-Scoring (Non-Parametric) ---\n\n")

# Create age strata
df$age_stratum <- cut(df$age,
  breaks = c(0, 30, 40, 50, 60, 70, 80, 120),
  labels = c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
  right = FALSE
)

for (ac in viable_adv) {
  ba_col  <- sub("_advance$", "", ac)
  adv_vals <- df[[ac]]
  strat_z_col <- paste0(ba_col, "_strat_z")
  df[[strat_z_col]] <- NA_real_

  strata_summary <- list()

  for (stratum in levels(df$age_stratum)) {
    idx <- df$age_stratum == stratum & !is.na(adv_vals) & is.finite(adv_vals)
    n_s <- sum(idx, na.rm = TRUE)
    if (n_s < 10) next

    x_s <- adv_vals[idx]
    mu_s  <- mean(x_s)
    sd_s  <- sd(x_s)

    if (sd_s > 0) {
      df[[strat_z_col]][idx] <- (x_s - mu_s) / sd_s
    }

    strata_summary[[stratum]] <- data.frame(
      clock = ba_col, age_stratum = stratum, n = n_s,
      mean_adv = round(mu_s, 3), sd_adv = round(sd_s, 3),
      stringsAsFactors = FALSE
    )
  }

  if (length(strata_summary) > 0) {
    st <- do.call(rbind, strata_summary)
    cat(sprintf("  %s:\n", ba_col))
    for (i in seq_len(nrow(st))) {
      cat(sprintf("    %-8s n = %5d  mean = %+.3f  SD = %.3f\n",
                  st$age_stratum[i], st$n[i], st$mean_adv[i], st$sd_adv[i]))
    }
  }
}

# ── Also compute NHANES-referenced Z-scores (dual reference) ─────────────

cat("\n--- Dual Reference: NHANES Z-scores (for international comparability) ---\n\n")

if (!is.null(bundle$pop_stats) && !is.null(bundle$resid_coefs)) {
  for (ac in viable_adv) {
    ps <- bundle$pop_stats[bundle$pop_stats$variable == ac, ]
    rc <- bundle$resid_coefs[bundle$resid_coefs$variable == ac, ]

    nhanes_z_col <- paste0(sub("_advance$", "", ac), "_nhanes_z")
    df[[nhanes_z_col]] <- NA_real_

    if (nrow(rc) == 1 && nrow(ps) == 1) {
      ok <- !is.na(df[[ac]]) & is.finite(df[[ac]]) & !is.na(df$age)
      predicted <- rc$intercept + rc$slope * df$age[ok]
      resid_nhanes <- df[[ac]][ok] - predicted
      df[[nhanes_z_col]][ok] <- (resid_nhanes - ps$mean) / ps$sd

      cat(sprintf("  %-40s NHANES ref: mu = %.3f, sd = %.3f\n",
                  sub("_advance$", "", ac), ps$mean, ps$sd))
    } else {
      cat(sprintf("  %-40s no NHANES reference available\n",
                  sub("_advance$", "", ac)))
    }
  }
} else {
  message(">> Bundle pop_stats/resid_coefs not available; skipping NHANES Z-scores")
}

###############################################################################
#   SECTION 4: RECENTERING QC — COMPARISON
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 4: RECENTERING QC — BEFORE vs AFTER           #\n")
cat("##########################################################\n\n")

# Compare original vs Method A (local residualization)
local_resid_cols <- paste0(sub("_advance$", "", viable_adv), "_local_resid")
local_resid_cols <- intersect(local_resid_cols, names(df))

if (length(local_resid_cols) > 0) {
  # Match original advancement columns to their local resid counterparts
  orig_for_compare <- paste0(sub("_local_resid$", "", local_resid_cols), "_advance")
  orig_for_compare <- intersect(orig_for_compare, names(df))

  if (length(orig_for_compare) == length(local_resid_cols)) {
    recenter_comparison <- compare_recentering(df, orig_for_compare, local_resid_cols)

    if (!is.null(recenter_comparison) && nrow(recenter_comparison) > 0) {
      cat("--- Before vs After Recentering (Method A) ---\n\n")
      print(recenter_comparison, row.names = FALSE)
      write.csv(recenter_comparison,
                file.path(OUTPUT_DIR, "recenter_comparison.csv"),
                row.names = FALSE)
    }
  }
}

# Density overlay plots for each viable clock
for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  local_col <- paste0(ba_col, "_local_resid")

  if (local_col %in% names(df)) {
    p <- plot_recenter_density(df, ac, local_col)
    fname <- paste0("recenter_density_", gsub("pheno_|phenoage_", "", ba_col), ".png")
    ggsave(file.path(OUTPUT_DIR, fname), p, width = 9, height = 5, dpi = 150)
    message(">> Saved: ", fname)
  }
}

###############################################################################
#   SECTION 5: POST-RECENTERING CLOCK QC
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 5: POST-RECENTERING CLOCK QC                  #\n")
cat("##########################################################\n\n")

# QC for local Z-scores
local_z_cols <- paste0(sub("_advance$", "", viable_adv), "_local_z")
local_z_cols <- intersect(local_z_cols, names(df))

if (length(local_z_cols) > 0) {
  post_qc_stats <- lapply(local_z_cols, function(zc) {
    x <- df[[zc]]
    ok <- !is.na(x) & is.finite(x)
    clock <- sub("_local_z$", "", zc)
    data.frame(
      clock      = clock,
      method     = "local_z",
      n          = sum(ok),
      mean       = round(mean(x[ok]), 4),
      sd         = round(sd(x[ok]), 4),
      median     = round(median(x[ok]), 4),
      pct_above0 = round(100 * mean(x[ok] > 0), 1),
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  post_qc <- do.call(rbind, post_qc_stats)
  cat("--- Post-Recentering Z-Score QC (Method A) ---\n")
  cat("  (Expected: mean ~ 0, SD ~ 1, pct_above0 ~ 50%)\n\n")
  print(post_qc, row.names = FALSE)

  write.csv(post_qc, file.path(OUTPUT_DIR, "qc_post_recenter.csv"), row.names = FALSE)
}

# QQ plots for recentered advancement
for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  local_col <- paste0(ba_col, "_local_resid")

  if (local_col %in% names(df)) {
    x <- df[[local_col]]
    x <- x[!is.na(x) & is.finite(x)]
    if (length(x) > 50) {
      p_qq <- plot_advancement_qq(x, title = paste0("QQ: ", ba_col, " (recentered)"))
      fname <- paste0("qq_recentered_", gsub("pheno_|phenoage_", "", ba_col), ".png")
      ggsave(file.path(OUTPUT_DIR, fname), p_qq, width = 7, height = 6, dpi = 150)
    }
  }
}

###############################################################################
#   SECTION 6: HEALTH PROXY VALIDATION
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 6: HEALTH PROXY VALIDATION                    #\n")
cat("##########################################################\n\n")

cat("  Without survival data, we validate recentered clocks against\n")
cat("  available health proxies from the clinical dataset.\n\n")

# Load raw clinical data for extra biomarkers as health proxies
health_proxies <- c(
  "hemoglobin", "hematocrit", "platelets", "mch", "mchc",
  "hdl", "ldl", "non_hdl", "lp_a",
  "ast", "alt", "ldh", "ferritin", "transferrin", "transferrin_saturation",
  "fibrinogen", "homocysteine", "TSH", "free_t4", "vitamin_d",
  "neutrophils", "eosinophils", "basophils", "monocytes", "eGFR"
)

# Merge health proxies into df
raw_clin$patient_id <- raw_clin$Protocolo
hp_available <- intersect(health_proxies, names(raw_clin))
if (length(hp_available) > 0 && "patient_id" %in% names(df)) {
  for (hp in hp_available) raw_clin[[hp]] <- as.numeric(raw_clin[[hp]])
  merge_cols <- c("patient_id", hp_available)
  merge_cols <- intersect(merge_cols, names(raw_clin))
  df <- merge(df, raw_clin[, merge_cols, drop = FALSE],
              by = "patient_id", all.x = TRUE, suffixes = c("", ".hp"))

  cat(sprintf(">> Health proxies available: %d\n", length(hp_available)))
  cat("   ", paste(hp_available, collapse = ", "), "\n\n")

  # Spearman correlations: local_z vs health proxies
  proxy_results <- list()
  for (zc in local_z_cols) {
    clock <- sub("_local_z$", "", zc)
    for (hp in hp_available) {
      ok <- !is.na(df[[zc]]) & is.finite(df[[zc]]) &
            !is.na(df[[hp]]) & is.finite(df[[hp]])
      n_ok <- sum(ok)
      if (n_ok < 30) next

      ct <- cor.test(df[[zc]][ok], df[[hp]][ok], method = "spearman")
      proxy_results[[length(proxy_results) + 1]] <- data.frame(
        clock   = clock,
        proxy   = hp,
        n       = n_ok,
        rho     = round(ct$estimate, 4),
        p_value = ct$p.value,
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }

  if (length(proxy_results) > 0) {
    proxy_table <- do.call(rbind, proxy_results)
    proxy_table$p_adj <- p.adjust(proxy_table$p_value, method = "BH")
    proxy_table$sig <- ifelse(proxy_table$p_adj < 0.001, "***",
                       ifelse(proxy_table$p_adj < 0.01, "**",
                       ifelse(proxy_table$p_adj < 0.05, "*", "")))

    # Show top correlations per clock
    cat("--- Top Health Proxy Correlations (recentered clocks) ---\n\n")
    for (clk in unique(proxy_table$clock)) {
      ct_sub <- proxy_table[proxy_table$clock == clk, ]
      ct_sub <- ct_sub[order(-abs(ct_sub$rho)), ]
      top_n  <- min(10, nrow(ct_sub))
      cat(sprintf("  %s (top %d):\n", clk, top_n))
      for (i in seq_len(top_n)) {
        cat(sprintf("    %-30s rho = %+.4f  p_adj = %.2e  %s  (n = %d)\n",
                    ct_sub$proxy[i], ct_sub$rho[i], ct_sub$p_adj[i],
                    ct_sub$sig[i], ct_sub$n[i]))
      }
      cat("\n")
    }

    write.csv(proxy_table, file.path(OUTPUT_DIR, "health_proxy_validation.csv"),
              row.names = FALSE)
    message(">> Saved: health_proxy_validation.csv")

    # Summarise: for each clock, count how many proxies are significantly correlated
    validation_summary <- proxy_table %>%
      group_by(clock) %>%
      summarise(
        n_proxies_tested = n(),
        n_sig_005        = sum(p_adj < 0.05),
        n_sig_001        = sum(p_adj < 0.01),
        median_abs_rho   = round(median(abs(rho)), 4),
        max_abs_rho      = round(max(abs(rho)), 4),
        .groups = "drop"
      )
    cat("--- Validation Summary ---\n\n")
    print(as.data.frame(validation_summary), row.names = FALSE)
    write.csv(validation_summary,
              file.path(OUTPUT_DIR, "validation_summary.csv"), row.names = FALSE)
  }
} else {
  cat(">> No health proxy biomarkers available for validation.\n")
}

###############################################################################
#   SECTION 7: METHOD COMPARISON — ALL THREE RECENTERING APPROACHES
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 7: METHOD COMPARISON                          #\n")
cat("##########################################################\n\n")

method_comparison <- list()

for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  cols_to_check <- list(
    original  = ac,
    local     = paste0(ba_col, "_local_z"),
    sex_adj   = paste0(ba_col, "_sexadj_z"),
    strat     = paste0(ba_col, "_strat_z"),
    nhanes    = paste0(ba_col, "_nhanes_z")
  )

  for (method_name in names(cols_to_check)) {
    col <- cols_to_check[[method_name]]
    if (!col %in% names(df)) next
    x <- df[[col]]
    ok <- !is.na(x) & is.finite(x)
    if (sum(ok) < 30) next

    method_comparison[[length(method_comparison) + 1]] <- data.frame(
      clock  = ba_col,
      method = method_name,
      n      = sum(ok),
      mean   = round(mean(x[ok]), 4),
      sd     = round(sd(x[ok]), 4),
      median = round(median(x[ok]), 4),
      pct_pos = round(100 * mean(x[ok] > 0), 1),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }
}

if (length(method_comparison) > 0) {
  mc_table <- do.call(rbind, method_comparison)
  cat("--- Method Comparison Summary ---\n\n")
  print(mc_table, row.names = FALSE)
  write.csv(mc_table, file.path(OUTPUT_DIR, "method_comparison.csv"), row.names = FALSE)
}

###############################################################################
#   SECTION 8: INTER-CLOCK CORRELATION (RECENTERED)
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 8: INTER-CLOCK CORRELATION (RECENTERED)       #\n")
cat("##########################################################\n\n")

if (length(local_z_cols) >= 2) {
  cor_vars <- local_z_cols
  cor_data <- df[, cor_vars, drop = FALSE]
  for (v in cor_vars) cor_data[[v]] <- as.numeric(cor_data[[v]])

  cor_mat <- cor(cor_data, use = "pairwise.complete.obs", method = "spearman")
  colnames(cor_mat) <- gsub("_local_z$", "", colnames(cor_mat))
  rownames(cor_mat) <- gsub("_local_z$", "", rownames(cor_mat))

  cat("--- Spearman Correlation Matrix (Recentered Clocks) ---\n\n")
  print(round(cor_mat, 3))

  write.csv(as.data.frame(cor_mat),
            file.path(OUTPUT_DIR, "interclock_correlation_recentered.csv"))

  # Heatmap
  if (requireNamespace("corrplot", quietly = TRUE)) {
    png(file.path(OUTPUT_DIR, "interclock_heatmap_recentered.png"),
        width = 800, height = 700, res = 120)
    corrplot::corrplot(cor_mat, method = "color", type = "upper",
                       tl.col = "black", tl.cex = 0.9,
                       addCoef.col = "black", number.cex = 0.8,
                       title = "Inter-Clock Correlation (Recentered)",
                       mar = c(0, 0, 2, 0))
    dev.off()
    message(">> Saved: interclock_heatmap_recentered.png")
  }
}

###############################################################################
#   SECTION 9: SAVE RECENTERED DATASET
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 9: SAVE RECENTERED DATASET                    #\n")
cat("##########################################################\n\n")

# Select output columns
out_cols <- c("patient_id", "age", "sex")

for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  out_cols <- c(out_cols,
    ba_col,                             # Raw biological age
    ac,                                 # Raw advancement (NHANES-ref)
    paste0(ba_col, "_local_resid"),     # Method A: local residual
    paste0(ba_col, "_local_z"),         # Method A: local Z-score
    paste0(ba_col, "_sexadj_resid"),    # Method B: sex-adjusted residual
    paste0(ba_col, "_sexadj_z"),        # Method B: sex-adjusted Z-score
    paste0(ba_col, "_strat_z"),         # Method C: age-stratified Z-score
    paste0(ba_col, "_nhanes_z")         # NHANES-referenced Z-score
  )
}

out_cols <- intersect(out_cols, names(df))
out_df <- df[, out_cols, drop = FALSE]

write.csv(out_df,
          file.path(OUTPUT_DIR, "clinical_recentered.csv"),
          row.names = FALSE)
message(">> Saved: clinical_recentered.csv (", ncol(out_df), " columns)")

# Save local regression coefficients for reproducibility
local_coefs <- list()
for (ac in viable_adv) {
  ba_col <- sub("_advance$", "", ac)
  ba_vals <- df[[ba_col]]
  ok <- !is.na(ba_vals) & is.finite(ba_vals) & !is.na(df$age)
  if (sum(ok) < 30) next

  fit <- lm(ba_vals[ok] ~ df$age[ok])
  local_coefs[[ba_col]] <- data.frame(
    clock     = ba_col,
    intercept = coef(fit)[1],
    slope     = coef(fit)[2],
    n         = sum(ok),
    mean_resid = mean(residuals(fit)),
    sd_resid   = sd(residuals(fit)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}
if (length(local_coefs) > 0) {
  coef_table <- do.call(rbind, local_coefs)
  write.csv(coef_table,
            file.path(OUTPUT_DIR, "local_regression_coefficients.csv"),
            row.names = FALSE)
  message(">> Saved: local_regression_coefficients.csv")
}

###############################################################################
#   FINAL SUMMARY
###############################################################################

cat("\n==========================================================\n")
cat("  RECENTERING PIPELINE COMPLETE\n")
cat("==========================================================\n")
cat("  Subjects:             ", format(nrow(df), big.mark = ","), "\n")
cat("  Clocks recentered:    ", length(viable_clocks), "\n")
cat("  Methods applied:      A (age-residual), B (age+sex), C (stratified)\n")
if (!is.null(nhanes_ref))
  cat("  NHANES comparison:    YES (", nrow(nhanes_ref), " reference subjects)\n")
cat("  Output directory:     ", OUTPUT_DIR, "\n")
cat("\n  Files saved:\n")
saved_files <- list.files(OUTPUT_DIR, full.names = FALSE)
for (f in saved_files) cat("    - ", f, "\n")

cat("\n  KEY OUTPUT COLUMNS (in clinical_recentered.csv):\n")
cat("    *_advance         Raw advancement (NHANES-referenced)\n")
cat("    *_local_resid     Method A: within-cohort residual\n")
cat("    *_local_z         Method A: within-cohort Z-score (RECOMMENDED)\n")
cat("    *_sexadj_z        Method B: sex-adjusted Z-score\n")
cat("    *_strat_z         Method C: age-stratified Z-score\n")
cat("    *_nhanes_z        NHANES-referenced Z-score (for comparability)\n")
cat("==========================================================\n")

message("\n>> Done. All outputs in: ", OUTPUT_DIR)
