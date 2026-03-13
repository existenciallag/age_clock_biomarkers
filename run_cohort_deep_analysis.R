###############################################################################
# run_cohort_deep_analysis.R — Age-stratified, sex-difference, outlier, and
#                               non-clock biomarker analysis
#
# Depends on: clinical_scored_lite.csv (from run_clinical_lite.R)
#             phenoage_master_PARTE1.csv (+ optional PARTE2.csv)
#
# This script does NOT modify run_clinical_lite.R.  It reads its output
# and the raw clinical data to perform deeper exploratory analyses:
#
#   1. Age-stratified biological age advancement
#   2. Sex differences in biological aging
#   3. Outlier analysis (fastest / slowest agers)
#   4. Non-clock biomarker correlations with BA advancement
#   5. Multi-biomarker tendency analysis
#
# Usage:  source("run_cohort_deep_analysis.R")
###############################################################################

cat("==========================================================\n")
cat("  COHORT DEEP ANALYSIS\n")
cat("  Age Stratification | Sex Differences | Outliers |\n")
cat("  Non-Clock Biomarker Correlations\n")
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
  library(corrplot)
})

# Try to load corrplot; if unavailable, define a fallback
if (!requireNamespace("corrplot", quietly = TRUE)) {
  message(">> Note: 'corrplot' not installed. Correlation matrices will use base heatmap.")
  USE_CORRPLOT <- FALSE
} else {
  USE_CORRPLOT <- TRUE
}

OUTPUT_DIR <- file.path(PROJECT_ROOT, "deep_analysis_output")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ── 1. Load scored data + raw clinical data ───────────────────────────────

scored_path <- file.path(PROJECT_ROOT, "clinical_scored_lite.csv")
if (!file.exists(scored_path)) {
  stop("clinical_scored_lite.csv not found. Run source('run_clinical_lite.R') first.")
}

message(">> Loading scored data...")
scored <- read.csv(scored_path, stringsAsFactors = FALSE)
message(">> Scored subjects: ", nrow(scored))

# Load raw clinical for non-clock biomarkers
message(">> Loading raw clinical data for extra biomarkers...")
f1 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE2.csv")
raw <- read.csv(f1, stringsAsFactors = FALSE)
if (file.exists(f2)) raw <- rbind(raw, read.csv(f2, stringsAsFactors = FALSE))

# Merge extra biomarkers into scored
# Identify non-clock biomarkers (not in any PhenoAge panel)
clock_raw_names <- c("albumin", "creatinine", "glucose", "crp", "log_crp",
                     "lymphocyte", "mcv", "rdw", "wbc", "alp",
                     "rbc", "ggt", "insulin", "triglycerides", "cholesterol",
                     "vitamin_b12", "hba1c", "bun", "uric_acid")

extra_biomarkers <- c(
  "hemoglobin", "hematocrit", "platelets", "mch", "mchc",
  "hdl", "ldl", "non_hdl", "lp_a",
  "ast", "alt", "ldh", "ferritin", "transferrin", "transferrin_saturation",
  "fibrinogen", "homocysteine",
  "TSH", "free_t4", "vitamin_d",
  "neutrophils", "eosinophils", "basophils", "monocytes",
  "eGFR"
)

# Keep only those that exist in raw data
extra_biomarkers <- intersect(extra_biomarkers, names(raw))
message(">> Extra (non-clock) biomarkers found: ", length(extra_biomarkers))
message("   ", paste(extra_biomarkers, collapse = ", "))

# Add Protocolo to scored for merging
raw$patient_id <- raw$Protocolo
merge_cols <- c("patient_id", extra_biomarkers)
merge_cols <- intersect(merge_cols, names(raw))

# Convert extra biomarkers to numeric
for (eb in extra_biomarkers) {
  raw[[eb]] <- as.numeric(raw[[eb]])
}

df <- merge(scored, raw[, merge_cols, drop = FALSE],
            by = "patient_id", all.x = TRUE, suffixes = c("", ".raw"))

# ── Identify available advancement columns ────────────────────────────────
adv_cols <- grep("_advance$", names(df), value = TRUE)
# Pick the primary advancement column (phenoage_orig preferred)
primary_adv <- if ("phenoage_orig_advance" %in% adv_cols) {
  "phenoage_orig_advance"
} else {
  adv_cols[1]
}

if (is.null(primary_adv) || length(adv_cols) == 0) {
  stop("No advancement columns found. Ensure run_clinical_lite.R scored at least one clock.")
}
message(">> Primary advancement column: ", primary_adv)
message(">> All advancement columns: ", paste(adv_cols, collapse = ", "))

# Clean primary advancement (remove extreme outliers for some analyses)
df$adv <- df[[primary_adv]]

###############################################################################
#                                                                             #
#   SECTION 1:  AGE-STRATIFIED BIOLOGICAL AGE ADVANCEMENT                     #
#                                                                             #
###############################################################################

cat("\n")
cat("##########################################################\n")
cat("#  SECTION 1: AGE-STRATIFIED ADVANCEMENT                 #\n")
cat("##########################################################\n\n")

# Create age groups
df$age_group <- cut(df$age,
  breaks = c(0, 25, 35, 45, 55, 65, 75, 85, 120),
  labels = c("<25", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85+"),
  right = FALSE
)

df$age_decade <- cut(df$age,
  breaks = c(0, 30, 40, 50, 60, 70, 80, 120),
  labels = c("<30", "30s", "40s", "50s", "60s", "70s", "80+"),
  right = FALSE
)

# --- 1a. Summary statistics by age group (all clocks) ---
cat("--- Advancement by Age Group (all clocks) ---\n\n")

age_strat_all <- list()
for (ac in adv_cols) {
  for (ag in levels(df$age_group)) {
    idx <- df$age_group == ag & !is.na(df[[ac]]) & is.finite(df[[ac]])
    n <- sum(idx, na.rm = TRUE)
    if (n >= 5) {
      vals <- df[[ac]][idx]
      age_strat_all[[paste(ac, ag)]] <- data.frame(
        clock     = gsub("_advance$", "", ac),
        age_group = ag,
        n         = n,
        mean_adv  = round(mean(vals), 2),
        sd_adv    = round(sd(vals), 2),
        median_adv = round(median(vals), 2),
        q25       = round(quantile(vals, 0.25), 2),
        q75       = round(quantile(vals, 0.75), 2),
        pct_accelerated = round(100 * mean(vals > 0), 1),
        stringsAsFactors = FALSE
      )
    }
  }
}
age_strat_table <- do.call(rbind, age_strat_all)
rownames(age_strat_table) <- NULL
print(age_strat_table)

write.csv(age_strat_table,
          file.path(OUTPUT_DIR, "age_stratified_advancement.csv"),
          row.names = FALSE)

# --- 1b. Kruskal-Wallis test: does advancement differ across age groups? ---
cat("\n--- Kruskal-Wallis Test: Advancement ~ Age Group ---\n")
for (ac in adv_cols) {
  ok <- !is.na(df[[ac]]) & is.finite(df[[ac]]) & !is.na(df$age_group)
  if (sum(ok) >= 30) {
    kw <- kruskal.test(df[[ac]][ok] ~ df$age_group[ok])
    cat(sprintf("  %-45s chi2 = %7.2f, df = %d, p = %s\n",
                gsub("_advance$", "", ac),
                kw$statistic, kw$parameter,
                ifelse(kw$p.value < 2.2e-16, "< 2.2e-16",
                       sprintf("%.2e", kw$p.value))))
  }
}

# --- 1c. Plot: Advancement by age group (boxplot + violin) ---
for (ac in adv_cols) {
  ok <- !is.na(df[[ac]]) & is.finite(df[[ac]]) & !is.na(df$age_group)
  if (sum(ok) < 50) next
  label <- gsub("pheno_|phenoage_|_advance", "", ac)

  p <- ggplot(df[ok, ], aes(x = age_group, y = .data[[ac]])) +
    geom_violin(fill = "lightblue", alpha = 0.4, colour = NA) +
    geom_boxplot(width = 0.15, fill = "steelblue", alpha = 0.7,
                 outlier.size = 0.3, outlier.alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
                 colour = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "darkred",
               linewidth = 0.6) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("BA Advancement by Age Group: ", label),
         subtitle = sprintf("Red diamond = mean | n = %s", format(sum(ok), big.mark = ",")),
         x = "Age Group", y = "Advancement (BA - CA, years)")

  ggsave(file.path(OUTPUT_DIR, paste0("age_strat_violin_", label, ".png")),
         p, width = 10, height = 6, dpi = 150)
  print(p)
}

# --- 1d. Plot: Mean advancement trajectory across age groups ---
if (nrow(age_strat_table) > 0) {
  traj <- age_strat_table %>%
    filter(grepl("phenoage_orig|levine_no_glucose|hema_integrated", clock))

  if (nrow(traj) > 3) {
    p_traj <- ggplot(traj, aes(x = age_group, y = mean_adv,
                                group = clock, colour = clock)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_adv - sd_adv / sqrt(n),
                        ymax = mean_adv + sd_adv / sqrt(n)),
                    width = 0.2, alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "gray40") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Mean BA Advancement Trajectory Across Age Groups",
           subtitle = "Error bars = SEM",
           x = "Age Group", y = "Mean Advancement (years)",
           colour = "Clock")

    ggsave(file.path(OUTPUT_DIR, "age_trajectory_clocks.png"),
           p_traj, width = 11, height = 6, dpi = 150)
    print(p_traj)
  }
}

###############################################################################
#                                                                             #
#   SECTION 2:  SEX DIFFERENCES IN BIOLOGICAL AGING                           #
#                                                                             #
###############################################################################

cat("\n")
cat("##########################################################\n")
cat("#  SECTION 2: SEX DIFFERENCES                            #\n")
cat("##########################################################\n\n")

# --- 2a. Summary by sex (all clocks) ---
cat("--- Advancement by Sex ---\n\n")

sex_stats <- list()
for (ac in adv_cols) {
  for (s in c("F", "M")) {
    idx <- df$sex == s & !is.na(df[[ac]]) & is.finite(df[[ac]])
    n <- sum(idx, na.rm = TRUE)
    if (n >= 10) {
      vals <- df[[ac]][idx]
      sex_stats[[paste(ac, s)]] <- data.frame(
        clock = gsub("_advance$", "", ac),
        sex   = s,
        n     = n,
        mean_age = round(mean(df$age[idx], na.rm = TRUE), 1),
        mean_adv = round(mean(vals), 2),
        sd_adv   = round(sd(vals), 2),
        median_adv = round(median(vals), 2),
        pct_accelerated = round(100 * mean(vals > 0), 1),
        stringsAsFactors = FALSE
      )
    }
  }
}
sex_table <- do.call(rbind, sex_stats)
rownames(sex_table) <- NULL
print(sex_table)

write.csv(sex_table, file.path(OUTPUT_DIR, "sex_differences_advancement.csv"),
          row.names = FALSE)

# --- 2b. Wilcoxon test: sex difference in advancement ---
cat("\n--- Wilcoxon Rank-Sum Test: Male vs Female Advancement ---\n")
sex_tests <- list()
for (ac in adv_cols) {
  f_vals <- df[[ac]][df$sex == "F" & !is.na(df[[ac]]) & is.finite(df[[ac]])]
  m_vals <- df[[ac]][df$sex == "M" & !is.na(df[[ac]]) & is.finite(df[[ac]])]

  if (length(f_vals) >= 20 && length(m_vals) >= 20) {
    wt <- wilcox.test(m_vals, f_vals, conf.int = TRUE)
    # Also t-test for effect size
    tt <- t.test(m_vals, f_vals)
    cohens_d <- (mean(m_vals) - mean(f_vals)) /
      sqrt((var(m_vals) * (length(m_vals)-1) + var(f_vals) * (length(f_vals)-1)) /
           (length(m_vals) + length(f_vals) - 2))

    sex_tests[[ac]] <- data.frame(
      clock    = gsub("_advance$", "", ac),
      mean_M   = round(mean(m_vals), 2),
      mean_F   = round(mean(f_vals), 2),
      diff_MF  = round(mean(m_vals) - mean(f_vals), 2),
      cohens_d = round(cohens_d, 3),
      wilcox_p = ifelse(wt$p.value < 2.2e-16, "< 2.2e-16",
                        sprintf("%.4e", wt$p.value)),
      ttest_p  = ifelse(tt$p.value < 2.2e-16, "< 2.2e-16",
                        sprintf("%.4e", tt$p.value)),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %-40s diff(M-F) = %+.2f  Cohen's d = %.3f  p = %s\n",
                gsub("_advance$", "", ac),
                mean(m_vals) - mean(f_vals), cohens_d,
                sex_tests[[ac]]$wilcox_p))
  }
}

if (length(sex_tests) > 0) {
  sex_test_table <- do.call(rbind, sex_tests)
  rownames(sex_test_table) <- NULL
  write.csv(sex_test_table,
            file.path(OUTPUT_DIR, "sex_tests_advancement.csv"),
            row.names = FALSE)
}

# --- 2c. Sex differences by age group (interaction) ---
cat("\n--- Sex × Age Group Interaction ---\n")
sex_age_stats <- list()
for (ac in adv_cols) {
  for (ag in levels(df$age_group)) {
    for (s in c("F", "M")) {
      idx <- df$sex == s & df$age_group == ag &
             !is.na(df[[ac]]) & is.finite(df[[ac]])
      n <- sum(idx, na.rm = TRUE)
      if (n >= 5) {
        vals <- df[[ac]][idx]
        sex_age_stats[[paste(ac, ag, s)]] <- data.frame(
          clock     = gsub("_advance$", "", ac),
          age_group = ag,
          sex       = s,
          n         = n,
          mean_adv  = round(mean(vals), 2),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(sex_age_stats) > 0) {
  sex_age_table <- do.call(rbind, sex_age_stats)
  rownames(sex_age_table) <- NULL
  write.csv(sex_age_table,
            file.path(OUTPUT_DIR, "sex_age_interaction.csv"),
            row.names = FALSE)
}

# --- 2d. Plots: sex differences ---
for (ac in adv_cols) {
  ok <- !is.na(df[[ac]]) & is.finite(df[[ac]]) & df$sex %in% c("F", "M")
  if (sum(ok) < 50) next
  label <- gsub("pheno_|phenoage_|_advance", "", ac)

  # Density by sex
  p_dens <- ggplot(df[ok, ], aes(x = .data[[ac]], fill = sex)) +
    geom_density(alpha = 0.4, colour = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30") +
    scale_fill_manual(values = c("F" = "firebrick3", "M" = "steelblue")) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Advancement Distribution by Sex: ", label),
         x = "Advancement (BA - CA, years)", y = "Density", fill = "Sex")

  ggsave(file.path(OUTPUT_DIR, paste0("sex_density_", label, ".png")),
         p_dens, width = 9, height = 5, dpi = 150)
  print(p_dens)

  # Boxplot by sex × age group
  if (!is.null(df$age_group)) {
    p_int <- ggplot(df[ok, ], aes(x = age_group, y = .data[[ac]], fill = sex)) +
      geom_boxplot(alpha = 0.6, outlier.size = 0.3, position = position_dodge(0.8)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "darkred") +
      scale_fill_manual(values = c("F" = "firebrick3", "M" = "steelblue")) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0("Sex × Age Group: ", label),
           x = "Age Group", y = "Advancement (years)", fill = "Sex")

    ggsave(file.path(OUTPUT_DIR, paste0("sex_age_boxplot_", label, ".png")),
           p_int, width = 11, height = 6, dpi = 150)
    print(p_int)
  }
}

###############################################################################
#                                                                             #
#   SECTION 3:  OUTLIER ANALYSIS — FASTEST & SLOWEST AGERS                    #
#                                                                             #
###############################################################################

cat("\n")
cat("##########################################################\n")
cat("#  SECTION 3: OUTLIER ANALYSIS                           #\n")
cat("##########################################################\n\n")

# Use primary advancement for outlier classification
ok_adv <- !is.na(df$adv) & is.finite(df$adv)
df_valid <- df[ok_adv, ]

if (nrow(df_valid) >= 100) {

  # Define outlier thresholds using percentiles
  q05  <- quantile(df_valid$adv, 0.05, na.rm = TRUE)
  q10  <- quantile(df_valid$adv, 0.10, na.rm = TRUE)
  q90  <- quantile(df_valid$adv, 0.90, na.rm = TRUE)
  q95  <- quantile(df_valid$adv, 0.95, na.rm = TRUE)
  iqr_val <- IQR(df_valid$adv, na.rm = TRUE)
  q25  <- quantile(df_valid$adv, 0.25, na.rm = TRUE)
  q75  <- quantile(df_valid$adv, 0.75, na.rm = TRUE)

  # Classify
  df_valid$aging_group <- "Normal"
  df_valid$aging_group[df_valid$adv <= q10] <- "Slow Ager (≤P10)"
  df_valid$aging_group[df_valid$adv <= q05] <- "Very Slow Ager (≤P5)"
  df_valid$aging_group[df_valid$adv >= q90] <- "Fast Ager (≥P90)"
  df_valid$aging_group[df_valid$adv >= q95] <- "Very Fast Ager (≥P95)"

  df_valid$aging_group <- factor(df_valid$aging_group,
    levels = c("Very Slow Ager (≤P5)", "Slow Ager (≤P10)", "Normal",
               "Fast Ager (≥P90)", "Very Fast Ager (≥P95)"))

  # Broader grouping for comparison
  df_valid$aging_class <- "Normal (P10-P90)"
  df_valid$aging_class[df_valid$adv <= q10] <- "Slow Ager (≤P10)"
  df_valid$aging_class[df_valid$adv >= q90] <- "Fast Ager (≥P90)"
  df_valid$aging_class <- factor(df_valid$aging_class,
    levels = c("Slow Ager (≤P10)", "Normal (P10-P90)", "Fast Ager (≥P90)"))

  cat(sprintf("  Advancement thresholds:\n"))
  cat(sprintf("    P5  = %+.2f  |  P10 = %+.2f  |  P25 = %+.2f\n", q05, q10, q25))
  cat(sprintf("    P50 = %+.2f  |  P75 = %+.2f\n",
              median(df_valid$adv), q75))
  cat(sprintf("    P90 = %+.2f  |  P95 = %+.2f\n", q90, q95))

  # --- 3a. Demographics of outlier groups ---
  cat("\n--- Outlier Group Demographics ---\n\n")
  outlier_demo <- df_valid %>%
    group_by(aging_class) %>%
    summarise(
      n = n(),
      mean_age = round(mean(age, na.rm = TRUE), 1),
      sd_age   = round(sd(age, na.rm = TRUE), 1),
      pct_female = round(100 * mean(sex == "F", na.rm = TRUE), 1),
      mean_adv = round(mean(adv, na.rm = TRUE), 2),
      sd_adv   = round(sd(adv, na.rm = TRUE), 2),
      .groups = "drop"
    )
  print(as.data.frame(outlier_demo))

  write.csv(as.data.frame(outlier_demo),
            file.path(OUTPUT_DIR, "outlier_demographics.csv"),
            row.names = FALSE)

  # --- 3b. Compare extra biomarkers across outlier groups ---
  cat("\n--- Non-Clock Biomarker Comparison: Fast vs Normal vs Slow Agers ---\n")
  cat("    (Kruskal-Wallis + pairwise Wilcoxon, BH-adjusted)\n\n")

  biomarker_outlier_tests <- list()
  for (eb in extra_biomarkers) {
    ok_eb <- !is.na(df_valid[[eb]]) & is.finite(df_valid[[eb]]) & !is.na(df_valid$aging_class)
    n_ok <- sum(ok_eb)

    # Need at least 10 in each group
    group_counts <- table(df_valid$aging_class[ok_eb])
    if (all(group_counts >= 10)) {
      kw <- kruskal.test(df_valid[[eb]][ok_eb] ~ df_valid$aging_class[ok_eb])

      # Group means
      means <- tapply(df_valid[[eb]][ok_eb], df_valid$aging_class[ok_eb], mean)

      # Pairwise: Fast vs Slow
      fast_v <- df_valid[[eb]][ok_eb & df_valid$aging_class == "Fast Ager (≥P90)"]
      slow_v <- df_valid[[eb]][ok_eb & df_valid$aging_class == "Slow Ager (≤P10)"]
      pw <- wilcox.test(fast_v, slow_v)

      biomarker_outlier_tests[[eb]] <- data.frame(
        biomarker   = eb,
        n_total     = n_ok,
        mean_slow   = round(means["Slow Ager (≤P10)"], 3),
        mean_normal = round(means["Normal (P10-P90)"], 3),
        mean_fast   = round(means["Fast Ager (≥P90)"], 3),
        kw_chi2     = round(kw$statistic, 2),
        kw_p        = kw$p.value,
        fast_vs_slow_p = pw$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(biomarker_outlier_tests) > 0) {
    outlier_bm_table <- do.call(rbind, biomarker_outlier_tests)
    rownames(outlier_bm_table) <- NULL

    # Adjust p-values
    outlier_bm_table$kw_p_adj <- p.adjust(outlier_bm_table$kw_p, method = "BH")
    outlier_bm_table$fast_vs_slow_p_adj <- p.adjust(outlier_bm_table$fast_vs_slow_p, method = "BH")

    # Sort by significance
    outlier_bm_table <- outlier_bm_table[order(outlier_bm_table$kw_p_adj), ]

    # Format for display
    display_tbl <- outlier_bm_table
    display_tbl$kw_p <- sprintf("%.2e", display_tbl$kw_p)
    display_tbl$kw_p_adj <- sprintf("%.2e", display_tbl$kw_p_adj)
    display_tbl$fast_vs_slow_p <- sprintf("%.2e", display_tbl$fast_vs_slow_p)
    display_tbl$fast_vs_slow_p_adj <- sprintf("%.2e", display_tbl$fast_vs_slow_p_adj)

    cat("  Top biomarkers differing between fast and slow agers:\n\n")
    print(display_tbl, row.names = FALSE)

    write.csv(outlier_bm_table,
              file.path(OUTPUT_DIR, "outlier_biomarker_comparison.csv"),
              row.names = FALSE)
  }

  # --- 3c. Plots: Outlier analysis ---

  # Advancement distribution with outlier zones
  p_out_hist <- ggplot(df_valid, aes(x = adv)) +
    geom_histogram(bins = 80, fill = "gray70", colour = "white", alpha = 0.8) +
    geom_vline(xintercept = c(q05, q10, q90, q95),
               linetype = c("dotted", "dashed", "dashed", "dotted"),
               colour = c("blue", "blue", "red", "red"),
               linewidth = 0.7) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
    annotate("text", x = q05, y = Inf, label = "P5", colour = "blue",
             vjust = 2, hjust = 1.2, size = 3.5) +
    annotate("text", x = q10, y = Inf, label = "P10", colour = "blue",
             vjust = 2, hjust = -0.2, size = 3.5) +
    annotate("text", x = q90, y = Inf, label = "P90", colour = "red",
             vjust = 2, hjust = 1.2, size = 3.5) +
    annotate("text", x = q95, y = Inf, label = "P95", colour = "red",
             vjust = 2, hjust = -0.2, size = 3.5) +
    theme_minimal(base_size = 13) +
    labs(title = "BA Advancement Distribution with Outlier Thresholds",
         subtitle = sprintf("n = %s | primary clock: %s",
                            format(nrow(df_valid), big.mark = ","),
                            gsub("_advance$", "", primary_adv)),
         x = "Advancement (BA - CA, years)", y = "Count")

  ggsave(file.path(OUTPUT_DIR, "outlier_histogram.png"),
         p_out_hist, width = 10, height = 6, dpi = 150)
  print(p_out_hist)

  # Boxplot of significant biomarkers by aging class
  if (length(biomarker_outlier_tests) > 0) {
    sig_bm <- outlier_bm_table$biomarker[outlier_bm_table$kw_p_adj < 0.05]
    sig_bm <- head(sig_bm, 12)  # top 12

    if (length(sig_bm) > 0) {
      for (eb in sig_bm) {
        ok_eb <- !is.na(df_valid[[eb]]) & is.finite(df_valid[[eb]])
        if (sum(ok_eb) < 50) next

        p_bm <- ggplot(df_valid[ok_eb, ], aes(x = aging_class, y = .data[[eb]],
                                                fill = aging_class)) +
          geom_violin(alpha = 0.3, colour = NA) +
          geom_boxplot(width = 0.2, alpha = 0.7, outlier.size = 0.3) +
          scale_fill_manual(values = c("Slow Ager (≤P10)" = "steelblue",
                                       "Normal (P10-P90)" = "gray60",
                                       "Fast Ager (≥P90)" = "firebrick")) +
          theme_minimal(base_size = 12) +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 20, hjust = 1)) +
          labs(title = paste0("Non-Clock Biomarker: ", eb),
               subtitle = sprintf("KW adj-p = %s",
                 sprintf("%.2e", outlier_bm_table$kw_p_adj[outlier_bm_table$biomarker == eb])),
               x = "", y = eb)

        ggsave(file.path(OUTPUT_DIR, paste0("outlier_bm_", eb, ".png")),
               p_bm, width = 8, height = 5, dpi = 150)
        print(p_bm)
      }
    }
  }

} else {
  message(">> Skipping outlier analysis: insufficient valid advancement data.")
}


###############################################################################
#                                                                             #
#   SECTION 4:  NON-CLOCK BIOMARKER CORRELATIONS WITH BA ADVANCEMENT          #
#                                                                             #
###############################################################################

cat("\n")
cat("##########################################################\n")
cat("#  SECTION 4: NON-CLOCK BIOMARKER CORRELATIONS           #\n")
cat("##########################################################\n\n")

# --- 4a. Spearman correlations with primary advancement ---
cat("--- Spearman Correlations: Extra Biomarkers vs BA Advancement ---\n\n")

corr_results <- list()
for (eb in extra_biomarkers) {
  ok_both <- !is.na(df[[eb]]) & is.finite(df[[eb]]) &
             !is.na(df$adv) & is.finite(df$adv)
  n_ok <- sum(ok_both)

  if (n_ok >= 30) {
    ct <- cor.test(df$adv[ok_both], df[[eb]][ok_both], method = "spearman")
    # Also Pearson for comparison
    ct_p <- cor.test(df$adv[ok_both], df[[eb]][ok_both], method = "pearson")

    corr_results[[eb]] <- data.frame(
      biomarker   = eb,
      n           = n_ok,
      rho         = round(ct$estimate, 4),
      rho_p       = ct$p.value,
      pearson_r   = round(ct_p$estimate, 4),
      pearson_p   = ct_p$p.value,
      stringsAsFactors = FALSE
    )
  }
}

if (length(corr_results) > 0) {
  corr_table <- do.call(rbind, corr_results)
  rownames(corr_table) <- NULL

  # BH adjustment
  corr_table$rho_p_adj <- p.adjust(corr_table$rho_p, method = "BH")
  corr_table$pearson_p_adj <- p.adjust(corr_table$pearson_p, method = "BH")

  # Sort by absolute rho
  corr_table <- corr_table[order(-abs(corr_table$rho)), ]

  # Display
  display_corr <- corr_table
  display_corr$rho_p <- sprintf("%.2e", display_corr$rho_p)
  display_corr$rho_p_adj <- sprintf("%.2e", display_corr$rho_p_adj)
  display_corr$pearson_p <- sprintf("%.2e", display_corr$pearson_p)
  display_corr$pearson_p_adj <- sprintf("%.2e", display_corr$pearson_p_adj)
  print(display_corr, row.names = FALSE)

  write.csv(corr_table,
            file.path(OUTPUT_DIR, "nonclock_biomarker_correlations.csv"),
            row.names = FALSE)

  # --- 4b. Age-adjusted partial correlations ---
  cat("\n--- Age-Adjusted Correlations (residuals after regressing out age) ---\n\n")

  partial_results <- list()
  for (eb in extra_biomarkers) {
    ok_all <- !is.na(df[[eb]]) & is.finite(df[[eb]]) &
              !is.na(df$adv) & is.finite(df$adv) &
              !is.na(df$age)
    n_ok <- sum(ok_all)
    if (n_ok >= 30) {
      # Residualize both advancement and biomarker against age
      resid_adv <- residuals(lm(adv ~ age, data = df[ok_all, ]))
      resid_bm  <- residuals(lm(df[[eb]][ok_all] ~ df$age[ok_all]))

      ct <- cor.test(resid_adv, resid_bm, method = "spearman")
      partial_results[[eb]] <- data.frame(
        biomarker     = eb,
        n             = n_ok,
        partial_rho   = round(ct$estimate, 4),
        partial_p     = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(partial_results) > 0) {
    partial_table <- do.call(rbind, partial_results)
    rownames(partial_table) <- NULL
    partial_table$partial_p_adj <- p.adjust(partial_table$partial_p, method = "BH")
    partial_table <- partial_table[order(-abs(partial_table$partial_rho)), ]

    display_part <- partial_table
    display_part$partial_p <- sprintf("%.2e", display_part$partial_p)
    display_part$partial_p_adj <- sprintf("%.2e", display_part$partial_p_adj)
    print(display_part, row.names = FALSE)

    write.csv(partial_table,
              file.path(OUTPUT_DIR, "nonclock_biomarker_partial_correlations.csv"),
              row.names = FALSE)
  }

  # --- 4c. Plot: Correlation bar chart ---
  top_corr <- head(corr_table, 20)
  top_corr$direction <- ifelse(top_corr$rho > 0, "Positive", "Negative")
  top_corr$sig <- ifelse(corr_table$rho_p_adj[match(top_corr$biomarker, corr_table$biomarker)] < 0.05,
                         "Significant", "NS")

  p_bar <- ggplot(top_corr, aes(x = reorder(biomarker, rho), y = rho,
                                 fill = direction)) +
    geom_col(alpha = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c("Positive" = "firebrick", "Negative" = "steelblue")) +
    theme_minimal(base_size = 13) +
    labs(title = "Non-Clock Biomarker Correlations with BA Advancement",
         subtitle = paste0("Spearman rho | Clock: ", gsub("_advance$", "", primary_adv)),
         x = "", y = "Spearman rho", fill = "Direction")

  ggsave(file.path(OUTPUT_DIR, "nonclock_correlation_barplot.png"),
         p_bar, width = 10, height = 7, dpi = 150)
  print(p_bar)

  # --- 4d. Scatter plots for top correlated biomarkers ---
  top_scatter <- head(corr_table[abs(corr_table$rho) >= 0.05, ], 8)

  for (i in seq_len(nrow(top_scatter))) {
    eb <- top_scatter$biomarker[i]
    ok_both <- !is.na(df[[eb]]) & is.finite(df[[eb]]) &
               !is.na(df$adv) & is.finite(df$adv)
    if (sum(ok_both) < 30) next

    p_sc <- ggplot(df[ok_both, ], aes(x = .data[[eb]], y = adv)) +
      geom_point(alpha = 0.05, size = 0.5, colour = "steelblue") +
      geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
      geom_smooth(method = "loess", colour = "firebrick", se = FALSE,
                  linewidth = 0.8, linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "gray40") +
      theme_minimal(base_size = 12) +
      labs(title = paste0(eb, " vs BA Advancement"),
           subtitle = sprintf("rho = %.3f, n = %s",
                              top_scatter$rho[i],
                              format(sum(ok_both), big.mark = ",")),
           x = eb, y = "Advancement (years)")

    ggsave(file.path(OUTPUT_DIR, paste0("scatter_adv_", eb, ".png")),
           p_sc, width = 8, height = 6, dpi = 150)
    print(p_sc)
  }
}

###############################################################################
#                                                                             #
#   SECTION 5:  MULTI-BIOMARKER CORRELATION MATRIX & TENDENCIES               #
#                                                                             #
###############################################################################

cat("\n")
cat("##########################################################\n")
cat("#  SECTION 5: MULTI-BIOMARKER CORRELATION MATRIX         #\n")
cat("##########################################################\n\n")

# Build a matrix of all extra biomarkers + all advancements
matrix_vars <- c(adv_cols, extra_biomarkers)
matrix_vars <- intersect(matrix_vars, names(df))

# Need at least 50 non-NA for a variable
keep_vars <- sapply(matrix_vars, function(v) {
  sum(!is.na(df[[v]]) & is.finite(df[[v]])) >= 50
})
matrix_vars <- matrix_vars[keep_vars]

if (length(matrix_vars) >= 5) {

  # Compute pairwise Spearman correlation matrix
  mat_data <- df[, matrix_vars, drop = FALSE]
  for (v in matrix_vars) mat_data[[v]] <- as.numeric(mat_data[[v]])

  cor_mat <- cor(mat_data, use = "pairwise.complete.obs", method = "spearman")

  # Save
  write.csv(cor_mat, file.path(OUTPUT_DIR, "full_correlation_matrix.csv"))

  # Plot correlation heatmap
  png(file.path(OUTPUT_DIR, "correlation_heatmap.png"),
      width = 1200, height = 1200, res = 100)
  if (USE_CORRPLOT) {
    corrplot(cor_mat, method = "color", type = "lower",
             tl.cex = 0.7, tl.col = "black",
             col = colorRampPalette(c("steelblue", "white", "firebrick"))(200),
             addCoef.col = "black", number.cex = 0.5,
             title = "Spearman Correlation Matrix: Advancements + Extra Biomarkers",
             mar = c(0, 0, 2, 0))
  } else {
    heatmap(cor_mat, scale = "none", symm = TRUE,
            col = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
            main = "Spearman Correlation Matrix")
  }
  dev.off()
  message(">> Saved: correlation_heatmap.png")

  # --- 5b. Focus heatmap: extra biomarkers vs advancement columns only ---
  adv_in_mat <- intersect(adv_cols, matrix_vars)
  extra_in_mat <- intersect(extra_biomarkers, matrix_vars)

  if (length(adv_in_mat) >= 1 && length(extra_in_mat) >= 3) {
    cross_cor <- cor_mat[extra_in_mat, adv_in_mat, drop = FALSE]

    png(file.path(OUTPUT_DIR, "cross_correlation_heatmap.png"),
        width = 900, height = 1000, res = 100)
    if (USE_CORRPLOT) {
      corrplot(cross_cor, method = "color", is.corr = FALSE,
               tl.cex = 0.8, tl.col = "black",
               col = colorRampPalette(c("steelblue", "white", "firebrick"))(200),
               addCoef.col = "black", number.cex = 0.6,
               cl.lim = c(min(cross_cor, na.rm=TRUE), max(cross_cor, na.rm=TRUE)),
               title = "Non-Clock Biomarkers vs Clock Advancements",
               mar = c(0, 0, 2, 0))
    } else {
      heatmap(cross_cor, scale = "none", Colv = NA,
              col = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
              main = "Non-Clock Biomarkers vs Clock Advancements")
    }
    dev.off()
    message(">> Saved: cross_correlation_heatmap.png")
  }

  # --- 5c. Biomarker tendency profiles by aging class ---
  if (exists("df_valid") && "aging_class" %in% names(df_valid)) {
    cat("\n--- Standardized Biomarker Profiles by Aging Class ---\n\n")

    profile_data <- list()
    for (eb in extra_in_mat) {
      ok_eb <- !is.na(df_valid[[eb]]) & is.finite(df_valid[[eb]])
      if (sum(ok_eb) < 50) next

      # Z-score relative to full cohort
      mu <- mean(df_valid[[eb]][ok_eb], na.rm = TRUE)
      s  <- sd(df_valid[[eb]][ok_eb], na.rm = TRUE)
      if (s == 0) next

      for (ac in levels(df_valid$aging_class)) {
        idx <- ok_eb & df_valid$aging_class == ac
        if (sum(idx) < 10) next
        profile_data[[paste(eb, ac)]] <- data.frame(
          biomarker = eb,
          aging_class = ac,
          mean_z = round((mean(df_valid[[eb]][idx]) - mu) / s, 3),
          n = sum(idx),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(profile_data) > 0) {
      profile_table <- do.call(rbind, profile_data)
      rownames(profile_table) <- NULL

      write.csv(profile_table,
                file.path(OUTPUT_DIR, "biomarker_zprofiles_by_aging_class.csv"),
                row.names = FALSE)

      # Plot: radar-style parallel coordinates
      profile_wide <- pivot_wider(profile_table,
                                  names_from = aging_class,
                                  values_from = mean_z)

      # Compute absolute max z for ordering
      profile_wide$max_diff <- abs(
        profile_wide[["Fast Ager (≥P90)"]] -
        profile_wide[["Slow Ager (≤P10)"]]
      )
      profile_wide <- profile_wide[order(-profile_wide$max_diff), ]
      top_profile <- head(profile_wide, 15)

      # Long format for plotting
      top_long <- pivot_longer(
        top_profile[, c("biomarker", "Slow Ager (≤P10)", "Normal (P10-P90)", "Fast Ager (≥P90)")],
        cols = -biomarker, names_to = "aging_class", values_to = "mean_z"
      )
      top_long$aging_class <- factor(top_long$aging_class,
        levels = c("Slow Ager (≤P10)", "Normal (P10-P90)", "Fast Ager (≥P90)"))
      top_long$biomarker <- factor(top_long$biomarker,
        levels = rev(top_profile$biomarker))

      p_prof <- ggplot(top_long, aes(x = mean_z, y = biomarker,
                                      colour = aging_class)) +
        geom_vline(xintercept = 0, colour = "gray80") +
        geom_point(size = 3) +
        geom_line(aes(group = biomarker), colour = "gray70", linewidth = 0.3) +
        scale_colour_manual(values = c("Slow Ager (≤P10)" = "steelblue",
                                       "Normal (P10-P90)" = "gray50",
                                       "Fast Ager (≥P90)" = "firebrick")) +
        theme_minimal(base_size = 12) +
        labs(title = "Standardized Biomarker Profiles: Fast vs Slow Agers",
             subtitle = "Mean Z-score relative to full cohort | Top 15 by Fast-Slow difference",
             x = "Mean Z-score", y = "", colour = "Aging Class")

      ggsave(file.path(OUTPUT_DIR, "biomarker_profiles_by_aging_class.png"),
             p_prof, width = 11, height = 7, dpi = 150)
      print(p_prof)
    }
  }
}

###############################################################################
#                                                                             #
#   SECTION 6:  SEX-STRATIFIED BIOMARKER CORRELATIONS                         #
#                                                                             #
###############################################################################

cat("\n")
cat("##########################################################\n")
cat("#  SECTION 6: SEX-STRATIFIED BIOMARKER CORRELATIONS      #\n")
cat("##########################################################\n\n")

sex_corr_results <- list()
for (s in c("F", "M")) {
  df_sex <- df[df$sex == s, ]
  for (eb in extra_biomarkers) {
    ok_both <- !is.na(df_sex[[eb]]) & is.finite(df_sex[[eb]]) &
               !is.na(df_sex$adv) & is.finite(df_sex$adv)
    n_ok <- sum(ok_both)
    if (n_ok >= 30) {
      ct <- cor.test(df_sex$adv[ok_both], df_sex[[eb]][ok_both], method = "spearman")
      sex_corr_results[[paste(eb, s)]] <- data.frame(
        biomarker = eb, sex = s, n = n_ok,
        rho = round(ct$estimate, 4),
        p = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(sex_corr_results) > 0) {
  sex_corr_table <- do.call(rbind, sex_corr_results)
  rownames(sex_corr_table) <- NULL
  sex_corr_table$p_adj <- p.adjust(sex_corr_table$p, method = "BH")
  sex_corr_table <- sex_corr_table[order(-abs(sex_corr_table$rho)), ]

  # Display
  display_sex <- sex_corr_table
  display_sex$p <- sprintf("%.2e", display_sex$p)
  display_sex$p_adj <- sprintf("%.2e", display_sex$p_adj)
  print(head(display_sex, 30), row.names = FALSE)

  write.csv(sex_corr_table,
            file.path(OUTPUT_DIR, "sex_stratified_biomarker_correlations.csv"),
            row.names = FALSE)

  # Plot: Compare rho between sexes
  wide_sex <- pivot_wider(sex_corr_table[, c("biomarker", "sex", "rho")],
                          names_from = sex, values_from = rho,
                          names_prefix = "rho_")

  if ("rho_F" %in% names(wide_sex) && "rho_M" %in% names(wide_sex)) {
    p_sex_rho <- ggplot(wide_sex, aes(x = rho_F, y = rho_M, label = biomarker)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray60") +
      geom_hline(yintercept = 0, colour = "gray80") +
      geom_vline(xintercept = 0, colour = "gray80") +
      geom_point(size = 2.5, colour = "purple4", alpha = 0.7) +
      geom_text(size = 2.8, hjust = -0.15, vjust = -0.3, colour = "gray30") +
      theme_minimal(base_size = 13) +
      labs(title = "Sex-Specific Correlations of Non-Clock Biomarkers with BA Advancement",
           subtitle = "Points off diagonal = sex-differential biomarkers",
           x = "Spearman rho (Female)", y = "Spearman rho (Male)")

    ggsave(file.path(OUTPUT_DIR, "sex_rho_comparison.png"),
           p_sex_rho, width = 9, height = 8, dpi = 150)
    print(p_sex_rho)
  }
}


###############################################################################
#   FINAL SUMMARY & SAVE                                                      #
###############################################################################

cat("\n")
cat("==========================================================\n")
cat("  DEEP ANALYSIS COMPLETE\n")
cat("==========================================================\n")
cat("  Subjects analysed:    ", format(nrow(df), big.mark = ","), "\n")
cat("  Primary clock:        ", gsub("_advance$", "", primary_adv), "\n")
cat("  Extra biomarkers:     ", length(extra_biomarkers), "\n")
cat("  Age groups:           ", nlevels(df$age_group), "\n")
cat("  Output directory:     ", OUTPUT_DIR, "\n")
cat("\n  Files saved:\n")

saved_files <- list.files(OUTPUT_DIR, full.names = FALSE)
for (f in saved_files) cat("    - ", f, "\n")

cat("==========================================================\n")

message("\n>> Done. All outputs in: ", OUTPUT_DIR)
