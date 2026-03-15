###############################################################################
# qc_recenter.R — QC & statistical comparison: Argentine cohort vs NHANES
#
# Provides functions used by run_recenter.R to:
#   1. Compute per-biomarker distributional statistics (mean, SD, median, IQR)
#   2. Compare Argentine vs NHANES biomarker distributions (Welch t, KS, SMD)
#   3. Compute clock-level QC metrics (N, coverage, BA-age correlation)
#   4. Compare advancement distributions before/after recentering
#   5. Generate diagnostic plots
#
# All functions are pure (no side effects) unless they return ggplot objects.
###############################################################################

# ── 1. Biomarker distribution summary ────────────────────────────────────────

#' Summarise distribution of a set of biomarker columns
#'
#' @param data   data.frame
#' @param vars   character vector of column names
#' @param label  character tag added as "cohort" column
#' @return data.frame with one row per variable
biomarker_summary <- function(data, vars, label = "cohort") {
  vars <- intersect(vars, names(data))
  out <- lapply(vars, function(v) {
    x <- data[[v]]
    x <- x[!is.na(x) & is.finite(x)]
    if (length(x) == 0) return(NULL)
    data.frame(
      cohort   = label,
      variable = v,
      n        = length(x),
      mean     = mean(x),
      sd       = sd(x),
      median   = median(x),
      q25      = quantile(x, 0.25),
      q75      = quantile(x, 0.75),
      min      = min(x),
      max      = max(x),
      skew     = if (length(x) > 2) {
        m3 <- mean((x - mean(x))^3); m3 / (sd(x)^3)
      } else NA_real_,
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  do.call(rbind, out)
}

# ── 2. Two-cohort biomarker comparison ───────────────────────────────────────

#' Compare biomarker distributions between two cohorts
#'
#' For each variable present in both data frames, computes:
#'   - Standardised Mean Difference (SMD = (mean_B - mean_A) / pooled_SD)
#'   - Welch t-test p-value
#'   - Two-sample Kolmogorov-Smirnov p-value
#'
#' @param data_a  data.frame (reference, e.g. NHANES IV)
#' @param data_b  data.frame (target, e.g. Argentine cohort)
#' @param vars    character vector of column names to compare
#' @param label_a label for reference cohort
#' @param label_b label for target cohort
#' @return data.frame with comparison statistics
compare_biomarkers <- function(data_a, data_b, vars,
                               label_a = "NHANES", label_b = "Argentina") {
  vars <- intersect(vars, intersect(names(data_a), names(data_b)))
  out <- lapply(vars, function(v) {
    xa <- data_a[[v]]; xa <- xa[!is.na(xa) & is.finite(xa)]
    xb <- data_b[[v]]; xb <- xb[!is.na(xb) & is.finite(xb)]
    if (length(xa) < 10 || length(xb) < 10) return(NULL)

    pooled_sd <- sqrt((var(xa) * (length(xa) - 1) +
                       var(xb) * (length(xb) - 1)) /
                      (length(xa) + length(xb) - 2))

    smd <- if (pooled_sd > 0) (mean(xb) - mean(xa)) / pooled_sd else NA_real_

    tt <- tryCatch(t.test(xa, xb), error = function(e) NULL)
    ks <- tryCatch(ks.test(xa, xb), error = function(e) NULL)

    data.frame(
      variable    = v,
      n_ref       = length(xa),
      mean_ref    = mean(xa),
      sd_ref      = sd(xa),
      n_target    = length(xb),
      mean_target = mean(xb),
      sd_target   = sd(xb),
      SMD         = round(smd, 4),
      welch_p     = if (!is.null(tt)) tt$p.value else NA_real_,
      ks_p        = if (!is.null(ks)) ks$p.value else NA_real_,
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  result <- do.call(rbind, out)
  if (!is.null(result)) {
    result$smd_flag <- ifelse(abs(result$SMD) > 0.5, "LARGE",
                       ifelse(abs(result$SMD) > 0.2, "moderate", "small"))
  }
  result
}

# ── 3. Clock-level QC for scored cohort ──────────────────────────────────────

#' QC summary table for a scored clinical dataset
#'
#' @param data       Scored data.frame (from clinical_scored_lite.csv)
#' @param adv_cols   Character vector of advancement columns
#' @return data.frame with QC metrics per clock
clinical_qc_table <- function(data, adv_cols = NULL) {
  if (is.null(adv_cols))
    adv_cols <- grep("_advance$", names(data), value = TRUE)

  out <- lapply(adv_cols, function(ac) {
    ba_col <- sub("_advance$", "", ac)
    x_adv <- data[[ac]]
    x_ba  <- if (ba_col %in% names(data)) data[[ba_col]] else NULL
    ok <- !is.na(x_adv) & is.finite(x_adv)

    r_ba_age <- if (!is.null(x_ba)) {
      ok_ba <- !is.na(x_ba) & is.finite(x_ba)
      if (sum(ok_ba) > 30) cor(data$age[ok_ba], x_ba[ok_ba]) else NA_real_
    } else NA_real_

    data.frame(
      clock      = ba_col,
      n_scored   = sum(ok),
      n_total    = nrow(data),
      coverage   = round(mean(ok), 4),
      mean_adv   = if (sum(ok) > 0) mean(x_adv[ok]) else NA_real_,
      sd_adv     = if (sum(ok) > 1) sd(x_adv[ok])   else NA_real_,
      median_adv = if (sum(ok) > 0) median(x_adv[ok]) else NA_real_,
      r_ba_age   = round(r_ba_age, 4),
      pct_accel  = if (sum(ok) > 0) round(100 * mean(x_adv[ok] > 0), 1) else NA_real_,
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  do.call(rbind, out)
}

# ── 4. Compare advancement before/after recentering ─────────────────────────

#' Summary statistics for advancement pre- vs post-recentering
#'
#' @param data        data.frame with both original and recentered columns
#' @param orig_cols   character vector of original advancement columns
#' @param local_cols  character vector of recentered advancement columns
#' @return data.frame
compare_recentering <- function(data, orig_cols, local_cols) {
  stopifnot(length(orig_cols) == length(local_cols))

  out <- lapply(seq_along(orig_cols), function(i) {
    oc <- orig_cols[i]; lc <- local_cols[i]
    xo <- data[[oc]]; xl <- data[[lc]]
    ok <- !is.na(xo) & is.finite(xo) & !is.na(xl) & is.finite(xl)
    n  <- sum(ok)
    if (n < 10) return(NULL)

    clock <- sub("_advance$", "", oc)
    data.frame(
      clock            = clock,
      n                = n,
      # Original (NHANES-referenced)
      orig_mean        = round(mean(xo[ok]), 3),
      orig_sd          = round(sd(xo[ok]), 3),
      orig_median      = round(median(xo[ok]), 3),
      orig_pct_accel   = round(100 * mean(xo[ok] > 0), 1),
      # Local recentered
      local_mean       = round(mean(xl[ok]), 3),
      local_sd         = round(sd(xl[ok]), 3),
      local_median     = round(median(xl[ok]), 3),
      local_pct_accel  = round(100 * mean(xl[ok] > 0), 1),
      # Correlation (should be ~1.0 for linear recentering)
      cor_orig_local   = round(cor(xo[ok], xl[ok]), 4),
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  do.call(rbind, out)
}

# ── 5. Diagnostic plots ────────────────────────────────────────────────────

#' Density overlay: original vs recentered advancement
#'
#' @param data      data.frame
#' @param orig_col  original advancement column name
#' @param local_col recentered advancement column name
#' @param title     plot title
#' @return ggplot object
plot_recenter_density <- function(data, orig_col, local_col, title = NULL) {
  ok <- !is.na(data[[orig_col]]) & is.finite(data[[orig_col]]) &
        !is.na(data[[local_col]]) & is.finite(data[[local_col]])
  pdata <- data.frame(
    Original   = data[[orig_col]][ok],
    Recentered = data[[local_col]][ok]
  )
  pdata_long <- tidyr::pivot_longer(pdata, cols = everything(),
                                    names_to = "method", values_to = "advancement")
  clock_label <- sub("_advance$", "", orig_col)
  if (is.null(title)) title <- paste0("Advancement Distribution: ", clock_label)

  ggplot(pdata_long, aes(x = advancement, fill = method, colour = method)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "gray30") +
    scale_fill_manual(values = c("Original" = "steelblue", "Recentered" = "firebrick")) +
    scale_colour_manual(values = c("Original" = "steelblue", "Recentered" = "firebrick")) +
    theme_minimal(base_size = 13) +
    labs(title = title,
         subtitle = sprintf("n = %s", format(sum(ok), big.mark = ",")),
         x = "Advancement (years)", y = "Density",
         fill = "Method", colour = "Method")
}

#' QQ-plot of advancement against normal distribution
#'
#' @param x     numeric vector of advancement values
#' @param title plot title
#' @return ggplot object
plot_advancement_qq <- function(x, title = "QQ Plot: Advancement") {
  x <- x[!is.na(x) & is.finite(x)]
  df <- data.frame(sample = sort(x))
  n <- length(x)
  df$theoretical <- qnorm(ppoints(n))

  ggplot(df, aes(x = theoretical, y = sample)) +
    geom_point(alpha = 0.3, size = 0.8, colour = "steelblue") +
    geom_abline(slope = sd(x), intercept = mean(x),
                colour = "red", linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = title, x = "Theoretical Quantiles", y = "Sample Quantiles")
}

#' SMD bar chart for biomarker comparison
#'
#' @param comp_table  output of compare_biomarkers()
#' @param title       plot title
#' @return ggplot object
plot_smd_chart <- function(comp_table, title = "Standardised Mean Differences: Argentina vs NHANES") {
  ct <- comp_table[order(abs(comp_table$SMD), decreasing = TRUE), ]
  ct$variable <- factor(ct$variable, levels = ct$variable)

  ggplot(ct, aes(x = SMD, y = variable, fill = smd_flag)) +
    geom_col(alpha = 0.8) +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", colour = "gray50") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", colour = "red") +
    scale_fill_manual(values = c("small" = "steelblue",
                                 "moderate" = "orange",
                                 "LARGE" = "firebrick")) +
    theme_minimal(base_size = 12) +
    labs(title = title,
         subtitle = "Dashed = |0.2| threshold; Dotted = |0.5| threshold",
         x = "SMD (Argentina − NHANES) / pooled SD",
         y = NULL, fill = "Magnitude")
}
