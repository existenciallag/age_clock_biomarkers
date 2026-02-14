###############################################################################
# survival.R — Cox proportional-hazards models, HR tables, forest plots
#
# Two families of models:
#   1. Raw models    — Surv(time, status) ~ age + z_<clock>
#                      (BA Z-scored but NOT residualized)
#   2. Residual models — Surv(time, status) ~ age + <clock>_resid_z
#                      (residualized then Z-scored)
#
# All HRs are "per +1 SD" to enable direct comparison across clocks.
#
# Requires: survival package
###############################################################################

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required for survival.R")
}

# ── Fit a single Cox model ─────────────────────────────────────────────────

#' Fit a Cox model for one predictor, adjusted for age (+/- gender)
#'
#' @param data        data.frame with time, status, age, and the predictor.
#' @param predictor   Name of the predictor column (should be Z-scored).
#' @param adjust_gender Logical; include gender as covariate?
#' @param min_n       Minimum sample size (rows after na.omit).
#' @return A one-row data.frame with HR, CI, p, C-index, AIC — or NULL.
fit_cox_single <- function(data, predictor, adjust_gender = FALSE,
                           min_n = COX_MIN_N) {

  covars <- c("age", predictor)
  if (adjust_gender && "gender" %in% names(data)) covars <- c(covars, "gender")

  cols <- c("time", "status", covars)
  cols <- intersect(cols, names(data))

  df <- data[, cols, drop = FALSE]
  df <- df[complete.cases(df), ]

  if (nrow(df) < min_n) return(NULL)

  fml <- as.formula(paste0("survival::Surv(time, status) ~ ",
                           paste(covars, collapse = " + ")))
  fit <- survival::coxph(fml, data = df)
  s   <- summary(fit)


  # The predictor is the second coefficient (first = age)
  idx <- which(rownames(s$coefficients) == predictor)
  if (length(idx) == 0) idx <- 2

  data.frame(
    clock  = predictor,
    N      = nrow(df),
    HR     = s$coefficients[idx, "exp(coef)"],
    LCL    = s$conf.int[idx, "lower .95"],
    UCL    = s$conf.int[idx, "upper .95"],
    p      = s$coefficients[idx, "Pr(>|z|)"],
    Cindex = s$concordance[1],
    AIC    = AIC(fit),
    stringsAsFactors = FALSE
  )
}

# ── Batch: raw BA Cox models ──────────────────────────────────────────────

#' Run Cox models on Z-scored BA columns (not residualized)
#'
#' @param data       Assembled + Z-scored data.frame.
#' @param clock_vars BA column names (without the z_ prefix).
#' @param adjust_gender Logical.
#' @return data.frame of HR results, one row per clock.
cox_raw <- function(data, clock_vars, adjust_gender = FALSE) {

  z_vars <- paste0("z_", clock_vars)
  z_vars <- intersect(z_vars, names(data))

  results <- lapply(z_vars, function(v) {
    fit_cox_single(data, v, adjust_gender = adjust_gender)
  })

  out <- do.call(rbind, results)
  if (!is.null(out)) {
    out$model <- "Raw (+1 SD)"
    out$clock <- sub("^z_", "", out$clock)
    out <- out[order(out$HR, decreasing = TRUE), ]
  }
  out
}

# ── Batch: residual age Cox models ────────────────────────────────────────

#' Run Cox models on residualized + Z-scored columns
#'
#' @param data       Assembled + residualized data.frame.
#' @param clock_vars Advancement column names (without _resid_z suffix).
#' @param adjust_gender Logical.
#' @return data.frame of HR results.
cox_residual <- function(data, clock_vars, adjust_gender = FALSE) {

  rz_vars <- paste0(clock_vars, "_resid_z")
  rz_vars <- intersect(rz_vars, names(data))

  results <- lapply(rz_vars, function(v) {
    fit_cox_single(data, v, adjust_gender = adjust_gender)
  })

  out <- do.call(rbind, results)
  if (!is.null(out)) {
    out$model <- "Residual (+1 SD)"
    out$clock <- sub("_resid_z$", "", out$clock)
    out <- out[order(out$HR, decreasing = TRUE), ]
  }
  out
}

# ── Combined HR table ─────────────────────────────────────────────────────

#' Build a combined HR table (raw + residual) for all clocks
#'
#' @param data       Fully processed data.frame.
#' @param clock_vars Advancement column names to evaluate.
#' @param adjust_gender Logical.
#' @return data.frame with columns: clock, model, N, HR, LCL, UCL, p, Cindex, AIC
hr_table <- function(data, clock_vars, adjust_gender = FALSE) {
  raw <- cox_raw(data, clock_vars, adjust_gender)
  res <- cox_residual(data, clock_vars, adjust_gender)
  rbind(raw, res)
}

# ── Forest plot ────────────────────────────────────────────────────────────

#' Forest plot of hazard ratios
#'
#' @param hr_df     data.frame from hr_table() or cox_residual().
#' @param title     Plot title.
#' @param subtitle  Plot subtitle.
#' @param colour    Point/bar colour.
#' @return A ggplot object.
forest_plot <- function(hr_df,
                        title    = "Hazard Ratios by Biological Clock",
                        subtitle = "Per +1 SD (adjusted for age)",
                        colour   = FOREST_COLOUR) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for forest_plot()")
  }

  # Clean up clock labels
  hr_df$label <- gsub("pheno_|phenoage_|_advance|_orig", "", hr_df$clock)

  p <- ggplot2::ggplot(
    hr_df,
    ggplot2::aes(
      x    = stats::reorder(label, HR),
      y    = HR,
      ymin = LCL,
      ymax = UCL
    )
  ) +
    ggplot2::geom_pointrange(size = 0.8, colour = colour) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x = "",
      y = "Hazard Ratio (95% CI)"
    )

  # Facet by model type if both raw and residual present
  if ("model" %in% names(hr_df) && length(unique(hr_df$model)) > 1) {
    p <- p + ggplot2::facet_wrap(~ model)
  }

  p
}
