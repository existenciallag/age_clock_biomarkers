###############################################################################
# stratify.R — Risk stratification, Kaplan-Meier curves, categorical Cox
#
# Assigns subjects to risk groups (tertiles by default) based on their
# residual biological age, then evaluates group-level survival.
#
# Clinical interpretation:
#   Tertile 1 = low biological risk  (younger than expected)
#   Tertile 2 = intermediate risk
#   Tertile 3 = high biological risk (older than expected)
###############################################################################

# ── Assign risk tertiles ──────────────────────────────────────────────────

#' Assign risk groups based on a residual age column
#'
#' @param data       data.frame.
#' @param resid_col  Name of the residualized column (e.g. "phenoage_orig_advance_resid").
#' @param n_groups   Number of quantile groups (default 3 = tertiles).
#' @param labels     Labels for the groups.
#' @param col_name   Name for the new factor column.
#' @return data.frame with the new risk-group column added.
add_risk_groups <- function(data,
                            resid_col,
                            n_groups = 3,
                            labels   = c("Low biological risk",
                                         "Intermediate risk",
                                         "High biological risk"),
                            col_name = "risk_group") {

  stopifnot(resid_col %in% names(data))

  data[[col_name]] <- dplyr::ntile(data[[resid_col]], n_groups)
  data[[col_name]] <- factor(
    data[[col_name]],
    levels = seq_len(n_groups),
    labels = labels[seq_len(n_groups)]
  )

  data
}

# ── Kaplan-Meier by risk group ────────────────────────────────────────────

#' Fit and plot Kaplan-Meier curves by risk group
#'
#' @param data     data.frame with time, status, and risk_group columns.
#' @param group_col Name of the risk-group factor column.
#' @param colours  Colour vector for the KM curves.
#' @param title    Plot title.
#' @return The survfit object (invisibly).
plot_km <- function(data,
                    group_col = "risk_group",
                    colours   = KM_COLOURS,
                    title     = "Survival by Residual Biological Age") {

  fml <- as.formula(paste0("survival::Surv(time, status) ~ ", group_col))
  km  <- survival::survfit(fml, data = data)

  plot(
    km,
    col  = colours,
    lwd  = 2,
    xlab = "Follow-up time (months)",
    ylab = "Survival probability",
    main = title
  )

  legend(
    "bottomleft",
    legend = levels(data[[group_col]]),
    col    = colours,
    lwd    = 2,
    bty    = "n"
  )

  invisible(km)
}

# ── Categorical Cox model ─────────────────────────────────────────────────

#' Cox model with risk group as a categorical predictor
#'
#' Tests whether higher risk groups have significantly worse survival
#' after adjusting for chronological age.
#'
#' @param data      data.frame.
#' @param group_col Name of the risk-group column.
#' @return A coxph summary object.
cox_categorical <- function(data, group_col = "risk_group") {

  fml <- as.formula(paste0("survival::Surv(time, status) ~ age + ", group_col))
  fit <- survival::coxph(fml, data = data)
  summary(fit)
}

# ── Multi-clock risk stratification ───────────────────────────────────────

#' Run risk stratification for multiple clocks at once
#'
#' For each residualized clock, adds risk groups and fits a categorical
#' Cox model. Returns a summary data.frame.
#'
#' @param data        data.frame with residualized columns.
#' @param resid_cols  Character vector of _resid column names.
#' @param n_groups    Number of quantile groups.
#' @return A data.frame: clock, group, N, HR, LCL, UCL, p.
stratify_multi <- function(data, resid_cols, n_groups = 3) {

  resid_cols <- intersect(resid_cols, names(data))

  results <- lapply(resid_cols, function(rc) {
    grp_col <- paste0(rc, "_grp")
    tmp <- add_risk_groups(data, rc, n_groups = n_groups, col_name = grp_col)

    # Need time/status/age
    cols_needed <- c("time", "status", "age", grp_col)
    if (!all(cols_needed %in% names(tmp))) return(NULL)

    sub <- tmp[complete.cases(tmp[, cols_needed]), ]
    if (nrow(sub) < COX_MIN_N) return(NULL)

    fml <- as.formula(paste0("survival::Surv(time, status) ~ age + ", grp_col))
    fit <- tryCatch(survival::coxph(fml, data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    s <- summary(fit)
    # Extract coefficients for the risk-group levels (skip age)
    idx <- grep(grp_col, rownames(s$coefficients))
    if (length(idx) == 0) return(NULL)

    data.frame(
      clock = rc,
      group = rownames(s$coefficients)[idx],
      N     = nrow(sub),
      HR    = s$coefficients[idx, "exp(coef)"],
      LCL   = s$conf.int[idx, "lower .95"],
      UCL   = s$conf.int[idx, "upper .95"],
      p     = s$coefficients[idx, "Pr(>|z|)"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  do.call(rbind, results)
}

# ── Density plot of residual age distributions ────────────────────────────

#' Faceted density plots for residual biological age
#'
#' @param data       data.frame.
#' @param resid_cols Character vector of residual columns.
#' @return A ggplot object.
plot_resid_density <- function(data, resid_cols) {
  resid_cols <- intersect(resid_cols, names(data))

  long <- tidyr::pivot_longer(
    data[, resid_cols, drop = FALSE],
    cols      = dplyr::everything(),
    names_to  = "clock",
    values_to = "residual_age"
  )

  # Clean labels
  long$clock <- gsub("pheno_|phenoage_|_advance|_resid_z|_resid|_orig", "",
                     long$clock)

  ggplot2::ggplot(long, ggplot2::aes(x = residual_age)) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.4) +
    ggplot2::facet_wrap(~ clock, scales = "free") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Distribution of Residual Biological Age",
      x     = "Residual age (years)",
      y     = "Density"
    )
}
