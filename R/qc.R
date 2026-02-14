###############################################################################
# qc.R — Quality control diagnostics for assembled clock data
#
# Produces a QC summary table and determines which clocks have enough
# valid data to be used in survival analysis.
###############################################################################

# ── QC diagnostic table ────────────────────────────────────────────────────

#' Build a QC summary for all advancement / HD columns
#'
#' @param data      Assembled data.frame from assemble_data().
#' @param clock_vars Character vector of column names to audit.
#'                   Defaults to all advancement + HD columns.
#' @return A data.frame with one row per clock:
#'   variable, n_total, n_nonNA, n_finite, sd, frac_valid
qc_table <- function(data, clock_vars = NULL) {

  if (is.null(clock_vars)) {
    clock_vars <- c(
      list_advance_columns(data),
      intersect(c("hd", "hd_log"), names(data))
    )
  }

  clock_vars <- intersect(clock_vars, names(data))

  out <- lapply(clock_vars, function(v) {
    x <- data[[v]]
    data.frame(
      variable   = v,
      n_total    = nrow(data),
      n_nonNA    = sum(!is.na(x)),
      n_finite   = sum(is.finite(x), na.rm = TRUE),
      sd         = sd(x, na.rm = TRUE),
      frac_valid = mean(is.finite(x), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out)
  result[order(result$frac_valid), ]
}

# ── Determine valid clocks ─────────────────────────────────────────────────

#' Filter QC table to clocks meeting minimum thresholds
#'
#' @param qc     Output of qc_table().
#' @param min_n  Minimum number of finite observations (default QC_MIN_N).
#' @param min_sd Minimum SD (default QC_MIN_SD).
#' @return Character vector of valid clock variable names.
valid_clocks <- function(qc, min_n = QC_MIN_N, min_sd = QC_MIN_SD) {
  qc$variable[qc$n_finite > min_n & qc$sd > min_sd]
}

# ── QC correlation check: PhenoAge orig vs age ─────────────────────────────

#' Print the Pearson correlation between PhenoAge original and chronological
#' age as a sanity check (expected ~0.94 on NHANES IV).
#'
#' @param data Assembled data.frame.
#' @return The correlation value (invisibly).
qc_pheno_age_cor <- function(data) {
  r <- cor(data$phenoage_orig, data$age, use = "complete.obs")
  message(">> QC: cor(PhenoAge_orig, age) = ", round(r, 4),
          "  [expected ~0.94]")
  invisible(r)
}

# ── QC overlap analysis ───────────────────────────────────────────────────

#' Report how many subjects have complete data across a set of clocks
#'
#' @param data       Assembled data.frame.
#' @param clock_vars Character vector of columns to check.
#' @return A data.frame: n_complete, n_total, frac_complete.
qc_overlap <- function(data, clock_vars = NULL) {
  if (is.null(clock_vars)) clock_vars <- list_ba_columns(data)
  clock_vars <- intersect(clock_vars, names(data))

  complete <- complete.cases(data[, clock_vars, drop = FALSE])
  data.frame(
    n_complete    = sum(complete),
    n_total       = nrow(data),
    frac_complete = mean(complete),
    clocks_tested = length(clock_vars)
  )
}
