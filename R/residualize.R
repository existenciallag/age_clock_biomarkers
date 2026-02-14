###############################################################################
# residualize.R — Residualization and Z-score standardization
#
# Residual biological age = the component of BA that is independent of
# chronological age.  This is obtained by regressing BA (or BA advancement)
# on age and taking the residual.
#
# Residualization is preferred over raw advancement (BA − CA) because the
# latter can retain a residual correlation with chronological age, biasing
# Cox models.
#
# After residualization, Z-scoring (mean = 0, SD = 1) makes hazard ratios
# interpretable as "per 1 SD of biological deviation".
#
# Reference: Levine ME (2013); Cohen AA et al. (2014)
###############################################################################

# ── Single-variable residualization ────────────────────────────────────────

#' Residualize a variable against chronological age
#'
#' Uses lm(y ~ age) with na.exclude so the output vector keeps the same
#' length as the input (NAs preserved in place).
#'
#' @param y   Numeric vector (BA or BA advancement).
#' @param age Numeric vector of chronological ages.
#' @return Numeric vector of residuals (same length as y).
residualize <- function(y, age) {
  resid(lm(y ~ age, na.action = na.exclude))
}

# ── Z-score standardization ───────────────────────────────────────────────

#' Standardize to mean 0, SD 1
#' @param x Numeric vector.
#' @return Numeric vector of Z-scores.
zscore <- function(x) as.numeric(scale(x))

# ── Batch residualization + Z-scoring ─────────────────────────────────────

#' Residualize and Z-score a set of clock advancement columns
#'
#' For each column in `clock_vars`, two new columns are added to `data`:
#'   - <varname>_resid   : raw residual from lm(var ~ age)
#'   - <varname>_resid_z : Z-scored residual
#'
#' @param data       Assembled data.frame (must contain an `age` column).
#' @param clock_vars Character vector of advancement column names to process.
#'                   Only columns that exist in `data` are processed.
#' @return The input data.frame with new _resid and _resid_z columns appended.
add_residuals <- function(data, clock_vars) {

  clock_vars <- intersect(clock_vars, names(data))

  for (v in clock_vars) {
    resid_col   <- paste0(v, "_resid")
    resid_z_col <- paste0(v, "_resid_z")

    data[[resid_col]]   <- residualize(data[[v]], data$age)
    data[[resid_z_col]] <- zscore(data[[resid_col]])
  }

  message(">> Residualized ", length(clock_vars), " clocks")
  data
}

# ── Also Z-score the raw BA columns (for raw Cox models) ──────────────────

#' Z-score raw BA / advancement columns
#'
#' Adds z_<varname> columns.
#'
#' @param data       data.frame.
#' @param clock_vars Column names to Z-score.
#' @return data.frame with z_* columns appended.
add_zscores <- function(data, clock_vars) {

  clock_vars <- intersect(clock_vars, names(data))

  for (v in clock_vars) {
    data[[paste0("z_", v)]] <- zscore(data[[v]])
  }

  data
}

# ── Convenience: list residual columns ────────────────────────────────────

list_resid_columns <- function(data) {
  grep("_resid$", names(data), value = TRUE)
}

list_resid_z_columns <- function(data) {
  grep("_resid_z$", names(data), value = TRUE)
}
