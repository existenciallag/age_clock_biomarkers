###############################################################################
# export.R — Export trained clock objects for clinical deployment
#
# Saves everything needed to score new patients WITHOUT retraining:
#   - Trained model objects (fit parameters from BioAge)
#   - Population reference statistics (means, SDs for Z-scoring)
#   - Panel definitions
#   - HR results for context
#
# The exported .rds file is self-contained: load it on any machine with
# base R + BioAge to score new lab data.
###############################################################################

# ── Extract population reference statistics ───────────────────────────────

#' Compute means and SDs from the NHANES IV training population
#'
#' These are needed to Z-score new patients against the same reference.
#'
#' @param data        Assembled NHANES IV data.frame.
#' @param clock_vars  Advancement column names to extract stats for.
#' @return A data.frame: variable, mean, sd, n.
extract_population_stats <- function(data, clock_vars) {
  clock_vars <- intersect(clock_vars, names(data))

  out <- lapply(clock_vars, function(v) {
    x <- data[[v]]
    data.frame(
      variable = v,
      mean     = mean(x, na.rm = TRUE),
      sd       = sd(x, na.rm = TRUE),
      n        = sum(!is.na(x)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Extract the age regression coefficients used for residualization
#'
#' For each advancement clock, fit lm(advance ~ age) on NHANES IV and
#' store intercept + slope so new patients can be residualized without
#' access to the training data.
#'
#' @param data       Assembled data.frame.
#' @param clock_vars Advancement column names.
#' @return A data.frame: variable, intercept, slope.
extract_resid_coefficients <- function(data, clock_vars) {
  clock_vars <- intersect(clock_vars, names(data))

  out <- lapply(clock_vars, function(v) {
    df <- data[, c(v, "age"), drop = FALSE]
    df <- df[complete.cases(df), ]
    if (nrow(df) < 100) return(NULL)

    fit <- lm(as.formula(paste(v, "~ age")), data = df)
    data.frame(
      variable  = v,
      intercept = coef(fit)[1],
      slope     = coef(fit)[2],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  do.call(rbind, out)
}

# ── Build the deployment bundle ───────────────────────────────────────────

#' Create a self-contained deployment object
#'
#' @param clocks       Output of train_all_clocks().
#' @param data         Assembled + processed data.frame.
#' @param clock_vars   Advancement column names used in analysis.
#' @param hr_results   HR table (from hr_table(), optional).
#' @param qc           QC table (from qc_table(), optional).
#' @return A list with all components needed for clinical scoring.
build_deployment_bundle <- function(clocks, data, clock_vars,
                                    hr_results = NULL, qc = NULL) {

  bundle <- list(
    # Model objects — contain the trained coefficients
    models = list(
      pheno_orig   = clocks$pheno_orig$fit,
      pheno_global = clocks$pheno_global$fit,
      kdm          = if (!is.null(clocks$kdm)) clocks$kdm$fit else NULL,
      hd           = if (!is.null(clocks$hd))  clocks$hd$fit  else NULL
    ),

    # Sub-clock model objects
    subclock_models = lapply(clocks$subclocks, function(x) x$fit),

    # Panel definitions (so we know which biomarkers each clock needs)
    panels = list(
      pheno_orig   = BIOMARKERS_PHENO_ORIG,
      pheno_global = PANEL_GLOBAL,
      kdm_hd       = BIOMARKERS_KDM_HD,
      subclocks    = PANELS_SUB[names(clocks$subclocks)]
    ),

    # System groupings
    panel_system = PANEL_SYSTEM,

    # Population reference stats for Z-scoring
    pop_stats = extract_population_stats(data, clock_vars),

    # Age regression coefficients for residualization
    resid_coefs = extract_resid_coefficients(data, clock_vars),

    # HR results for clinical context
    hr_results = hr_results,

    # QC table
    qc = qc,

    # Metadata
    meta = list(
      created    = Sys.time(),
      R_version  = R.version.string,
      n_subjects = nrow(data),
      clocks_trained = c(
        "phenoage_orig", "phenoage_global",
        if (!is.null(clocks$kdm)) "kdm",
        if (!is.null(clocks$hd))  "hd",
        paste0("pheno_", names(clocks$subclocks))
      )
    )
  )

  bundle
}

# ── Save / load ───────────────────────────────────────────────────────────

#' Save deployment bundle to disk
#'
#' @param bundle  Output of build_deployment_bundle().
#' @param path    File path (.rds).
save_bundle <- function(bundle, path = "bioage_deployment_bundle.rds") {
  saveRDS(bundle, path)
  message(">> Deployment bundle saved to: ", path)
  message("   Contains ", length(bundle$meta$clocks_trained), " trained clocks")
  message("   Trained on ", bundle$meta$n_subjects, " NHANES IV subjects")
}

#' Load a deployment bundle
#'
#' @param path File path (.rds).
#' @return The deployment bundle list.
load_bundle <- function(path = "bioage_deployment_bundle.rds") {
  bundle <- readRDS(path)
  message(">> Bundle loaded: ", length(bundle$meta$clocks_trained), " clocks, ",
          "created ", bundle$meta$created)
  bundle
}
