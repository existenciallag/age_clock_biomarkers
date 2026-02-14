###############################################################################
# clinical_score.R — Score new patients using a trained deployment bundle
#
# This module applies pre-trained BioAge models to new laboratory data
# WITHOUT retraining.  It is designed for clinical validation on local
# lab populations where mortality outcomes are unavailable.
#
# Workflow:
#   1. Load the deployment bundle (from export.R)
#   2. Read new patient data (CSV / data.frame)
#   3. Map local biomarker names to NHANES names
#   4. Apply trained PhenoAge / KDM / HD models
#   5. Compute advancement, residual age, and Z-scores
#   6. Produce a per-patient clinical report table
#
# The key insight: BioAge's phenoage_calc() / kdm_calc() functions can
# score new data using pre-trained fit objects.  If those are not exported,
# we replicate the PhenoAge scoring formula from Levine (2018).
###############################################################################

# ── Biomarker name mapping ─────────────────────────────────────────────────

#' Default mapping from common clinical lab names → NHANES variable names
#'
#' Users should modify this to match their local LIMS export format.
#' Left = your lab name, Right = NHANES name expected by BioAge.
DEFAULT_NAME_MAP <- c(
  # Clinical name           = NHANES name
  "albumin_g_L"             = "albumin_gL",
  "albumin_gdL"             = "albumin",
  "alkaline_phosphatase"    = "alp",
  "log_crp"                 = "lncrp",
  "crp_log"                 = "lncrp",
  "total_cholesterol"       = "totchol",
  "creatinine_umol_log"     = "lncreat_umol",
  "creatinine_log"          = "lncreat",
  "hba1c_pct"               = "hba1c",
  "systolic_bp"             = "sbp",
  "blood_urea_nitrogen"     = "bun",
  "uric_acid"               = "uap",
  "lymphocyte_pct"          = "lymph",
  "mean_corpuscular_volume" = "mcv",
  "white_blood_cells"       = "wbc",
  "red_cell_distribution"   = "rdw",
  "red_blood_cells"         = "rbc",
  "gamma_gt"                = "ggt",
  "triglycerides"           = "trig",
  "vitamin_b12"             = "vitaminB12",
  "fasting_insulin"         = "insulin"
)

# ── Map column names ──────────────────────────────────────────────────────

#' Rename columns in a patient data.frame to NHANES conventions
#'
#' @param patient_data  data.frame of new patient lab results.
#' @param name_map      Named character vector: clinical name → NHANES name.
#'                      Defaults to DEFAULT_NAME_MAP.
#' @return data.frame with renamed columns (unmatched columns left as-is).
map_biomarker_names <- function(patient_data, name_map = DEFAULT_NAME_MAP) {
  current <- names(patient_data)
  for (i in seq_along(name_map)) {
    idx <- match(names(name_map)[i], current)
    if (!is.na(idx)) {
      current[idx] <- name_map[i]
    }
  }
  names(patient_data) <- current
  patient_data
}

# ── PhenoAge scoring from first principles ─────────────────────────────────
#
# PhenoAge is computed in two steps (Levine 2018):
#   Step 1: Gompertz mortality score
#     xb = b0 + sum(bi * xi)  where xi = biomarkers
#     M  = 1 - exp( -exp(xb) * (exp(120 * gamma) - 1) / gamma )
#   Step 2: Invert Gompertz to get PhenoAge
#     PhenoAge = 141.50225 + ln(-0.00553 * ln(1 - M)) / 0.090165
#
# The coefficients (b0, bi, gamma) come from the trained fit object.

#' Score PhenoAge for new data using trained coefficients
#'
#' @param new_data   data.frame with biomarker columns.
#' @param fit        Fit object from the deployment bundle.
#' @param biomarkers Character vector of biomarker names for this clock.
#' @return Numeric vector of PhenoAge values (one per row).
score_phenoage <- function(new_data, fit, biomarkers) {

  # If BioAge provides a scoring function, use it
  if (exists("phenoage_calc", where = asNamespace("BioAge"), inherits = FALSE)) {
    return(BioAge::phenoage_calc(new_data, fit))
  }

  # Manual scoring from Levine (2018) formula
  # The fit object should contain $coefficients and Gompertz parameters
  # Attempt to extract from the fit structure
  if (!is.null(fit$coefficients)) {
    coefs <- fit$coefficients
    # Construct linear predictor
    bm <- intersect(biomarkers, names(new_data))
    if (length(bm) == 0) {
      warning("No matching biomarkers found in new data")
      return(rep(NA_real_, nrow(new_data)))
    }

    mat <- as.matrix(new_data[, bm, drop = FALSE])

    # Match coefficient names to biomarker columns
    coef_names <- names(coefs)
    matched <- intersect(coef_names, bm)

    if (length(matched) < length(bm)) {
      warning("Some biomarkers missing from fit coefficients: ",
              paste(setdiff(bm, matched), collapse = ", "))
    }

    xb <- as.numeric(mat[, matched, drop = FALSE] %*% coefs[matched])

    # Add intercept if present
    if ("(Intercept)" %in% coef_names) {
      xb <- xb + coefs["(Intercept)"]
    }

    # Gompertz parameters from Levine (2018)
    gamma <- 0.090165
    intercept_phenoage <- 141.50225
    slope_phenoage     <- -0.00553

    # Mortality probability
    M <- 1 - exp(-exp(xb) * (exp(120 * gamma) - 1) / gamma)

    # Invert to PhenoAge
    phenoage <- intercept_phenoage + log(slope_phenoage * log(1 - M)) / gamma

    return(phenoage)
  }

  warning("Could not score PhenoAge: fit object structure not recognized")
  rep(NA_real_, nrow(new_data))
}

# ── Score all clocks for new patients ─────────────────────────────────────

#' Apply all trained clocks to a new patient dataset
#'
#' @param patient_data  data.frame with columns: patient_id, age, + biomarkers.
#' @param bundle        Deployment bundle from load_bundle().
#' @param name_map      Optional biomarker name mapping.
#' @return A data.frame with one row per patient, containing:
#'   - BA values for each successfully scored clock
#'   - Advancement (BA − CA) for each clock
#'   - Residualized age (using saved regression coefficients)
#'   - Z-scores (using saved population statistics)
score_patients <- function(patient_data, bundle, name_map = NULL) {

  if (!is.null(name_map)) {
    patient_data <- map_biomarker_names(patient_data, name_map)
  }

  # Ensure required columns
  if (!"age" %in% names(patient_data)) {
    stop("patient_data must contain an 'age' column")
  }
  if (!"patient_id" %in% names(patient_data)) {
    patient_data$patient_id <- seq_len(nrow(patient_data))
  }

  result <- data.frame(
    patient_id = patient_data$patient_id,
    age        = patient_data$age,
    stringsAsFactors = FALSE
  )

  # ---- Score canonical clocks ---------------------------------------------

  # PhenoAge original
  tryCatch({
    ba <- score_phenoage(patient_data, bundle$models$pheno_orig,
                         bundle$panels$pheno_orig)
    result$phenoage_orig <- ba
    result$phenoage_orig_advance <- ba - patient_data$age
  }, error = function(e) message("  ! PhenoAge orig scoring failed: ", e$message))

  # PhenoAge global
  tryCatch({
    ba <- score_phenoage(patient_data, bundle$models$pheno_global,
                         bundle$panels$pheno_global)
    result$phenoage_global <- ba
    result$phenoage_global_advance <- ba - patient_data$age
  }, error = function(e) message("  ! PhenoAge global scoring failed: ", e$message))

  # ---- Score sub-clocks ---------------------------------------------------
  for (nm in names(bundle$subclock_models)) {
    col_ba  <- paste0("pheno_", nm)
    col_adv <- paste0("pheno_", nm, "_advance")

    tryCatch({
      ba <- score_phenoage(patient_data, bundle$subclock_models[[nm]],
                           bundle$panels$subclocks[[nm]])
      result[[col_ba]]  <- ba
      result[[col_adv]] <- ba - patient_data$age
    }, error = function(e) {
      message("  ! Sub-clock '", nm, "' scoring failed: ", e$message)
    })
  }

  # ---- Residualize using saved coefficients ─────────────────────────────
  if (!is.null(bundle$resid_coefs)) {
    for (i in seq_len(nrow(bundle$resid_coefs))) {
      v   <- bundle$resid_coefs$variable[i]
      b0  <- bundle$resid_coefs$intercept[i]
      b1  <- bundle$resid_coefs$slope[i]

      if (v %in% names(result)) {
        predicted <- b0 + b1 * result$age
        result[[paste0(v, "_resid")]] <- result[[v]] - predicted
      }
    }
  }

  # ---- Z-score using saved population stats ──────────────────────────────
  if (!is.null(bundle$pop_stats)) {
    for (i in seq_len(nrow(bundle$pop_stats))) {
      v    <- bundle$pop_stats$variable[i]
      mu   <- bundle$pop_stats$mean[i]
      s    <- bundle$pop_stats$sd[i]

      resid_col <- paste0(v, "_resid")
      if (resid_col %in% names(result) && s > 0) {
        result[[paste0(v, "_resid_z")]] <- (result[[resid_col]] - mu) / s
      }
    }
  }

  message(">> Scored ", nrow(result), " patients across ",
          sum(grepl("_advance$", names(result))), " clocks")

  result
}

# ── Clinical report table ─────────────────────────────────────────────────

#' Generate a per-patient clinical summary table
#'
#' Produces a human-readable table suitable for inclusion in a lab report.
#'
#' @param scored_data  Output of score_patients().
#' @param bundle       Deployment bundle (for HR context).
#' @return A data.frame with patient_id, age, and per-clock columns:
#'   BA, advancement, residual, z-score, risk_level.
clinical_report_table <- function(scored_data, bundle = NULL) {

  adv_cols   <- grep("_advance$", names(scored_data), value = TRUE)
  resid_cols <- grep("_resid$",   names(scored_data), value = TRUE)

  # Build a long-format report
  rows <- list()
  for (ac in adv_cols) {
    clock_name <- sub("_advance$", "", ac)
    ba_col     <- clock_name
    resid_col  <- paste0(ac, "_resid")
    z_col      <- paste0(ac, "_resid_z")

    for (i in seq_len(nrow(scored_data))) {
      ba   <- if (ba_col %in% names(scored_data)) scored_data[[ba_col]][i] else NA
      adv  <- scored_data[[ac]][i]
      res  <- if (resid_col %in% names(scored_data)) scored_data[[resid_col]][i] else NA
      zsc  <- if (z_col %in% names(scored_data)) scored_data[[z_col]][i] else NA

      # Risk level based on z-score
      risk <- if (is.na(zsc)) {
        NA_character_
      } else if (zsc <= -1) {
        "Low"
      } else if (zsc <= 1) {
        "Average"
      } else {
        "Elevated"
      }

      # System label
      short <- sub("^pheno_|^phenoage_", "", clock_name)
      system <- if (short %in% names(PANEL_SYSTEM)) {
        PANEL_SYSTEM[short]
      } else {
        "Global"
      }

      rows[[length(rows) + 1]] <- data.frame(
        patient_id  = scored_data$patient_id[i],
        age         = scored_data$age[i],
        system      = system,
        clock       = clock_name,
        bio_age     = round(ba, 1),
        advancement = round(adv, 1),
        residual    = round(res, 2),
        z_score     = round(zsc, 2),
        risk_level  = risk,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}
