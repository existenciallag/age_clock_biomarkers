###############################################################################
# train.R — Train biological age clocks using the BioAge package
#
# Wraps BioAge's *_nhanes() functions to train canonical clocks (PhenoAge,
# KDM, HD) and an arbitrary set of PhenoAge sub-clocks.
#
# All training happens on NHANES III; projection is onto NHANES IV.
# The BioAge package handles both steps internally — we never touch
# or modify its embedded reference datasets.
#
# Returns a named list of trained clock objects, each containing:
#   $data  — NHANES IV projected data.frame
#   $fit   — Model fit object (coefficients, parameters)
###############################################################################

# NOTE: This module requires R/config.R to be sourced first.
# Use run_pipeline.R to orchestrate all modules in the correct order.

# ── Core training function ─────────────────────────────────────────────────

#' Train all biological age clocks
#'
#' @param panels_sub  Named list of sub-clock biomarker vectors.
#'                    Defaults to PANELS_SUB from config.R.
#' @param train_kdm   Logical; train KDM clock?
#' @param train_hd    Logical; train HD clock?
#' @param verbose     Logical; print progress messages.
#'
#' @return A named list with components:
#'   \item{pheno_orig}{PhenoAge original (12 biomarkers)}
#'   \item{pheno_global}{PhenoAge global panel}
#'   \item{kdm}{KDM biological age (if train_kdm = TRUE)}
#'   \item{hd}{Homeostatic Dysregulation (if train_hd = TRUE)}
#'   \item{subclocks}{Named list of successfully trained sub-clocks}
#'   \item{failed}{Character vector of sub-clock names that failed}
train_all_clocks <- function(panels_sub = PANELS_SUB,
                             train_kdm  = TRUE,
                             train_hd   = TRUE,
                             verbose    = TRUE) {

  if (!requireNamespace("BioAge", quietly = TRUE)) {
    stop("BioAge package is required. Install with:\n",
         "  devtools::install_github('dayoonkwon/BioAge')")
  }

  result <- list()


  # ---- Canonical clocks ---------------------------------------------------
  if (verbose) message(">> Training PhenoAge original (12 biomarkers)")
  result$pheno_orig <- BioAge::phenoage_nhanes(BIOMARKERS_PHENO_ORIG)

  if (verbose) message(">> Training PhenoAge global panel")
  result$pheno_global <- BioAge::phenoage_nhanes(PANEL_GLOBAL)

  if (train_kdm) {
    if (verbose) message(">> Training KDM biological age")
    result$kdm <- BioAge::kdm_nhanes(BIOMARKERS_KDM_HD)
  }

  if (train_hd) {
    if (verbose) message(">> Training Homeostatic Dysregulation (HD)")
    result$hd <- BioAge::hd_nhanes(BIOMARKERS_KDM_HD)
  }

  # ---- Sub-clocks ---------------------------------------------------------
  if (verbose) message(">> Training ", length(panels_sub), " sub-clocks")

  trained   <- list()
  failed    <- character(0)

  for (nm in names(panels_sub)) {
    bm <- panels_sub[[nm]]
    if (verbose) message("   Training sub-clock: ", nm,
                         " [", paste(bm, collapse=", "), "]")
    clock <- tryCatch(
      phenoage_nhanes_extended(bm),
      error = function(e) {
        if (verbose) message("   ! Sub-clock '", nm, "' failed: ", e$message)
        NULL
      }
    )
    if (!is.null(clock)) {
      trained[[nm]] <- clock
    } else {
      failed <- c(failed, nm)
    }
  }

  result$subclocks <- trained
  result$failed    <- failed

  if (verbose) {
    message(">> Training complete: ",
            length(trained), " sub-clocks OK, ",
            length(failed), " failed")
  }

  result
}

# ── Custom phenoage_nhanes with derived variables ─────────────────────────

#' Extended phenoage_nhanes that adds derived columns to NHANES data
#'
#' The standard BioAge::phenoage_nhanes() only creates albumin (g/L),
#' glucose (mmol/L), creat (μmol/L), lncreat (log μmol/L).
#' This wrapper adds additional derived variables (e.g., lnglucose)
#' so sub-clocks can use log-transformed biomarkers.
#'
#' @param biomarkers Character vector of biomarker names.
#' @return Same structure as BioAge::phenoage_nhanes().
phenoage_nhanes_extended <- function(biomarkers) {

  # Derived variables that need custom NHANES columns
  needs_custom <- any(biomarkers %in% c("lnglucose"))

  if (!needs_custom) {
    # Standard path — no custom columns needed
    return(BioAge::phenoage_nhanes(biomarkers))
  }

  # Custom path: replicate phenoage_nhanes logic with extra mutate steps
  train <- BioAge::phenoage_calc(
    data = BioAge::NHANES3 %>%
      dplyr::filter(age >= 20 & age <= 84) %>%
      dplyr::mutate(
        albumin   = albumin_gL,
        glucose   = glucose_mmol,
        creat     = creat_umol,
        lncreat   = lncreat_umol,
        lnglucose = ifelse(glucose_mmol > 0, log(glucose_mmol), NA)
      ),
    biomarkers = biomarkers,
    fit = NULL
  )

  test <- BioAge::phenoage_calc(
    data = BioAge::NHANES4 %>%
      dplyr::filter(age >= 20) %>%
      dplyr::mutate(
        albumin   = albumin_gL,
        glucose   = glucose_mmol,
        creat     = creat_umol,
        lncreat   = lncreat_umol,
        lnglucose = ifelse(glucose_mmol > 0, log(glucose_mmol), NA)
      ),
    biomarkers = biomarkers,
    fit = train$fit
  )

  dat <- dplyr::left_join(
    BioAge::NHANES4,
    test$data[, c("sampleID", "phenoage", "phenoage_advance")],
    by = "sampleID"
  )

  result <- list(data = dat, fit = train$fit)
  class(result) <- append(class(result), "phenoage")
  result
}

# ── Convenience: train a single custom sub-clock ───────────────────────────

#' Train a single PhenoAge sub-clock
#'
#' @param biomarkers Character vector of biomarker names.
#' @param name       Optional label for messages.
#' @return BioAge phenoage_nhanes result, or NULL on failure.
train_subclock <- function(biomarkers, name = "custom") {
  tryCatch(
    phenoage_nhanes_extended(biomarkers),
    error = function(e) {
      message("Sub-clock '", name, "' failed: ", e$message)
      NULL
    }
  )
}
