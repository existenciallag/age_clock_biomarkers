###############################################################################
# clinical_score.R — Score new patients using trained models or published
#                    Levine 2018 PhenoAge coefficients
#
# Two scoring modes:
#   A) Published PhenoAge: uses hardcoded Levine 2018 coefficients (9 biomarkers)
#      → No bundle needed, exact replication of the original paper
#   B) Bundle-based: uses trained coefficients from BioAge pipeline
#      → For sub-clocks and custom panels
#
# Both modes produce: biological age, advancement (BA - CA).
# Bundle mode additionally provides residualization and Z-scoring.
###############################################################################

# ── Biomarker name mapping ─────────────────────────────────────────────────

DEFAULT_NAME_MAP <- c(
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

map_biomarker_names <- function(patient_data, name_map = DEFAULT_NAME_MAP) {
  current <- names(patient_data)
  for (i in seq_along(name_map)) {
    idx <- match(names(name_map)[i], current)
    if (!is.na(idx)) current[idx] <- name_map[i]
  }
  names(patient_data) <- current
  patient_data
}

# ══════════════════════════════════════════════════════════════════════════
# A) PUBLISHED LEVINE 2018 PHENOAGE — hardcoded coefficients
# ══════════════════════════════════════════════════════════════════════════

#' Score PhenoAge using exact published Levine (2018) coefficients
#'
#' This does NOT require a trained bundle. Uses the 9 biomarkers + age
#' from the original paper with the published Cox/Gompertz parameters.
#'
#' @param data  data.frame with columns: albumin, alp, glucose, lncrp,
#'              lncreat, lymph, mcv, rdw, wbc, age
#' @return Numeric vector of PhenoAge values (NA where any biomarker missing)
score_phenoage_levine2018 <- function(data) {

  required <- c("albumin", "alp", "glucose", "lncrp", "lncreat",
                 "lymph", "mcv", "rdw", "wbc", "age")
  available <- intersect(required, names(data))
  missing   <- setdiff(required, available)

  if (length(missing) > 0) {
    message("  Levine PhenoAge: missing columns: ", paste(missing, collapse = ", "))
    if (length(missing) > 2) {
      message("  Too many missing — returning NA")
      return(rep(NA_real_, nrow(data)))
    }
  }

  n <- nrow(data)
  coefs <- LEVINE_2018_COEFS

  # Build linear predictor
  xb <- rep(coefs["intercept"], n)

  bm_names <- c("albumin", "alp", "glucose", "lncrp", "lncreat",
                 "lymph", "mcv", "rdw", "wbc", "age")
  for (bm in bm_names) {
    if (bm %in% names(data)) {
      xb <- xb + coefs[bm] * as.numeric(data[[bm]])
    } else {
      xb <- rep(NA_real_, n)  # Can't compute without all variables
      break
    }
  }

  # Gompertz mortality probability (120 months)
  gamma <- LEVINE_GAMMA_MORT
  M <- 1 - exp(-exp(xb) * (exp(120 * gamma) - 1) / gamma)

  # Invert to PhenoAge
  phenoage <- LEVINE_PHENOAGE_A +
    log(LEVINE_PHENOAGE_B * log(1 - M)) / LEVINE_PHENOAGE_C

  as.numeric(phenoage)
}

# ══════════════════════════════════════════════════════════════════════════
# B) BUNDLE-BASED SCORING — uses trained fit objects
# ══════════════════════════════════════════════════════════════════════════

#' Extract coefficient vector from a BioAge fit object
#'
#' BioAge's phenoage_nhanes() uses flexsurv::flexsurvreg() internally.
#' The fit object may be a flexsurvreg, coxph, or custom list.
#' This function tries every possible extraction method.
extract_fit_coefs <- function(fit) {

  # Method 1: try coef() generically (works for coxph, flexsurvreg, glm, etc.)
  coefs <- tryCatch(coef(fit), error = function(e) NULL)
  if (!is.null(coefs) && is.numeric(coefs) && length(coefs) > 1) {
    return(coefs)
  }

  # Method 2: flexsurvreg stores coefficients in $res or $res.t
  if (!is.null(fit$res) && is.matrix(fit$res)) {
    coefs <- fit$res[, "est"]
    if (!is.null(coefs) && length(coefs) > 1) return(coefs)
  }
  if (!is.null(fit$res.t) && is.matrix(fit$res.t)) {
    coefs <- fit$res.t[, "est"]
    if (!is.null(coefs) && length(coefs) > 1) return(coefs)
  }

  # Method 3: direct $coefficients (numeric vector)
  if (!is.null(fit$coefficients)) {
    if (is.numeric(fit$coefficients) && length(fit$coefficients) > 1) {
      return(fit$coefficients)
    }
    # Could be a matrix: extract first column
    if (is.matrix(fit$coefficients)) {
      coefs <- fit$coefficients[, 1]
      names(coefs) <- rownames(fit$coefficients)
      return(coefs)
    }
    # Could be a data.frame
    if (is.data.frame(fit$coefficients)) {
      coefs <- fit$coefficients[, 1]
      names(coefs) <- rownames(fit$coefficients)
      return(coefs)
    }
  }

  # Method 4: $coef
  if (!is.null(fit$coef) && is.numeric(fit$coef) && length(fit$coef) > 1) {
    return(fit$coef)
  }

  # Method 5: recursive search for named numeric vectors > length 2
  if (is.list(fit)) {
    for (nm in names(fit)) {
      v <- fit[[nm]]
      if (is.numeric(v) && length(v) > 2 && !is.null(names(v))) {
        return(v)
      }
      if (is.list(v)) {
        for (nm2 in names(v)) {
          v2 <- v[[nm2]]
          if (is.numeric(v2) && length(v2) > 2 && !is.null(names(v2))) {
            return(v2)
          }
        }
      }
    }
  }

  # DEBUG: print the structure so user can report it
  message("  DEBUG: fit object structure:")
  message("    class: ", paste(class(fit), collapse = ", "))
  message("    names: ", paste(names(fit), collapse = ", "))
  for (nm in names(fit)) {
    v <- fit[[nm]]
    message("    $", nm, ": class=", paste(class(v), collapse = ","),
            " length=", length(v),
            if (is.numeric(v) && length(v) <= 20) paste0(" values=", paste(round(v, 4), collapse = ",")) else "")
  }

  NULL
}

#' Extract Gompertz gamma from a fit object
extract_fit_gamma <- function(fit) {
  if (!is.null(fit$gamma) && is.numeric(fit$gamma)) return(fit$gamma)
  if (is.list(fit)) {
    for (nm in names(fit)) {
      if (nm == "gamma" && is.numeric(fit[[nm]])) return(fit[[nm]])
    }
  }
  # Default: Levine 2018 Gompertz gamma
  LEVINE_GAMMA_MORT
}

#' Score PhenoAge using a trained fit object from the bundle
#'
#' @param new_data    data.frame with biomarker columns
#' @param fit         Fit object from the deployment bundle
#' @param biomarkers  Character vector of biomarker names for this clock
#' @return Numeric vector of PhenoAge values
score_phenoage <- function(new_data, fit, biomarkers) {

  n <- nrow(new_data)

  # Extract coefficients
  coefs <- extract_fit_coefs(fit)
  if (is.null(coefs)) {
    warning("Could not extract coefficients from fit object (class: ",
            paste(class(fit), collapse = "/"), ")")
    return(rep(NA_real_, n))
  }

  # Match biomarkers: must exist in BOTH data and coefficients
  bm_in_data  <- intersect(biomarkers, names(new_data))
  bm_in_coefs <- intersect(biomarkers, names(coefs))
  bm_matched  <- intersect(bm_in_data, bm_in_coefs)

  if (length(bm_matched) == 0) {
    warning("No biomarker match. Data has: ",
            paste(bm_in_data, collapse = ", "),
            ". Coefs have: ", paste(names(coefs), collapse = ", "))
    return(rep(NA_real_, n))
  }

  if (length(bm_matched) < length(biomarkers)) {
    missing <- setdiff(biomarkers, bm_matched)
    message("  Note: ", length(missing), " biomarkers unavailable: ",
            paste(missing, collapse = ", "))
  }

  # Linear predictor: xb = X %*% beta
  mat <- as.matrix(new_data[, bm_matched, drop = FALSE])
  xb  <- as.numeric(mat %*% coefs[bm_matched])

  # Add age component if it's in the coefficients but NOT in biomarkers list
  if ("age" %in% names(coefs) && !"age" %in% biomarkers) {
    if ("age" %in% names(new_data)) {
      xb <- xb + coefs["age"] * as.numeric(new_data$age)
    }
  }

  # Gompertz gamma from fit
  gamma <- extract_fit_gamma(fit)

  # Mortality probability (120 months)
  M <- 1 - exp(-exp(xb) * (exp(120 * gamma) - 1) / gamma)

  # PhenoAge inversion
  phenoage <- LEVINE_PHENOAGE_A +
    log(LEVINE_PHENOAGE_B * log(1 - M)) / LEVINE_PHENOAGE_C

  as.numeric(phenoage)
}

# ══════════════════════════════════════════════════════════════════════════
# MAIN: Score all clocks for new patients
# ══════════════════════════════════════════════════════════════════════════

#' Apply all trained clocks to a new patient dataset
#'
#' @param patient_data  data.frame with columns: patient_id, age, + biomarkers
#' @param bundle        Deployment bundle from load_bundle()
#' @param name_map      Optional biomarker name mapping
#' @return data.frame with BA values, advancement, residuals, Z-scores
score_patients <- function(patient_data, bundle, name_map = NULL) {

  if (!is.null(name_map)) {
    patient_data <- map_biomarker_names(patient_data, name_map)
  }

  if (!"age" %in% names(patient_data)) stop("patient_data must have 'age'")
  if (!"patient_id" %in% names(patient_data)) {
    patient_data$patient_id <- seq_len(nrow(patient_data))
  }

  result <- data.frame(
    patient_id = patient_data$patient_id,
    age        = patient_data$age,
    stringsAsFactors = FALSE
  )

  # ---- 1. Published Levine PhenoAge (9 biomarkers + age) ------------------
  tryCatch({
    ba <- score_phenoage_levine2018(patient_data)
    n_ok <- sum(!is.na(ba))
    if (n_ok > 0) {
      result$phenoage_levine <- ba
      result$phenoage_levine_advance <- ba - patient_data$age
      message(">> PhenoAge Levine 2018: ", n_ok, " / ", nrow(result), " scored")
    } else {
      message(">> PhenoAge Levine 2018: 0 scored (missing biomarkers)")
    }
  }, error = function(e) {
    message("  ! Levine PhenoAge failed: ", conditionMessage(e))
  })

  # ---- 2. Bundle-based canonical clocks -----------------------------------
  for (clock_name in c("pheno_orig", "pheno_global")) {
    col_ba  <- paste0("phenoage_", sub("pheno_", "", clock_name))
    col_adv <- paste0(col_ba, "_advance")

    model <- bundle$models[[clock_name]]
    panel <- bundle$panels[[clock_name]]

    if (is.null(model) || is.null(panel)) next

    tryCatch({
      ba <- score_phenoage(patient_data, model, panel)
      n_ok <- sum(!is.na(ba))
      if (n_ok > 0) {
        result[[col_ba]]  <- ba
        result[[col_adv]] <- ba - patient_data$age
        message(">> ", col_ba, ": ", n_ok, " scored")
      } else {
        message(">> ", col_ba, ": 0 scored")
      }
    }, error = function(e) {
      message("  ! ", col_ba, " failed: ", conditionMessage(e))
    })
  }

  # ---- 3. Sub-clocks ------------------------------------------------------
  for (nm in names(bundle$subclock_models)) {
    col_ba  <- paste0("pheno_", nm)
    col_adv <- paste0("pheno_", nm, "_advance")

    model <- bundle$subclock_models[[nm]]
    panel <- bundle$panels$subclocks[[nm]]

    if (is.null(model) || is.null(panel)) {
      message("  ! Sub-clock '", nm, "': no model or panel in bundle")
      next
    }

    tryCatch({
      ba <- score_phenoage(patient_data, model, panel)
      n_ok <- sum(!is.na(ba))
      if (n_ok > 0) {
        result[[col_ba]]  <- ba
        result[[col_adv]] <- ba - patient_data$age
        message(">> pheno_", nm, ": ", n_ok, " scored")
      } else {
        message(">> pheno_", nm, ": 0 scored (check biomarker availability)")
      }
    }, error = function(e) {
      message("  ! Sub-clock '", nm, "' failed: ", conditionMessage(e))
    })
  }

  # ---- 4. Residualize using saved coefficients ----------------------------
  if (!is.null(bundle$resid_coefs)) {
    for (i in seq_len(nrow(bundle$resid_coefs))) {
      v  <- bundle$resid_coefs$variable[i]
      b0 <- bundle$resid_coefs$intercept[i]
      b1 <- bundle$resid_coefs$slope[i]

      if (v %in% names(result)) {
        predicted <- b0 + b1 * result$age
        result[[paste0(v, "_resid")]] <- result[[v]] - predicted
      }
    }
  }

  # ---- 5. Z-score using saved population stats ----------------------------
  if (!is.null(bundle$pop_stats)) {
    for (i in seq_len(nrow(bundle$pop_stats))) {
      v  <- bundle$pop_stats$variable[i]
      mu <- bundle$pop_stats$mean[i]
      s  <- bundle$pop_stats$sd[i]

      resid_col <- paste0(v, "_resid")
      if (resid_col %in% names(result) && s > 0) {
        result[[paste0(v, "_resid_z")]] <- (result[[resid_col]] - mu) / s
      }
    }
  }

  n_clocks <- sum(grepl("_advance$", names(result)))
  message(">> Scored ", nrow(result), " patients across ", n_clocks, " clocks")

  result
}

# ── Clinical report table ─────────────────────────────────────────────────

clinical_report_table <- function(scored_data, bundle = NULL) {

  adv_cols <- grep("_advance$", names(scored_data), value = TRUE)

  rows <- list()
  for (ac in adv_cols) {
    clock_name <- sub("_advance$", "", ac)
    ba_col     <- clock_name
    resid_col  <- paste0(ac, "_resid")
    z_col      <- paste0(ac, "_resid_z")

    for (i in seq_len(nrow(scored_data))) {
      ba  <- if (ba_col %in% names(scored_data)) scored_data[[ba_col]][i] else NA
      adv <- scored_data[[ac]][i]
      res <- if (resid_col %in% names(scored_data)) scored_data[[resid_col]][i] else NA
      zsc <- if (z_col %in% names(scored_data)) scored_data[[z_col]][i] else NA

      risk <- if (is.na(zsc)) {
        NA_character_
      } else if (zsc <= -1) {
        "Low"
      } else if (zsc <= 1) {
        "Average"
      } else {
        "Elevated"
      }

      short  <- sub("^pheno_|^phenoage_", "", clock_name)
      system <- if (short %in% names(PANEL_SYSTEM)) PANEL_SYSTEM[short] else "Global"

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
