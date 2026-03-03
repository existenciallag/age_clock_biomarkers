###############################################################################
# config.R — Panel definitions, system mappings, and toolkit constants
#
# This file is the single source of truth for biomarker panel composition.
# To add or remove a biomarker, edit the relevant list here and re-run
# the pipeline.  No other file needs to change.
#
# Conventions
#   - PhenoAge-style clocks use albumin_gL / lncreat_umol (NHANES units)
#   - KDM / HD clocks use albumin / lncreat (BioAge internal units)
#   - Sub-clocks are always trained via phenoage_nhanes()
###############################################################################

# ── 1. Canonical biomarker panels ──────────────────────────────────────────

#' Original 12-biomarker PhenoAge panel (Levine et al.)
BIOMARKERS_PHENO_ORIG <- c(
  "albumin_gL", "alp", "lncrp", "totchol",
  "lncreat_umol", "hba1c", "sbp", "bun",
  "uap", "lymph", "mcv", "wbc"
)

#' Panel shared by KDM and HD clocks (BioAge internal names)
BIOMARKERS_KDM_HD <- c(
  "albumin", "alp", "lncrp", "totchol",
  "lncreat", "hba1c", "sbp", "bun",
  "uap", "lymph", "mcv", "wbc"
)

#' Global PhenoAge panel — same as original by default.
#' Override to test reduced global panels.
PANEL_GLOBAL <- BIOMARKERS_PHENO_ORIG

# ── 2. Sub-clock panel definitions ─────────────────────────────────────────
#' Each entry maps a short name → character vector of biomarkers.
#' All sub-clocks are trained as PhenoAge models on NHANES III.

PANELS_SUB <- list(
  hepatic_enzime_insulin    = c("alp", "ggt", "insulin"),
  hepatic_lipid             = c("trig", "totchol"),
  hema_integrated           = c("rdw", "mcv", "rbc", "wbc", "lymph"),
  micronutrient_methylation = c("vitaminB12", "hba1c", "rdw"),
  renal_A                   = c("lncreat_umol", "bun"),
  # Full PhenoAge panel minus SBP — for cohorts lacking blood pressure data
  full_no_sbp               = c("albumin_gL", "alp", "lncrp", "totchol",
                                 "lncreat_umol", "hba1c", "bun",
                                 "uap", "lymph", "mcv", "wbc")
)

# ── 3. System-level grouping ───────────────────────────────────────────────
#' Maps each sub-clock name to a physiological system label.
#' Used for forest-plot grouping and clinical report sections.

PANEL_SYSTEM <- c(
  hepatic_enzime_insulin    = "Hepatic",
  hepatic_lipid             = "Hepatic",
  hema_integrated           = "Hematologic",
  micronutrient_methylation = "Micronutrient / Methylation",
  renal_A                   = "Renal",
  full_no_sbp               = "Global (no SBP)"
)

# ── 4. QC thresholds ──────────────────────────────────────────────────────

#' Minimum number of finite observations to consider a clock valid for
#' survival analysis.
QC_MIN_N <- 2000

#' Minimum standard deviation — guards against near-constant clocks.
QC_MIN_SD <- 0

#' Minimum sample size for a Cox model to be reported.
COX_MIN_N <- 200

# ── 5. Mortality column detection ─────────────────────────────────────────
#' Patterns used to auto-detect mortality status and follow-up columns.
#' BioAge already renames mortstat -> status, permth_exm -> time.
#' Checked in order; first match wins.

MORT_STATUS_PATTERNS <- c("^status$", "mortstat", "dies", "dead", "outcome")
MORT_TIME_PATTERNS   <- c("^time$", "permth_exm", "permth_int", "followup", "surv")

# ── 6. Health outcome columns (from BioAge NHANES IV) ───────────────────
#' These are used by BioAge::table_health() and our health correlation step.
HEALTH_OUTCOMES <- c("health", "adl", "lnwalk", "grip_scaled")

#' SES variables for BioAge::table_ses()
SES_VARS <- c("edu", "annual_income", "poverty_ratio")

# ── 7. Palette defaults ───────────────────────────────────────────────────
HEATMAP_PALETTE <- c("blue", "white", "red")
FOREST_COLOUR   <- "steelblue"
KM_COLOURS      <- c("darkgreen", "gray40", "darkred")
