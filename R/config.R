###############################################################################
# config.R — Panel definitions, system mappings, and toolkit constants
#
# This file is the single source of truth for biomarker panel composition.
# To add or remove a biomarker, edit the relevant list here and re-run
# the pipeline.  No other file needs to change.
#
# Conventions
#   - PhenoAge Levine 2018 uses 9 biomarkers + age (age added internally)
#   - BioAge variable names match NHANES III/IV column names
#   - Sub-clocks are always trained via phenoage_nhanes()
###############################################################################

# ── 1. Canonical biomarker panels ──────────────────────────────────────────

#' Original Levine 2018 PhenoAge: 9 biomarkers
#' (chronological age is added by the algorithm, not listed here)
#'
#' Units in NHANES:
#'   albumin    = g/dL
#'   alp        = U/L  (alkaline phosphatase)
#'   glucose    = mg/dL (serum glucose)
#'   lncrp      = log(CRP mg/dL)
#'   lncreat    = log(creatinine mg/dL)
#'   lymph      = lymphocyte %
#'   mcv        = fL   (mean cell volume)
#'   rdw        = %    (red cell distribution width)
#'   wbc        = 1000 cells/uL
BIOMARKERS_PHENO_LEVINE <- c(
  "albumin", "alp", "glucose", "lncrp", "lncreat",
  "lymph", "mcv", "rdw", "wbc"
)

#' Alias: keep backward compatibility
BIOMARKERS_PHENO_ORIG <- BIOMARKERS_PHENO_LEVINE

#' Panel shared by KDM and HD clocks (same BioAge internal names)
#' NOTE: KDM/HD use a broader 12-biomarker set
BIOMARKERS_KDM_HD <- c(
  "albumin", "alp", "lncrp", "totchol",
  "lncreat", "hba1c", "sbp", "bun",
  "uap", "lymph", "mcv", "wbc"
)

#' Global PhenoAge panel — same as Levine by default.
PANEL_GLOBAL <- BIOMARKERS_PHENO_LEVINE

# ── 2. Sub-clock panel definitions ─────────────────────────────────────────
#' Each entry maps a short name → character vector of biomarkers.
#' All sub-clocks are trained as PhenoAge models on NHANES III.

PANELS_SUB <- list(
  hepatic_enzime_insulin    = c("alp", "ggt", "insulin"),
  hepatic_lipid             = c("trig", "totchol"),
  hema_integrated           = c("rdw", "mcv", "rbc", "wbc", "lymph"),
  micronutrient_methylation = c("vitaminB12", "hba1c", "rdw"),
  renal_A                   = c("lncreat", "bun"),
  # Levine-9 minus glucose — for cohorts lacking glucose data
  levine_no_glucose         = c("albumin", "alp", "lncrp", "lncreat",
                                 "lymph", "mcv", "rdw", "wbc")
)

# ── 3. System-level grouping ───────────────────────────────────────────────
#' Maps each sub-clock name to a physiological system label.

PANEL_SYSTEM <- c(
  hepatic_enzime_insulin    = "Hepatic",
  hepatic_lipid             = "Hepatic",
  hema_integrated           = "Hematologic",
  micronutrient_methylation = "Micronutrient / Methylation",
  renal_A                   = "Renal",
  levine_no_glucose         = "Global (no Glucose)"
)

# ── 4. Published Levine 2018 PhenoAge coefficients ────────────────────────
#' These are the EXACT published coefficients from Levine (2018) Table S2.
#' Used for direct scoring without a trained fit object.
#'
#' Mortality score: xb = sum(beta_i * x_i)
#' Then: M = 1 - exp(-exp(xb) * (exp(120*gamma) - 1) / gamma)
#' Then: PhenoAge = 141.50225 + ln(-0.00553 * ln(1-M)) / 0.090165

LEVINE_2018_COEFS <- c(
  intercept = -19.9067,
  albumin   =  -0.0336,   # g/dL
  alp       =   0.00188,  # U/L
  glucose   =   0.00553,  # mg/dL (NOTE: called "serum glucose" in paper)
  lncrp     =   0.0954,   # log(CRP mg/dL)
  lncreat   =   0.0095,   # log(creatinine mg/dL)
  lymph     =  -0.0120,   # lymphocyte %
  mcv       =   0.0268,   # fL
  rdw       =   0.3306,   # %
  wbc       =   0.0554,   # 1000 cells/uL
  age       =   0.0804    # years
)

LEVINE_GAMMA_MORT   <- 0.0076927   # Gompertz mortality rate
LEVINE_PHENOAGE_A   <- 141.50225   # PhenoAge inversion intercept
LEVINE_PHENOAGE_B   <- -0.00553    # PhenoAge inversion factor
LEVINE_PHENOAGE_C   <- 0.090165    # PhenoAge inversion slope

# ── 5. QC thresholds ──────────────────────────────────────────────────────

QC_MIN_N <- 2000
QC_MIN_SD <- 0
COX_MIN_N <- 200

# ── 6. Mortality column detection ─────────────────────────────────────────

MORT_STATUS_PATTERNS <- c("^status$", "mortstat", "dies", "dead", "outcome")
MORT_TIME_PATTERNS   <- c("^time$", "permth_exm", "permth_int", "followup", "surv")

# ── 7. Health outcome columns (from BioAge NHANES IV) ───────────────────

HEALTH_OUTCOMES <- c("health", "adl", "lnwalk", "grip_scaled")
SES_VARS <- c("edu", "annual_income", "poverty_ratio")

# ── 8. Palette defaults ───────────────────────────────────────────────────

HEATMAP_PALETTE <- c("blue", "white", "red")
FOREST_COLOUR   <- "steelblue"
KM_COLOURS      <- c("darkgreen", "gray40", "darkred")
