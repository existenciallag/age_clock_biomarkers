###############################################################################
# run_clinical.R — Score new lab patients using a trained deployment bundle
#
# This script is designed for clinical laboratory deployment:
#   1. Loads the pre-trained model bundle (from run_pipeline.R)
#   2. Reads new patient lab data (CSV)
#   3. Maps biomarker names to NHANES conventions
#   4. Scores biological age across all available clocks
#   5. Computes residuals and Z-scores against the NHANES reference
#   6. Produces a clinical report table
#
# NO RETRAINING occurs — all model parameters come from the bundle.
# This is appropriate when clinical outcomes (mortality) are not
# available in the local population, per the validation design.
#
# Usage:
#   1. Edit PATIENT_DATA_PATH and NAME_MAP below
#   2. source("run_clinical.R")
#
# Output:
#   - clinical_report.csv : per-patient, per-clock report
#   - scored_data.csv     : wide-format scored dataset
###############################################################################

cat("==========================================================\n")
cat("  BIOLOGICAL AGE — CLINICAL SCORING\n")
cat("  Applying pre-trained NHANES model to local lab data\n")
cat("==========================================================\n\n")

# ── 0. Load dependencies ─────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
})

# Load toolkit modules
source("R/config.R")
source("R/export.R")
source("R/clinical_score.R")

# ── 1. Configuration ─────────────────────────────────────────────────────

# Path to the deployment bundle (created by run_pipeline.R)
BUNDLE_PATH <- "bioage_deployment_bundle.rds"

# Path to patient lab data (CSV)
# Expected columns: patient_id, age, gender, + biomarker columns
PATIENT_DATA_PATH <- "patient_data.csv"

# Biomarker name mapping: your lab column names → NHANES names
# Modify this to match your LIMS export format.
# Only include names that DIFFER from NHANES conventions.
# If your lab already uses NHANES names, set NAME_MAP <- NULL.
NAME_MAP <- c(
  # "your_lab_name" = "nhanes_name"
  # Example:
  # "ALB_g_L"    = "albumin_gL",
  # "ALP_U_L"    = "alp",
  # "CRP_log"    = "lncrp",
  # "CHOL_total" = "totchol",
  # "CREAT_log"  = "lncreat_umol",
  # "HbA1c"      = "hba1c",
  # "SBP_mmHg"   = "sbp",
  # "BUN_mg_dL"  = "bun",
  # "Uric_acid"  = "uap",
  # "Lymph_pct"  = "lymph",
  # "MCV_fL"     = "mcv",
  # "WBC_K_uL"   = "wbc"
)

# ── 2. Load bundle ───────────────────────────────────────────────────────

message(">> Loading deployment bundle...")
bundle <- load_bundle(BUNDLE_PATH)

cat("\nTrained clocks available:\n")
cat("  ", paste(bundle$meta$clocks_trained, collapse = "\n   "), "\n")

cat("\nPopulation reference (NHANES IV):\n")
cat("  N subjects:  ", bundle$meta$n_subjects, "\n")
cat("  R version:   ", bundle$meta$R_version, "\n")
cat("  Created:     ", format(bundle$meta$created), "\n\n")

# ── 3. Load patient data ─────────────────────────────────────────────────

message(">> Loading patient data from: ", PATIENT_DATA_PATH)

if (!file.exists(PATIENT_DATA_PATH)) {
  stop("Patient data file not found: ", PATIENT_DATA_PATH, "\n",
       "Create a CSV with columns: patient_id, age, gender, + biomarkers")
}

patients <- read.csv(PATIENT_DATA_PATH, stringsAsFactors = FALSE)
message(">> Loaded ", nrow(patients), " patients, ",
        ncol(patients), " columns")

# Quick audit: which expected biomarkers are present?
all_expected <- unique(unlist(bundle$panels))
if (!is.null(NAME_MAP) && length(NAME_MAP) > 0) {
  # Check after mapping
  test_names <- names(patients)
  for (i in seq_along(NAME_MAP)) {
    idx <- match(names(NAME_MAP)[i], test_names)
    if (!is.na(idx)) test_names[idx] <- NAME_MAP[i]
  }
  present  <- intersect(test_names, all_expected)
  missing  <- setdiff(all_expected, test_names)
} else {
  present  <- intersect(names(patients), all_expected)
  missing  <- setdiff(all_expected, names(patients))
}

cat("\nBiomarker audit:\n")
cat("  Present: ", length(present), "/", length(all_expected), "\n")
if (length(missing) > 0) {
  cat("  Missing: ", paste(missing, collapse = ", "), "\n")
  cat("  (Clocks needing these biomarkers will return NA)\n")
}

# ── 4. Score patients ────────────────────────────────────────────────────

message("\n>> Scoring patients...")

scored <- score_patients(
  patient_data = patients,
  bundle       = bundle,
  name_map     = if (length(NAME_MAP) > 0) NAME_MAP else NULL
)

# ── 5. Generate clinical report ──────────────────────────────────────────

message(">> Generating clinical report table...")

report <- clinical_report_table(scored, bundle)

# ── 6. Display results ───────────────────────────────────────────────────

cat("\n")
cat("==========================================================\n")
cat("  CLINICAL REPORT — PER-PATIENT SUMMARY\n")
cat("==========================================================\n\n")

print(report, row.names = FALSE)

# ── 7. Save outputs ──────────────────────────────────────────────────────

write.csv(report, "clinical_report.csv", row.names = FALSE)
message("\n>> Saved: clinical_report.csv")

write.csv(scored, "scored_data.csv", row.names = FALSE)
message(">> Saved: scored_data.csv")

# ── Summary ──────────────────────────────────────────────────────────────

n_clocks_scored <- sum(grepl("_advance$", names(scored)))

cat("\n")
cat("==========================================================\n")
cat("  SCORING COMPLETE\n")
cat("==========================================================\n")
cat("  Patients scored:  ", nrow(scored), "\n")
cat("  Clocks applied:   ", n_clocks_scored, "\n")
cat("  Report rows:      ", nrow(report), "\n")
cat("  Output files:\n")
cat("    - clinical_report.csv (long format, per clock)\n")
cat("    - scored_data.csv     (wide format, all scores)\n")
cat("==========================================================\n")

# Interpretation guide
cat("\n--- Interpretation Guide ---\n")
cat("  advancement > 0 : Biologically older than chronological age\n")
cat("  advancement < 0 : Biologically younger than chronological age\n")
cat("  z_score > +1    : Elevated aging (>1 SD above NHANES reference)\n")
cat("  z_score < -1    : Reduced aging  (<1 SD below NHANES reference)\n")
cat("  risk_level:\n")
cat("    Low     = z-score <= -1  (protective)\n")
cat("    Average = -1 < z-score <= 1\n")
cat("    Elevated = z-score > 1  (accelerated aging)\n")
