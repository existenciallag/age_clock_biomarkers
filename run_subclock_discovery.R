###############################################################################
# run_subclock_discovery.R — Sub-clock feasibility, fast/slow ager profiling,
#                            and new panel candidate analysis
#
# PURPOSE
# -------
# 1. Cross-reference Argentine cohort biomarkers vs NHANES IV availability
# 2. Identify which NEW sub-clocks could be trained (shared analytes with
#    sufficient coverage in BOTH populations)
# 3. Profile fast vs slow biological agers using ALL available biomarkers
#    (not just clock components) — characterizing the "accelerated aging
#    phenotype" in this population
# 4. Evaluate functional medicine reference ranges as potential sensitivity
#    enhancement for aging classification
#
# REQUIREMENTS
# ------------
#   - clinical_recentered.csv (from run_recenter.R)
#   - phenoage_master_PARTE1.csv (+ optional PARTE2.csv)
#   - bioage_deployment_bundle.rds (from run_pipeline.R)
#
# USAGE:  source("run_subclock_discovery.R")
###############################################################################

cat("==========================================================\n")
cat("  SUB-CLOCK DISCOVERY & FAST/SLOW AGER PROFILING\n")
cat("  Shared analytes | New panel candidates | Ager phenotyping\n")
cat("==========================================================\n\n")

# ── 0. Setup ──────────────────────────────────────────────────────────────

PROJECT_ROOT <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) getwd())
PROJECT_ROOT <- normalizePath(PROJECT_ROOT)
setwd(PROJECT_ROOT)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(BioAge)
})

source(file.path(PROJECT_ROOT, "R", "config.R"), local = FALSE)
source(file.path(PROJECT_ROOT, "R", "export.R"), local = FALSE)

OUTPUT_DIR <- file.path(PROJECT_ROOT, "discovery_output")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────

message(">> Loading data...")

# Recentered scores
recenter_path <- file.path(PROJECT_ROOT, "recenter_output", "clinical_recentered.csv")
scored_path   <- file.path(PROJECT_ROOT, "clinical_scored_lite.csv")
bundle_path   <- file.path(PROJECT_ROOT, "bioage_deployment_bundle.rds")

if (file.exists(recenter_path)) {
  scores <- read.csv(recenter_path, stringsAsFactors = FALSE)
  message(">> Loaded recentered data: ", nrow(scores), " subjects")
} else if (file.exists(scored_path)) {
  scores <- read.csv(scored_path, stringsAsFactors = FALSE)
  message(">> Loaded scored data (not recentered): ", nrow(scores), " subjects")
} else {
  stop("No scored data found. Run run_clinical_lite.R and run_recenter.R first.")
}

# Raw clinical
f1 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE2.csv")
raw <- read.csv(f1, stringsAsFactors = FALSE)
if (file.exists(f2)) raw <- rbind(raw, read.csv(f2, stringsAsFactors = FALSE))

raw$patient_id <- raw$Protocolo
raw$age <- as.numeric(raw$age)
raw <- raw[!is.na(raw$age) & raw$age >= 1 & raw$age <= 120, ]
message(">> Raw clinical records: ", nrow(raw))

# Bundle
bundle <- NULL
if (file.exists(bundle_path)) {
  bundle <- readRDS(bundle_path)
  message(">> Bundle loaded: ", length(bundle$meta$clocks_trained), " clocks")
}

# NHANES IV reference
message(">> Extracting NHANES IV reference...")
nhanes <- tryCatch({
  res <- phenoage_nhanes(BIOMARKERS_PHENO_LEVINE)
  res$data
}, error = function(e) {
  message(">> NHANES IV extraction failed: ", e$message)
  NULL
})
if (!is.null(nhanes)) message(">> NHANES IV: ", nrow(nhanes), " subjects")

###############################################################################
#   SECTION 1: COMPLETE BIOMARKER INVENTORY — CROSS-REFERENCE
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 1: BIOMARKER CROSS-REFERENCE                  #\n")
cat("##########################################################\n\n")

# ── Argentine cohort: all numeric columns with > 0.5% coverage ──

arg_all_vars <- names(raw)
arg_numeric <- sapply(arg_all_vars, function(v) {
  x <- suppressWarnings(as.numeric(raw[[v]]))
  sum(!is.na(x)) / nrow(raw)
})
arg_numeric <- arg_numeric[arg_numeric > 0.005]
# Exclude metadata columns
exclude_cols <- c("Protocolo", "source_file", "date_admission", "date_result",
                  "note", "note_pending", "patient_id",
                  grep("eligible_|morphology|lipemia|flag|log_", names(arg_numeric),
                       value = TRUE))
arg_biomarkers <- setdiff(names(arg_numeric), exclude_cols)

cat("--- Argentine Cohort: Available Biomarkers ---\n\n")
arg_inv <- data.frame(
  variable = arg_biomarkers,
  n        = sapply(arg_biomarkers, function(v) {
    sum(!is.na(suppressWarnings(as.numeric(raw[[v]]))))
  }),
  stringsAsFactors = FALSE
)
arg_inv$pct <- round(100 * arg_inv$n / nrow(raw), 1)
arg_inv <- arg_inv[order(-arg_inv$pct), ]
rownames(arg_inv) <- NULL

# Tier classification
arg_inv$tier <- ifelse(arg_inv$pct >= 60, "A (>=60%)",
                ifelse(arg_inv$pct >= 20, "B (20-60%)",
                ifelse(arg_inv$pct >= 5,  "C (5-20%)", "D (<5%)")))

print(arg_inv, row.names = FALSE)
write.csv(arg_inv, file.path(OUTPUT_DIR, "argentine_biomarker_inventory.csv"),
          row.names = FALSE)

# ── NHANES IV: all biomarker columns ──

# Known NHANES IV biomarkers (from BioAge package data-raw/nhanes_all.R)
nhanes_biomarkers <- c(
  # Core blood chemistry
  "albumin", "alp", "bun", "creat", "glucose", "uap", "ttbl", "bap",
  "lncrp", "lncreat", "lnbun", "lnuap", "lnalp",
  "albumin_gL", "creat_umol", "lncreat_umol", "glucose_mmol", "glucose_fasting",
  # Hematology
  "wbc", "rbc", "mcv", "rdw", "lymph", "neut", "basopa", "eosnpa", "monopa",
  # Lipids
  "totchol", "hdl", "ldl", "trig",
  # Metabolic
  "hba1c", "lnhba1c", "insulin", "ggt",
  # Vitamins
  "vitaminA", "vitaminE", "vitaminB12", "vitaminC",
  # Other
  "cadmium", "cyst",
  # Physical
  "sbp", "dbp", "meanbp", "pulse", "bmi", "fev", "fev_1000",
  "waist", "height", "weight",
  # Body composition
  "grip_scaled", "grip_r", "grip_l"
)

# Check which actually exist in the NHANES data object
nhanes_available <- character(0)
if (!is.null(nhanes)) {
  nhanes_available <- intersect(nhanes_biomarkers, names(nhanes))
  nhanes_inv <- data.frame(
    variable = nhanes_available,
    n = sapply(nhanes_available, function(v) sum(!is.na(nhanes[[v]]))),
    stringsAsFactors = FALSE
  )
  nhanes_inv$pct <- round(100 * nhanes_inv$n / nrow(nhanes), 1)
  nhanes_inv <- nhanes_inv[order(-nhanes_inv$pct), ]
  rownames(nhanes_inv) <- NULL

  cat("\n--- NHANES IV: Available Biomarkers ---\n\n")
  print(nhanes_inv, row.names = FALSE)
  write.csv(nhanes_inv, file.path(OUTPUT_DIR, "nhanes_biomarker_inventory.csv"),
            row.names = FALSE)
}

# ── Cross-reference: Argentine → NHANES name mapping ──

cat("\n--- Cross-Reference: Shared Biomarkers ---\n\n")

# Manual mapping: Argentine column name → NHANES column name
# (where a direct or transformable correspondence exists)
arg_to_nhanes <- c(
  # Direct matches (same name)
  albumin         = "albumin",
  alp             = "alp",
  glucose         = "glucose",
  creatinine      = "creat",
  crp             = "lncrp",       # needs log transform
  # Hematology
  wbc             = "wbc",
  rbc             = "rbc",
  mcv             = "mcv",
  rdw             = "rdw",
  lymphocyte      = "lymph",
  neutrophils     = "neut",
  basophils       = "basopa",
  eosinophils     = "eosnpa",
  monocytes       = "monopa",
  # Lipids
  cholesterol     = "totchol",
  hdl             = "hdl",
  ldl             = "ldl",
  triglycerides   = "trig",
  # Metabolic
  hba1c           = "hba1c",
  insulin         = "insulin",
  ggt             = "ggt",
  # Vitamins
  vitamin_b12     = "vitaminB12",
  # Renal
  bun             = "bun",
  uric_acid       = "uap"
)

# Argentine-only biomarkers (no NHANES equivalent)
arg_only <- c("hemoglobin", "hematocrit", "platelets", "mch", "mchc",
              "non_hdl", "lp_a", "ast", "alt", "ldh",
              "ferritin", "transferrin", "transferrin_saturation",
              "fibrinogen", "homocysteine", "TSH", "free_t4", "vitamin_d",
              "band_neutrophils", "metamyelocytes", "myelocytes",
              "lymphocyte_abs")

# Build cross-reference table
xref <- data.frame(
  arg_name    = names(arg_to_nhanes),
  nhanes_name = unname(arg_to_nhanes),
  stringsAsFactors = FALSE
)
xref$arg_n <- sapply(xref$arg_name, function(v) {
  if (v %in% names(raw)) sum(!is.na(suppressWarnings(as.numeric(raw[[v]])))) else 0
})
xref$arg_pct <- round(100 * xref$arg_n / nrow(raw), 1)

xref$nhanes_n <- 0
xref$nhanes_pct <- 0
if (!is.null(nhanes)) {
  for (i in seq_len(nrow(xref))) {
    nv <- xref$nhanes_name[i]
    if (nv %in% names(nhanes)) {
      xref$nhanes_n[i]   <- sum(!is.na(nhanes[[nv]]))
      xref$nhanes_pct[i] <- round(100 * xref$nhanes_n[i] / nrow(nhanes), 1)
    }
  }
}

# Flag usability for clock training (need >2000 in NHANES, >1000 in Argentina)
xref$trainable <- xref$nhanes_n >= 2000 & xref$arg_n >= 1000
xref <- xref[order(-xref$arg_pct), ]
rownames(xref) <- NULL

print(xref, row.names = FALSE)
write.csv(xref, file.path(OUTPUT_DIR, "biomarker_crossreference.csv"), row.names = FALSE)

# Count trainable
n_trainable <- sum(xref$trainable)
cat(sprintf("\n>> Trainable biomarkers (NHANES >= 2000 & Argentina >= 1000): %d / %d\n",
            n_trainable, nrow(xref)))

###############################################################################
#   SECTION 2: NEW SUB-CLOCK CANDIDATE PANELS
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 2: NEW SUB-CLOCK CANDIDATES                   #\n")
cat("##########################################################\n\n")

# Trainable biomarkers (mapped to NHANES names)
trainable <- xref[xref$trainable, ]
trainable_nhanes <- trainable$nhanes_name

cat("Trainable biomarkers (NHANES names):\n")
cat("  ", paste(trainable_nhanes, collapse = ", "), "\n\n")

# ── Existing sub-clock panels ──
cat("--- Existing sub-clock panels ---\n")
for (nm in names(PANELS_SUB)) {
  bm <- PANELS_SUB[[nm]]
  in_both <- bm %in% trainable_nhanes
  cat(sprintf("  %-30s [%s] — %d/%d biomarkers trainable\n",
              nm, paste(bm, collapse = ", "),
              sum(in_both), length(bm)))
}

# ── Proposed NEW panels based on shared analytes ──

cat("\n--- Proposed NEW Sub-Clock Panels ---\n\n")

# Candidate panels organized by physiological system
# Only include biomarkers that are trainable (shared + sufficient N)
candidate_panels <- list()

# 1. CBC-Extended: full blood count beyond current hema_integrated
cbc_extended_bm <- intersect(
  c("rdw", "mcv", "rbc", "wbc", "lymph", "neut", "basopa", "eosnpa", "monopa"),
  trainable_nhanes
)
if (length(cbc_extended_bm) >= 3) {
  candidate_panels$cbc_extended <- cbc_extended_bm
}

# 2. Lipid-Metabolic: lipid panel + metabolic markers
lipid_meta_bm <- intersect(
  c("totchol", "hdl", "ldl", "trig", "hba1c"),
  trainable_nhanes
)
if (length(lipid_meta_bm) >= 3) {
  candidate_panels$lipid_metabolic <- lipid_meta_bm
}

# 3. Metabolic-Glycemic: insulin resistance markers
glycemic_bm <- intersect(
  c("glucose", "hba1c", "insulin"),
  trainable_nhanes
)
if (length(glycemic_bm) >= 2) {
  candidate_panels$glycemic <- glycemic_bm
}

# 4. Inflammatory-Hematologic: WBC differential + CRP
inflam_bm <- intersect(
  c("wbc", "neut", "lymph", "lncrp", "monopa"),
  trainable_nhanes
)
if (length(inflam_bm) >= 3) {
  candidate_panels$inflammatory <- inflam_bm
}

# 5. Lipid-Only: cholesterol + triglycerides + HDL
lipid_only_bm <- intersect(
  c("totchol", "hdl", "trig"),
  trainable_nhanes
)
if (length(lipid_only_bm) >= 2) {
  candidate_panels$lipid_triad <- lipid_only_bm
}

# 6. Metabolic-Renal: renal + metabolic crossover
renal_meta_bm <- intersect(
  c("creat", "bun", "uap", "albumin", "glucose"),
  trainable_nhanes
)
# Only if there's at least 2 with reasonable Argentine coverage
renal_meta_n <- sapply(renal_meta_bm, function(nv) {
  av <- xref$arg_name[xref$nhanes_name == nv]
  if (length(av) > 0 && av %in% names(raw))
    sum(!is.na(suppressWarnings(as.numeric(raw[[av]]))))
  else 0
})
renal_meta_bm <- renal_meta_bm[renal_meta_n >= 500]
if (length(renal_meta_bm) >= 2) {
  candidate_panels$renal_metabolic <- renal_meta_bm
}

# 7. CBC-Minimal: most available hematologic markers only
cbc_minimal_bm <- intersect(c("rdw", "mcv", "wbc", "lymph"), trainable_nhanes)
if (length(cbc_minimal_bm) >= 3) {
  candidate_panels$cbc_minimal <- cbc_minimal_bm
}

# Print candidates with estimated Argentine coverage
for (nm in names(candidate_panels)) {
  bm <- candidate_panels[[nm]]
  # Estimate coverage: patients with ALL biomarkers
  # Map back to Argentine column names
  arg_names <- sapply(bm, function(nv) {
    hit <- xref$arg_name[xref$nhanes_name == nv]
    if (length(hit) > 0) hit[1] else NA
  })
  arg_names <- arg_names[!is.na(arg_names)]

  if (length(arg_names) == length(bm)) {
    # Count complete cases
    complete <- complete.cases(raw[, arg_names, drop = FALSE])
    n_complete <- sum(complete)
    pct_complete <- round(100 * n_complete / nrow(raw), 1)
  } else {
    n_complete <- NA
    pct_complete <- NA
  }

  cat(sprintf("  %-25s [%s]\n", nm, paste(bm, collapse = ", ")))
  cat(sprintf("  %-25s Estimated N (Argentina) = %s (%s%%)\n",
              "", format(n_complete, big.mark = ","), pct_complete))

  # Check vs existing panels — any overlap?
  for (enm in names(PANELS_SUB)) {
    overlap <- intersect(bm, PANELS_SUB[[enm]])
    if (length(overlap) > 0) {
      cat(sprintf("  %-25s Overlaps with %s: %s\n",
                  "", enm, paste(overlap, collapse = ", ")))
    }
  }
  cat("\n")
}

write.csv(
  data.frame(
    panel = rep(names(candidate_panels), sapply(candidate_panels, length)),
    biomarker_nhanes = unlist(candidate_panels),
    stringsAsFactors = FALSE
  ),
  file.path(OUTPUT_DIR, "candidate_panels.csv"),
  row.names = FALSE
)

###############################################################################
#   SECTION 3: FAST vs SLOW AGER PROFILING
###############################################################################

cat("##########################################################\n")
cat("#  SECTION 3: FAST vs SLOW AGER PROFILING                #\n")
cat("##########################################################\n\n")

# Use the best recentered clock (highest coverage)
# Prefer sexadj_z columns (adjusted for age + sex) over local_z (age only)
z_cols <- grep("_sexadj_z$", names(scores), value = TRUE)
if (length(z_cols) == 0) {
  z_cols <- grep("_local_z$", names(scores), value = TRUE)
  message(">> NOTE: No sex-adjusted Z-scores found, falling back to local_z (age-only)")
}
if (length(z_cols) == 0) {
  # Fallback to raw advancement
  z_cols <- grep("_advance$", names(scores), value = TRUE)
}

# Pick the one with most data
z_coverage <- sapply(z_cols, function(zc) sum(!is.na(scores[[zc]])))
primary_z <- z_cols[which.max(z_coverage)]
message(">> Primary Z-score for profiling: ", primary_z,
        " (n = ", max(z_coverage), ")")
message(">> This Z-score is residualized for age",
        ifelse(grepl("sexadj", primary_z), " AND sex", " only (not sex-adjusted)"))

# Merge scores with raw biomarkers
df <- merge(scores, raw[, c("patient_id", setdiff(names(raw),
            c(names(scores), "Protocolo", "source_file",
              "date_admission", "date_result")))],
            by = "patient_id", all.x = TRUE)

# Convert all potential biomarkers to numeric
all_biomarkers <- c(
  # CBC
  "hemoglobin", "hematocrit", "platelets", "mch", "mchc",
  "neutrophils", "eosinophils", "basophils", "monocytes",
  "lymphocyte_abs", "band_neutrophils",
  # Lipids
  "cholesterol", "hdl", "ldl", "non_hdl", "triglycerides", "lp_a",
  # Metabolic
  "glucose", "hba1c", "insulin",
  # Liver
  "ast", "alt", "ggt", "ldh", "albumin", "alp",
  # Iron
  "ferritin", "transferrin", "transferrin_saturation",
  # Renal
  "creatinine", "bun", "uric_acid", "eGFR",
  # Inflammatory
  "crp", "fibrinogen", "homocysteine",
  # Thyroid
  "TSH", "free_t4",
  # Vitamins
  "vitamin_d", "vitamin_b12",
  # Hematology (from scored)
  "wbc", "rbc", "mcv", "rdw", "lymphocyte"
)
all_biomarkers <- intersect(all_biomarkers, names(df))

for (bm in all_biomarkers) {
  df[[bm]] <- suppressWarnings(as.numeric(df[[bm]]))
}

# ── Classify fast vs slow agers ──

ok <- !is.na(df[[primary_z]]) & is.finite(df[[primary_z]])
df_valid <- df[ok, ]
n_valid <- nrow(df_valid)
message(">> Valid subjects for profiling: ", n_valid)

# Quintile-based classification
quants <- quantile(df_valid[[primary_z]], probs = c(0.10, 0.25, 0.75, 0.90),
                   na.rm = TRUE)

df_valid$ager_group <- cut(df_valid[[primary_z]],
  breaks = c(-Inf, quants[1], quants[2], quants[3], quants[4], Inf),
  labels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")
)

cat(sprintf(">> Ager groups (based on %s):\n", primary_z))
tbl <- table(df_valid$ager_group)
for (g in names(tbl)) {
  cat(sprintf("   %-15s n = %6d (%5.1f%%)\n", g, tbl[g], 100 * tbl[g] / n_valid))
}

# ── 3a. Compare biomarker profiles across ager groups ──

cat("\n--- Biomarker Profiles by Ager Group ---\n\n")

# Only use biomarkers with >= 50 observations in the valid set
usable_bm <- character(0)
for (bm in all_biomarkers) {
  n_ok <- sum(!is.na(df_valid[[bm]]) & is.finite(df_valid[[bm]]))
  if (n_ok >= 50) usable_bm <- c(usable_bm, bm)
}
message(">> Usable biomarkers for profiling: ", length(usable_bm))

profile_results <- list()

for (bm in usable_bm) {
  x <- df_valid[[bm]]
  g <- df_valid$ager_group
  ok_bm <- !is.na(x) & is.finite(x) & !is.na(g)

  if (sum(ok_bm) < 30) next

  # Per-group statistics
  group_stats <- tapply(x[ok_bm], g[ok_bm], function(v) {
    c(n = length(v), mean = mean(v), sd = sd(v), median = median(v))
  })

  # Kruskal-Wallis test
  kw <- tryCatch(kruskal.test(x[ok_bm] ~ g[ok_bm]), error = function(e) NULL)

  # Spearman correlation with Z-score
  ct <- cor.test(df_valid[[primary_z]][ok_bm], x[ok_bm], method = "spearman")

  # Effect size: mean(P5) - mean(P1) / pooled SD
  p1 <- x[ok_bm & g == "P1_very_slow"]
  p5 <- x[ok_bm & g == "P5_very_fast"]
  if (length(p1) >= 5 && length(p5) >= 5) {
    pooled_sd <- sqrt((var(p1) * (length(p1)-1) + var(p5) * (length(p5)-1)) /
                      (length(p1) + length(p5) - 2))
    cohens_d <- if (pooled_sd > 0) (mean(p5) - mean(p1)) / pooled_sd else NA
  } else {
    cohens_d <- NA
  }

  profile_results[[bm]] <- data.frame(
    biomarker   = bm,
    n_total     = sum(ok_bm),
    rho_z       = round(ct$estimate, 4),
    rho_p       = ct$p.value,
    kw_chi2     = if (!is.null(kw)) round(kw$statistic, 2) else NA,
    kw_p        = if (!is.null(kw)) kw$p.value else NA,
    d_P5_vs_P1  = round(cohens_d, 3),
    mean_P1     = if (length(p1) >= 5) round(mean(p1), 3) else NA,
    mean_P5     = if (length(p5) >= 5) round(mean(p5), 3) else NA,
    direction   = ifelse(!is.na(cohens_d),
                         ifelse(cohens_d > 0.1, "HIGHER_in_fast",
                         ifelse(cohens_d < -0.1, "LOWER_in_fast", "~equal")),
                         NA),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

if (length(profile_results) > 0) {
  profile_table <- do.call(rbind, profile_results)
  profile_table$rho_padj <- p.adjust(profile_table$rho_p, method = "BH")
  profile_table$kw_padj  <- p.adjust(profile_table$kw_p, method = "BH")
  profile_table <- profile_table[order(-abs(profile_table$d_P5_vs_P1)), ]
  rownames(profile_table) <- NULL

  cat("  Top differentiating biomarkers (|Cohen's d| > 0.1, fast vs slow agers):\n\n")
  sig <- profile_table[!is.na(profile_table$d_P5_vs_P1) &
                        abs(profile_table$d_P5_vs_P1) > 0.1, ]
  if (nrow(sig) > 0) {
    print(sig[, c("biomarker", "n_total", "rho_z", "d_P5_vs_P1",
                  "mean_P1", "mean_P5", "direction", "kw_padj")],
          row.names = FALSE)
  }

  write.csv(profile_table,
            file.path(OUTPUT_DIR, "ager_biomarker_profiles.csv"),
            row.names = FALSE)
  message(">> Saved: ager_biomarker_profiles.csv")
}

###############################################################################
#   SECTION 4: FUNCTIONAL MEDICINE REFERENCE RANGES
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 4: FUNCTIONAL MEDICINE REFERENCE RANGES       #\n")
cat("##########################################################\n\n")

cat("  Conventional lab ranges detect disease.\n")
cat("  Functional medicine ranges detect suboptimal health.\n")
cat("  Question: Are functional ranges more sensitive for\n")
cat("  detecting accelerated biological aging?\n\n")

# Functional medicine optimal ranges (published consensus ranges)
# Format: list(biomarker = c(optimal_low, optimal_high, conventional_low, conventional_high))
functional_ranges <- list(
  # Hematology
  hemoglobin        = list(func = c(13.5, 14.5), conv = c(12.0, 17.5), sex_diff = TRUE,
                           func_F = c(12.5, 13.5), func_M = c(14.0, 15.0)),
  hematocrit        = list(func = c(37, 44),    conv = c(36, 52)),
  rbc               = list(func = c(4.2, 4.9),  conv = c(3.8, 5.8)),
  platelets         = list(func = c(200, 300),   conv = c(150, 400)),
  mcv               = list(func = c(85, 92),     conv = c(80, 100)),
  mch               = list(func = c(28, 32),     conv = c(27, 33)),
  mchc              = list(func = c(32, 35),     conv = c(31, 37)),
  rdw               = list(func = c(11.5, 13.0), conv = c(11.5, 14.5)),
  wbc               = list(func = c(5.0, 7.5),   conv = c(4.0, 11.0)),
  # Metabolic
  glucose           = list(func = c(75, 90),     conv = c(70, 110)),
  hba1c             = list(func = c(4.5, 5.3),   conv = c(4.0, 5.7)),
  insulin           = list(func = c(2, 5),       conv = c(2, 25)),
  # Lipids
  cholesterol       = list(func = c(150, 200),   conv = c(0, 200)),
  hdl               = list(func = c(55, 80),     conv = c(40, 100)),
  ldl               = list(func = c(0, 100),     conv = c(0, 130)),
  triglycerides     = list(func = c(50, 100),    conv = c(0, 150)),
  # Iron
  ferritin          = list(func = c(50, 150),    conv = c(10, 300), sex_diff = TRUE,
                           func_F = c(50, 100), func_M = c(50, 150)),
  transferrin_saturation = list(func = c(25, 35), conv = c(15, 50)),
  # Thyroid
  TSH               = list(func = c(1.0, 2.5),   conv = c(0.4, 4.5)),
  free_t4           = list(func = c(1.0, 1.5),   conv = c(0.8, 1.8)),
  # Vitamins
  vitamin_d         = list(func = c(50, 80),     conv = c(30, 100)),
  vitamin_b12       = list(func = c(500, 1000),  conv = c(200, 900)),
  # Inflammatory
  crp               = list(func = c(0, 0.5),     conv = c(0, 5.0)),
  homocysteine      = list(func = c(5, 8),       conv = c(5, 15)),
  # Liver
  ast               = list(func = c(10, 26),     conv = c(0, 40)),
  alt               = list(func = c(10, 26),     conv = c(0, 40)),
  ggt               = list(func = c(10, 30),     conv = c(0, 55)),
  albumin           = list(func = c(4.0, 5.0),   conv = c(3.5, 5.5)),
  # Renal
  creatinine        = list(func = c(0.8, 1.1),   conv = c(0.6, 1.3))
)

# ── For each biomarker: compare conventional vs functional classification ──

func_results <- list()

for (bm in names(functional_ranges)) {
  if (!bm %in% names(df_valid)) next
  x <- df_valid[[bm]]
  z <- df_valid[[primary_z]]
  ok <- !is.na(x) & is.finite(x) & !is.na(z)
  if (sum(ok) < 30) next

  fr <- functional_ranges[[bm]]
  fl <- fr$func[1]; fh <- fr$func[2]
  cl <- fr$conv[1]; ch <- fr$conv[2]

  # Classify: conventional normal, functional suboptimal, conventional abnormal
  conv_abnormal <- x[ok] < cl | x[ok] > ch
  func_suboptimal <- (x[ok] < fl | x[ok] > fh) & !conv_abnormal
  func_optimal <- x[ok] >= fl & x[ok] <= fh

  n_conv_abn  <- sum(conv_abnormal)
  n_func_sub  <- sum(func_suboptimal)
  n_func_opt  <- sum(func_optimal)

  # Mean Z-score in each group
  z_ok <- z[ok]
  z_conv_abn <- if (n_conv_abn >= 5) mean(z_ok[conv_abnormal]) else NA
  z_func_sub <- if (n_func_sub >= 5) mean(z_ok[func_suboptimal]) else NA
  z_func_opt <- if (n_func_opt >= 5) mean(z_ok[func_optimal]) else NA

  # Test: do functionally suboptimal subjects have higher aging Z-scores?
  # Compare functional-optimal vs functional-suboptimal (WITHIN conventional normal)
  conv_normal <- !conv_abnormal
  if (sum(conv_normal & func_suboptimal) >= 10 &&
      sum(conv_normal & func_optimal) >= 10) {
    wt <- wilcox.test(z_ok[conv_normal & func_suboptimal],
                      z_ok[conv_normal & func_optimal])
    # Effect size
    sub_z  <- z_ok[conv_normal & func_suboptimal]
    opt_z  <- z_ok[conv_normal & func_optimal]
    ps <- sqrt((var(sub_z)*(length(sub_z)-1) + var(opt_z)*(length(opt_z)-1)) /
               (length(sub_z) + length(opt_z) - 2))
    d_func <- if (ps > 0) (mean(sub_z) - mean(opt_z)) / ps else NA
  } else {
    wt <- NULL
    d_func <- NA
  }

  func_results[[bm]] <- data.frame(
    biomarker       = bm,
    n               = sum(ok),
    # Conventional classification
    n_conv_normal   = sum(!conv_abnormal),
    n_conv_abnormal = n_conv_abn,
    pct_conv_abn    = round(100 * n_conv_abn / sum(ok), 1),
    z_conv_abnormal = round(z_conv_abn, 3),
    # Functional suboptimal (within conventional normal)
    n_func_subopt   = n_func_sub,
    pct_func_subopt = round(100 * n_func_sub / sum(ok), 1),
    z_func_subopt   = round(z_func_sub, 3),
    # Functional optimal
    n_func_optimal  = n_func_opt,
    z_func_optimal  = round(z_func_opt, 3),
    # Sensitivity comparison
    d_subopt_vs_opt = round(d_func, 3),
    p_subopt_vs_opt = if (!is.null(wt)) wt$p.value else NA,
    # Is functional range more informative than conventional?
    func_sensitive  = ifelse(!is.na(d_func) && abs(d_func) > 0.1 &&
                             !is.null(wt) && wt$p.value < 0.05,
                             "YES", "no"),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

if (length(func_results) > 0) {
  func_table <- do.call(rbind, func_results)
  func_table$p_adj <- p.adjust(func_table$p_subopt_vs_opt, method = "BH")
  func_table <- func_table[order(-abs(func_table$d_subopt_vs_opt)), ]
  rownames(func_table) <- NULL

  cat("--- Functional vs Conventional Range Sensitivity ---\n\n")
  cat("  d_subopt_vs_opt: Cohen's d between functionally-suboptimal and\n")
  cat("  functionally-optimal subjects (WITHIN conventional normal range).\n")
  cat("  Positive d = suboptimal subjects age faster.\n\n")

  print(func_table[, c("biomarker", "n", "pct_conv_abn", "pct_func_subopt",
                        "z_func_optimal", "z_func_subopt",
                        "d_subopt_vs_opt", "func_sensitive")],
        row.names = FALSE)

  sensitive_bm <- func_table$biomarker[func_table$func_sensitive == "YES"]
  if (length(sensitive_bm) > 0) {
    cat(sprintf("\n>> Biomarkers where functional ranges detect aging differences\n"))
    cat(">> that conventional ranges miss (%d found):\n", length(sensitive_bm))
    cat("   ", paste(sensitive_bm, collapse = ", "), "\n")
  }

  write.csv(func_table,
            file.path(OUTPUT_DIR, "functional_range_sensitivity.csv"),
            row.names = FALSE)
  message(">> Saved: functional_range_sensitivity.csv")
}

###############################################################################
#   SECTION 5: MULTI-BIOMARKER AGING SIGNATURE
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 5: AGING SIGNATURE — MULTI-BIOMARKER          #\n")
cat("##########################################################\n\n")

# For the top differentiating biomarkers, build a composite profile
# showing the "accelerated aging phenotype"

if (exists("profile_table") && nrow(profile_table) > 0) {
  # Select top differentiators (|d| > 0.15 or significant KW)
  sig_bm <- profile_table$biomarker[
    !is.na(profile_table$d_P5_vs_P1) &
    abs(profile_table$d_P5_vs_P1) > 0.15 &
    !is.na(profile_table$kw_padj) &
    profile_table$kw_padj < 0.05
  ]

  if (length(sig_bm) >= 3) {
    cat(sprintf(">> Building aging signature from %d significant biomarkers\n\n",
                length(sig_bm)))

    # Radar/spider chart data: mean Z-score per ager group
    radar_data <- list()
    for (bm in sig_bm) {
      if (!bm %in% names(df_valid)) next
      x <- df_valid[[bm]]
      ok_bm <- !is.na(x) & is.finite(x) & !is.na(df_valid$ager_group)
      if (sum(ok_bm) < 30) next

      # Z-score the biomarker within the cohort
      x_z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      group_means <- tapply(x_z[ok_bm], df_valid$ager_group[ok_bm], mean)

      for (g in names(group_means)) {
        radar_data[[length(radar_data) + 1]] <- data.frame(
          biomarker = bm,
          group     = g,
          mean_z    = round(group_means[g], 3),
          stringsAsFactors = FALSE, row.names = NULL
        )
      }
    }

    if (length(radar_data) > 0) {
      radar_table <- do.call(rbind, radar_data)

      # Heatmap-style plot: biomarkers × ager groups
      radar_wide <- pivot_wider(radar_table, names_from = group,
                                values_from = mean_z)
      write.csv(radar_wide,
                file.path(OUTPUT_DIR, "aging_signature_heatmap_data.csv"),
                row.names = FALSE)

      # ggplot heatmap
      radar_table$group <- factor(radar_table$group,
        levels = c("P1_very_slow", "P2_slow", "P3_average",
                   "P4_fast", "P5_very_fast"))
      # Reorder biomarkers by P5-P1 difference
      bm_order <- radar_table %>%
        filter(group %in% c("P1_very_slow", "P5_very_fast")) %>%
        pivot_wider(names_from = group, values_from = mean_z) %>%
        mutate(diff = P5_very_fast - P1_very_slow) %>%
        arrange(diff)
      radar_table$biomarker <- factor(radar_table$biomarker,
                                       levels = bm_order$biomarker)

      p_heat <- ggplot(radar_table, aes(x = group, y = biomarker, fill = mean_z)) +
        geom_tile(colour = "white", linewidth = 0.5) +
        geom_text(aes(label = sprintf("%.2f", mean_z)), size = 3) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                             midpoint = 0, name = "Mean Z") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Accelerated Aging Biomarker Signature",
             subtitle = paste0("Based on ", primary_z, " quintiles"),
             x = "Ager Group", y = "Biomarker")

      ggsave(file.path(OUTPUT_DIR, "aging_signature_heatmap.png"),
             p_heat, width = 10, height = max(6, length(sig_bm) * 0.35 + 2),
             dpi = 150)
      message(">> Saved: aging_signature_heatmap.png")

      # Print the signature
      cat("--- Accelerated Aging Signature ---\n\n")
      print(as.data.frame(radar_wide), row.names = FALSE)
    }
  }
}

###############################################################################
#   SECTION 6: ARGENTINE-ONLY BIOMARKER ANALYSIS
###############################################################################

cat("\n##########################################################\n")
cat("#  SECTION 6: ARGENTINE-ONLY BIOMARKERS                  #\n")
cat("##########################################################\n\n")

cat("  These biomarkers exist in the Argentine dataset but NOT in NHANES.\n")
cat("  They cannot be used for PhenoAge clock training but CAN be used\n")
cat("  for validating and characterizing the aging phenotype.\n\n")

arg_only_available <- intersect(arg_only, names(df_valid))
arg_only_n <- sapply(arg_only_available, function(v) {
  sum(!is.na(df_valid[[v]]) & is.finite(suppressWarnings(as.numeric(df_valid[[v]]))))
})
arg_only_available <- arg_only_available[arg_only_n >= 30]

if (length(arg_only_available) > 0) {
  cat(sprintf("  Available Argentine-only biomarkers: %d\n", length(arg_only_available)))
  cat("  ", paste(arg_only_available, collapse = ", "), "\n\n")

  arg_only_results <- list()
  for (bm in arg_only_available) {
    x <- suppressWarnings(as.numeric(df_valid[[bm]]))
    z <- df_valid[[primary_z]]
    ok <- !is.na(x) & is.finite(x) & !is.na(z)
    if (sum(ok) < 30) next

    ct <- cor.test(z[ok], x[ok], method = "spearman")
    arg_only_results[[bm]] <- data.frame(
      biomarker = bm,
      n         = sum(ok),
      rho       = round(ct$estimate, 4),
      p_value   = ct$p.value,
      stringsAsFactors = FALSE, row.names = NULL
    )
  }

  if (length(arg_only_results) > 0) {
    aor <- do.call(rbind, arg_only_results)
    aor$p_adj <- p.adjust(aor$p_value, method = "BH")
    aor$sig <- ifelse(aor$p_adj < 0.001, "***",
               ifelse(aor$p_adj < 0.01, "**",
               ifelse(aor$p_adj < 0.05, "*", "")))
    aor <- aor[order(-abs(aor$rho)), ]
    rownames(aor) <- NULL

    cat("--- Correlation with Biological Aging (Argentine-only biomarkers) ---\n\n")
    print(aor, row.names = FALSE)
    write.csv(aor, file.path(OUTPUT_DIR, "argentine_only_biomarker_correlations.csv"),
              row.names = FALSE)
  }
}

###############################################################################
#   FINAL SUMMARY
###############################################################################

cat("\n==========================================================\n")
cat("  DISCOVERY ANALYSIS COMPLETE\n")
cat("==========================================================\n")
cat("  Subjects profiled:    ", format(n_valid, big.mark = ","), "\n")
cat("  Primary Z-score:      ", primary_z, "\n")
cat("  Shared biomarkers:    ", nrow(xref), "\n")
cat("  Trainable:            ", n_trainable, "\n")
cat("  Candidate panels:     ", length(candidate_panels), "\n")
cat("  Argentine-only:       ", length(arg_only_available), "\n")
cat("  Functional ranges:    ", length(func_results), " tested\n")
cat("  Output directory:     ", OUTPUT_DIR, "\n")
cat("\n  Files saved:\n")
saved <- list.files(OUTPUT_DIR, full.names = FALSE)
for (f in saved) cat("    - ", f, "\n")
cat("==========================================================\n")
message("\n>> Done. All outputs in: ", OUTPUT_DIR)
