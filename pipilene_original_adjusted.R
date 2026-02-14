##########################################################
# BIOLOGICAL AGING – CANONICAL BioAge PIPELINE
# NHANES III → NHANES IV
# Canonical validation + PhenoAge Global + Subclocks trained
##########################################################

library(BioAge)
library(dplyr)
library(purrr)

# ==========================================================
# 1. BIOMARKER SETS
# ==========================================================

biomarkers_pheno_orig <- c(
  "albumin_gL","alp","lncrp","totchol",
  "lncreat_umol","hba1c","sbp","bun",
  "uap","lymph","mcv","wbc"
)

biomarkers_kdm_hd <- c(
  "albumin","alp","lncrp","totchol",
  "lncreat","hba1c","sbp","bun",
  "uap","lymph","mcv","wbc"
)

# GLOBAL = original (sanity check)
panel_global <- biomarkers_pheno_orig


# ==========================================================
# 2. TRAIN CANONICAL CLOCKS (NHANES III)
# ==========================================================

message(">> Training canonical BioAge clocks")

hd      <- hd_nhanes(biomarkers_kdm_hd)
kdm     <- kdm_nhanes(biomarkers_kdm_hd)
pheno   <- phenoage_nhanes(biomarkers_pheno_orig)
pheno_g <- phenoage_nhanes(panel_global)


# ==========================================================
# 3. TRAIN SUBCLOCKS (not evaluated yet)
# ==========================================================

panels_sub <- list(
  hepatic_enzime_insulin    = c("alp","ggt","insulin"),
  hepatic_lipid            = c("trig","totchol"),
  hema_integrated          = c("rdw","mcv","rbc","wbc","lymph"),
  micronutrient_methyl     = c("vitaminB12","hba1c","rdw"),
  renal_A                  = c("lncreat_umol","bun")
)

message(">> Training sub-clocks")

subclocks <- map(
  panels_sub,
  ~ tryCatch(phenoage_nhanes(.x), error = function(e) NULL)
)

subclocks <- subclocks[!map_lgl(subclocks, is.null)]


# ==========================================================
# 4. ASSEMBLE NHANES IV DATA
# ==========================================================

data <- merge(hd$data, kdm$data) %>%
  merge(., pheno$data) %>%
  merge(., pheno_g$data %>%
          select(sampleID, phenoage) %>%
          rename(phenoage_global = phenoage))

# attach subclocks (not used yet)
for (nm in names(subclocks)) {
  tmp <- subclocks[[nm]]$data %>%
    select(sampleID, phenoage) %>%
    rename(!!paste0("pheno_", nm) := phenoage)
  data <- left_join(data, tmp, by = "sampleID")
}


# ==========================================================
# 4b. CREATE ADVANCEMENT FOR GLOBAL PHENOAGE
# ==========================================================

data$phenoage_global_advance <- data$phenoage_global - data$age


# ==========================================================
# 5. BA vs CHRONOLOGICAL AGE
# ==========================================================

agevar_ba <- c(
  "kdm0","phenoage0",
  "kdm","phenoage",
  "phenoage_global",
  "hd","hd_log"
)

label_ba <- c(
  "KDM\nBiological Age",
  "Levine\nPhenotypic Age",
  "Modified-KDM\nBiological Age",
  "Modified-Levine\nPhenotypic Age",
  "Global\nPhenotypic Age",
  "Homeostatic\nDysregulation",
  "Log\nHomeostatic\nDysregulation"
)

plot_ba(data, agevar_ba, label_ba)


# ==========================================================
# 6. BIOLOGICAL AGE ADVANCEMENT (BA − CA)
# ==========================================================

agevar_baa <- c(
  "kdm_advance0","phenoage_advance0",
  "kdm_advance","phenoage_advance",
  "phenoage_global_advance",
  "hd","hd_log"
)

label_baa <- c(
  "kdm_advance0"              = "KDM\nBiological Age\nAdvancement",
  "phenoage_advance0"        = "Levine\nPhenotypic Age\nAdvancement",
  "kdm_advance"              = "Modified-KDM\nBiological Age\nAdvancement",
  "phenoage_advance"         = "Modified-Levine\nPhenotypic Age\nAdvancement",
  "phenoage_global_advance"  = "Global\nPhenotypic Age\nAdvancement",
  "hd"                       = "Homeostatic\nDysregulation",
  "hd_log"                   = "Log\nHomeostatic\nDysregulation"
)

axis_type <- setNames(rep("float", length(agevar_baa)), agevar_baa)

plot_baa(data, agevar_baa, label_baa, axis_type)


# ==========================================================
# 7. MORTALITY & HEALTH
# ==========================================================

table_surv(data, agevar_baa, label_baa)

table_health(
  data,
  agevar  = agevar_baa,
  outcome = c("health","adl","lnwalk","grip_scaled"),
  label   = label_baa
)$table

##########################################################
# END – Canonical BioAge validation + extensions
##########################################################
