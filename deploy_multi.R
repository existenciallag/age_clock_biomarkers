############################################################
# BIOLOGICAL AGING – MULTI-CLOCK FRAMEWORK
# NHANES III → NHANES IV (BioAge native)
############################################################

library(BioAge)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

# ======================
# 1. BIOMARCADORES
# ======================

# --- Originales (baseline histórico) ---
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

# --- Nuevo PhenoAge global (sin uap) ---
panel_global <- c(
  "albumin_gL","alp","lncrp","totchol",
  "lncreat_umol","hba1c","sbp","bun",
  "lymph","mcv","wbc"
)

# --- Sub-relojes ---
panel_metabolic <- c(
  "glucose","hba1c","trig","hdl","ldl","totchol","ggt"
)

panel_renal <- c(
  "lncreat_umol","bun"
)

panel_hematopoietic <- c(
  "rbc","mcv","rdw","wbc"
)

panel_exposome <- c(
  "cadmium","ggt","albumin_gL"
)

# ======================
# 2. ENTRENAMIENTO (NHANES)
# ======================

message(">> Entrenando relojes originales")

pheno_orig <- phenoage_nhanes(biomarkers = biomarkers_pheno_orig)
kdm_orig   <- kdm_nhanes(biomarkers = biomarkers_kdm_hd)
hd_orig    <- hd_nhanes(biomarkers = biomarkers_kdm_hd)

message(">> Entrenando PhenoAge global (modificado)")

pheno_global <- phenoage_nhanes(biomarkers = panel_global)

message(">> Entrenando sub-relojes")

ph_metabolic <- phenoage_nhanes(biomarkers = panel_metabolic)
ph_renal     <- phenoage_nhanes(biomarkers = panel_renal)
ph_hema      <- phenoage_nhanes(biomarkers = panel_hematopoietic)
ph_exposome  <- phenoage_nhanes(biomarkers = panel_exposome)

# ======================
# 3. ENSAMBLADO DE DATASET
# ======================

data_all <- pheno_orig$data %>%
  select(sampleID, age, time, status,
         phenoage0, phenoage) %>%
  rename(phenoage_orig = phenoage) %>%
  
  # --- KDM ---
  left_join(
    kdm_orig$data %>%
      select(sampleID, kdm0, kdm),
    by = "sampleID"
  ) %>%
  
  # --- HD ---
  left_join(
    hd_orig$data %>%
      select(sampleID, hd, hd_log),
    by = "sampleID"
  ) %>%
  
  # --- PhenoAge global modificado ---
  left_join(
    pheno_global$data %>%
      select(sampleID, phenoage) %>%
      rename(phenoage_global = phenoage),
    by = "sampleID"
  ) %>%
  
  # --- Sub-relojes ---
  left_join(
    ph_metabolic$data %>%
      select(sampleID, phenoage) %>%
      rename(pheno_metabolic = phenoage),
    by = "sampleID"
  ) %>%
  
  left_join(
    ph_renal$data %>%
      select(sampleID, phenoage) %>%
      rename(pheno_renal = phenoage),
    by = "sampleID"
  ) %>%
  
  left_join(
    ph_hema$data %>%
      select(sampleID, phenoage) %>%
      rename(pheno_hema = phenoage),
    by = "sampleID"
  ) %>%
  
  left_join(
    ph_exposome$data %>%
      select(sampleID, phenoage) %>%
      rename(pheno_exposome = phenoage),
    by = "sampleID"
  )

