library(BioAge)
library(dplyr)

# ------------------------------------------------------------------
# DEFINÍ TU SET DE BIOMARCADORES (ejemplo laboratorio realista)
# ------------------------------------------------------------------

biomarkers_blood <- c(
  "albumin_gL",     # albúmina
  "alp",            # fosfatasa alcalina
  "lncrp",          # CRP log
  "totchol",        # colesterol total
  "lncreat_umol",   # creatinina log
  "hba1c",          # HbA1c
  "sbp",            # presión sistólica
  "bun",            # urea
  "lymph",          # linfocitos %
  "mcv",            # volumen corpuscular medio
  "wbc"             # glóbulos blancos
)

# ------------------------------------------------------------------
# “ENTRENAMIENTO” PHENOAGE EN NHANES III
# ------------------------------------------------------------------

phenoage_blood <- phenoage_nhanes(
  biomarkers = biomarkers_blood
)

# ------------------------------------------------------------------
# RESULTADOS
# ------------------------------------------------------------------

# Dataset proyectado (por defecto NHANES IV)
data_pheno <- phenoage_blood$data

# Modelo ajustado (ESTO ES LO IMPORTANTE)
fit_pheno <- phenoage_blood$fit

