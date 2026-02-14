############################################################
# BIOLOGICAL AGING – MULTI-CLOCK FRAMEWORK (EXPLORATORY)
# NHANES III → NHANES IV
############################################################

library(BioAge)
library(dplyr)
library(survival)
library(purrr)
library(tibble)

# ======================
# 1. RELOJES BASE (NO TOCAR)
# ======================

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

panel_global <- c(
  "albumin_gL","alp","lncrp","totchol",
  "lncreat_umol","hba1c","sbp","bun",
  "lymph","mcv","wbc"
)

message(">> Entrenando relojes base")

pheno_orig   <- phenoage_nhanes(biomarkers = biomarkers_pheno_orig)
kdm_orig     <- kdm_nhanes(biomarkers = biomarkers_kdm_hd)
hd_orig      <- hd_nhanes(biomarkers = biomarkers_kdm_hd)
pheno_global <- phenoage_nhanes(biomarkers = panel_global)

# ======================
# 2. DEFINICIÓN DE SUB-RELOJES (AMPLIADOS)
# ======================

panels <- list(
  
  # =====================
  # METABÓLICOS
  # =====================
  metabolic_A = c("glucose","hba1c","insulin","trig","hdl","ldl"),
  metabolic_B = c("ggt","trig","hdl","totchol","albumin_gL"),
  metabolic_C = c("hba1c","trig","hdl"),
  metabolic_integrated = c("hba1c","insulin","trig","hdl","ggt"),
  metabolic_inflamation = c("insulin","lncrp","albumin_gL"),
  metabolic_inflamation = c("insulin","lncrp","albumin_gL","rdw","wbc"),
  
  
  
  # =====================
  # RENALES
  # =====================
  renal_A = c("lncreat_umol","bun"),
  renal_B = c("lncreat_umol","bun","lncrp"),
  renal_inflammatory = c("lncreat_umol","bun","cyst","lncrp"),
  
  # =====================
  # HEMATOLÓGICOS
  # =====================
  hema_A = c("rbc","mcv","rdw"),
  hema_Aa = c("rbc","vitaminB12","rdw"),
  hema_B = c("wbc","lymph"),
  hema_C = c("rdw","wbc","lymph"),
  hema_integrated = c("rdw","mcv","rbc","wbc","lymph"),
  
  # =====================
  # INFLAMACIÓN / PROTEÍNAS
  # =====================
  inflammation_core = c("lncrp","albumin_gL"),
  inflammation_extended = c("lncrp","wbc","rdw","albumin_gL"),
  
  # =====================
  # HEPÁTICO / EXPOSOMA
  # =====================
  hepatic_detox = c("ggt","albumin_gL"),
  hepatic_lipid = c("trig","totchol"),
  exposome_A = c("cadmium","ggt","albumin_gL"),
  exposome_B = c("cadmium","lncrp","wbc"),
  
  # =====================
  # MICRONUTRIENTES
  # =====================
  micronutrient_antioxidant = c(
    "vitaminA","vitaminC","vitaminE"
  ),
  
  micronutrient_methylation = c(
    "vitaminB12","hba1c","rdw"
  ),
  
  micronutrient_integrated = c(
    "vitaminB12","vitaminC","vitaminE","lncrp"
  )
)

# ======================
# 3. ENTRENAMIENTO AUTOMÁTICO (ROBUSTO)
# ======================

message(">> Entrenando sub-relojes")

trained_clocks <- map(
  panels,
  ~ tryCatch(
    phenoage_nhanes(biomarkers = .x),
    error = function(e) {
      message("Panel descartado: ", paste(.x, collapse = ", "))
      NULL
    }
  )
)

# eliminar NULLs
trained_clocks <- trained_clocks[!map_lgl(trained_clocks, is.null)]

# ======================
# 4. DATASET BASE
# ======================

data_all <- pheno_orig$data %>%
  select(sampleID, age, time, status,
         phenoage0, phenoage) %>%
  rename(phenoage_orig = phenoage) %>%
  
  left_join(
    kdm_orig$data %>% select(sampleID, kdm0, kdm),
    by = "sampleID"
  ) %>%
  
  left_join(
    hd_orig$data %>% select(sampleID, hd, hd_log),
    by = "sampleID"
  ) %>%
  
  left_join(
    pheno_global$data %>%
      select(sampleID, phenoage) %>%
      rename(phenoage_global = phenoage),
    by = "sampleID"
  )

# ======================
# 5. AGREGAR SUB-RELOJES
# ======================

for (clock_name in names(trained_clocks)) {
  
  clock_data <- trained_clocks[[clock_name]]$data %>%
    select(sampleID, phenoage) %>%
    rename(!!paste0("pheno_", clock_name) := phenoage)
  
  data_all <- left_join(data_all, clock_data, by = "sampleID")
}

# ======================
# 6. FUNCIÓN DE COX (SEGURA)
# ======================

run_cox <- function(var, data) {
  
  df <- data %>%
    select(time, status, age, phenoage_global, all_of(var)) %>%
    na.omit()
  
  if (nrow(df) < 200) return(NULL)
  
  f1 <- as.formula(paste0("Surv(time, status) ~ ", var, " + age"))
  f2 <- as.formula(paste0("Surv(time, status) ~ ", var, " + age + phenoage_global"))
  
  m1 <- coxph(f1, data = df)
  m2 <- coxph(f2, data = df)
  
  tibble(
    clock = var,
    n = nrow(df),
    HR_age = exp(coef(m1)[1]),
    p_age = summary(m1)$coef[1,5],
    HR_adj_global = exp(coef(m2)[1]),
    p_adj_global = summary(m2)$coef[1,5],
    cindex = summary(m1)$concordance[1]
  )
}

# ======================
# 7. BARRIDO DE HR
# ======================

clock_vars <- names(data_all) %>% grep("^pheno_", ., value = TRUE)

results_hr <- map_dfr(
  clock_vars,
  run_cox,
  data = data_all
) %>%
  arrange(desc(HR_adj_global))

# ======================
# 8. RESULTADOS
# ======================

print(results_hr, n = Inf)


