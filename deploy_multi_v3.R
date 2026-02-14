###########################################################
# BIOLOGICAL AGING – MULTI-CLOCK FRAMEWORK (REDUCED SET)
# NHANES III → NHANES IV
###########################################################

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
# 2. DEFINICIÓN DE SUB-RELOJES (SOLO LOS SELECCIONADOS)
# ======================

panels <- list(
  hepatic_enzime_insulin = c("alp","ggt","insulin"),
  hepatic_lipid = c("trig","totchol"),
  hema_integrated = c("rdw","mcv","rbc","wbc","lymph"),
  micronutrient_methylation = c("vitaminB12","hba1c","rdw"),
  renal_A = c("lncreat_umol","bun")
)

# ======================
# 3. SISTEMA BIOLÓGICO
# ======================

panel_system <- c(
  hepatic_enzime_insulin = "Hepatic",
  hepatic_lipid = "Hepatic",
  hema_integrated = "Hematologic",
  micronutrient_methylation = "Micronutrients",
  renal_A = "Renal"
)

# ======================
# 4. ENTRENAMIENTO AUTOMÁTICO
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

trained_clocks <- trained_clocks[!map_lgl(trained_clocks, is.null)]

# ======================
# 5. DATASET BASE
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
# 6. AGREGAR SUB-RELOJES
# ======================

for (clock_name in names(trained_clocks)) {
  
  clock_data <- trained_clocks[[clock_name]]$data %>%
    select(sampleID, phenoage) %>%
    rename(!!paste0("pheno_", clock_name) := phenoage)
  
  data_all <- left_join(data_all, clock_data, by = "sampleID")
}

# ======================
# 7. FUNCIÓN DE COX
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
# 8. BARRIDO HR
# ======================

clock_vars <- names(data_all) %>% grep("^pheno_", ., value = TRUE)

results_hr <- map_dfr(
  clock_vars,
  run_cox,
  data = data_all
) %>%
  mutate(
    system = panel_system[gsub("^pheno_", "", clock)]
  ) %>%
  arrange(system, desc(HR_adj_global))

# ======================
# 9. RESULTADOS
# ======================

print(results_hr, n = Inf)
