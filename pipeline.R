###########################################################
# BIOLOGICAL AGING – MULTI-CLOCK FRAMEWORK (FINAL)
# BASED ON BioAge internal NHANES handling
###########################################################

library(BioAge)
library(dplyr)
library(survival)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(corrplot)

# ==========================================================
# 1. RELOJES BASE (NO TOCAR)
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

panel_global <- c(
  "albumin_gL","alp","lncrp","totchol",
  "lncreat_umol","hba1c","sbp","bun",
  "lymph","mcv","wbc"
)

message(">> Entrenando relojes base")

pheno_orig   <- phenoage_nhanes(biomarkers_pheno_orig)
kdm_orig     <- kdm_nhanes(biomarkers_kdm_hd)
hd_orig      <- hd_nhanes(biomarkers_kdm_hd)
pheno_global <- phenoage_nhanes(panel_global)

# ==========================================================
# 2. SUB-RELOJES
# ==========================================================

panels <- list(
  hepatic_enzime_insulin   = c("alp","ggt","insulin"),
  hepatic_lipid            = c("trig","totchol"),
  hema_integrated          = c("rdw","mcv","rbc","wbc","lymph"),
  micronutrient_methylation= c("vitaminB12","hba1c","rdw"),
  renal_A                  = c("lncreat_umol","bun")
)

panel_system <- c(
  hepatic_enzime_insulin   = "Hepatic",
  hepatic_lipid            = "Hepatic",
  hema_integrated          = "Hematologic",
  micronutrient_methylation= "Micronutrients",
  renal_A                  = "Renal"
)

message(">> Entrenando sub-relojes")

trained_clocks <- map(
  panels,
  ~ tryCatch(
    phenoage_nhanes(.x),
    error = function(e) NULL
  )
)

trained_clocks <- trained_clocks[!map_lgl(trained_clocks, is.null)]

# ==========================================================
# 3. DATASET FINAL (ESTE ERA EL CORRECTO)
# ==========================================================

data_all <- pheno_orig$data %>%
  select(sampleID, age, time, status, phenoage) %>%
  rename(phenoage_orig = phenoage) %>%
  
  left_join(
    kdm_orig$data %>% select(sampleID, kdm),
    by = "sampleID"
  ) %>%
  
  left_join(
    hd_orig$data %>% select(sampleID, hd),
    by = "sampleID"
  ) %>%
  
  left_join(
    pheno_global$data %>%
      select(sampleID, phenoage) %>%
      rename(phenoage_global = phenoage),
    by = "sampleID"
  )

for (clock_name in names(trained_clocks)) {
  clock_data <- trained_clocks[[clock_name]]$data %>%
    select(sampleID, phenoage) %>%
    rename(!!paste0("pheno_", clock_name) := phenoage)
  
  data_all <- left_join(data_all, clock_data, by = "sampleID")
}

# ==========================================================
# 4. LISTA DE RELOJES
# ==========================================================

clocks <- names(data_all) %>%
  grep("^pheno_|^kdm$|^hd$", ., value = TRUE)

# ==========================================================
# 5. ESCALADO 1 SD (CRÍTICO)
# ==========================================================

z <- function(x) as.numeric(scale(x))

data_all <- data_all %>%
  mutate(across(all_of(clocks), z, .names = "z_{.col}"))

# ==========================================================
# 6. COX – RELOJES CRUDOS
# ==========================================================

cox_raw <- map_df(
  clocks,
  function(v) {
    df <- data_all %>%
      select(time, status, age, all_of(paste0("z_", v))) %>%
      na.omit()
    
    if (nrow(df) < 200) return(NULL)
    
    f <- as.formula(paste0("Surv(time,status) ~ age + z_", v))
    m <- coxph(f, df)
    s <- summary(m)
    
    tibble(
      clock = v,
      model = "Raw (+1 SD)",
      HR = s$coef[2,"exp(coef)"],
      p  = s$coef[2,"Pr(>|z|)"],
      cindex = s$concordance[1]
    )
  }
)

# ==========================================================
# 7. RESIDUAL AGE (CORRECTO)
# ==========================================================

residualize <- function(y, age) {
  resid(lm(y ~ age, na.action = na.exclude))
}

data_all <- data_all %>%
  mutate(
    across(
      all_of(clocks),
      ~ residualize(.x, data_all$age),
      .names = "resid_{.col}"
    )
  ) %>%
  mutate(
    across(starts_with("resid_"), z, .names = "z_{.col}")
  )

# ==========================================================
# 8. COX – RESIDUAL AGE
# ==========================================================

cox_resid <- map_df(
  clocks,
  function(v) {
    df <- data_all %>%
      select(time, status, age, all_of(paste0("z_resid_", v))) %>%
      na.omit()
    
    if (nrow(df) < 200) return(NULL)
    
    f <- as.formula(paste0("Surv(time,status) ~ age + z_resid_", v))
    m <- coxph(f, df)
    s <- summary(m)
    
    tibble(
      clock = v,
      model = "Residual (+1 SD)",
      HR = s$coef[2,"exp(coef)"],
      p  = s$coef[2,"Pr(>|z|)"],
      cindex = s$concordance[1]
    )
  }
)

# ==========================================================
# 9. RESULTADOS FINALES
# ==========================================================

results <- bind_rows(cox_raw, cox_resid) %>%
  arrange(model, desc(HR))

print(results, n = Inf)

# ==========================================================
# 10. HEATMAP – SUBRELOJES RESIDUALES
# ==========================================================

resid_subs <- grep("^resid_pheno_", names(data_all), value = TRUE)

cor_mat <- data_all %>%
  select(all_of(resid_subs)) %>%
  cor(use = "pairwise.complete.obs")

corrplot(
  cor_mat,
  method = "color",
  tl.cex = 0.9,
  title = "Residual sub-clocks – clinical independence",
  mar = c(0,0,2,0)
)

###########################################################
# FIN – ESTE PIPELINE ES CONSISTENTE Y FUNCIONA
###########################################################


