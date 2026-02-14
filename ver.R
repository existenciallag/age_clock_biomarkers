##########################################################
# BIOLOGICAL AGING – FULL OPTIMIZED PIPELINE
# QC: PhenoAge Original (r=0.94) + Subclocks Analysis
##########################################################

# Cargar librerías necesarias
library(BioAge)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(stringr)

# ==========================================================
# 1. DEFINICIÓN DE PANELES (Basado en Paper y Objetivos)
# ==========================================================

# Panel para KDM y HD (Variables estándar en BioAge)
biomarkers_kdm_hd <- c("albumin","alp","lncrp","totchol","lncreat","hba1c","sbp","bun","uap","lymph","mcv","wbc")

# Panel PhenoAge Original (Referencia para QC r=0.94)
biomarkers_pheno_orig <- c("albumin_gL","alp","lncrp","totchol","lncreat_umol","hba1c","sbp","bun","uap","lymph","mcv","wbc")

# Sub-relojes para interpretación bioquímica específica
panels_sub <- list(
  hepatic_enzime  = c("alp","ggt","insulin"),
  hepatic_lipid   = c("trig","totchol"),
  hematologic     = c("rdw","mcv","rbc","wbc","lymph"),
  renal           = c("lncreat_umol","bun")
)

# ==========================================================
# 2. ENTRENAMIENTO DE RELOJES (NHANES III -> NHANES IV)
# ==========================================================
message(">> Entrenando relojes canónicos y sub-relojes...")

hd    <- hd_nhanes(biomarkers_kdm_hd)
kdm   <- kdm_nhanes(biomarkers_kdm_hd)
pheno <- phenoage_nhanes(biomarkers_pheno_orig)

# Entrenar sub-relojes dinámicamente
subclocks <- map(panels_sub, ~ tryCatch(phenoage_nhanes(.x), error = function(e) NULL))
subclocks <- subclocks[!map_lgl(subclocks, is.null)]

# ==========================================================
# 3. ENSAMBLAJE DE DATOS CON PROTECCIÓN DE ERRORES
# ==========================================================

# 3a. Identificar columnas de mortalidad (pueden variar según la versión del dataset)
all_cols <- names(pheno$data)
mort_col <- all_cols[grep("mortstat|dies|dead|outcome", all_cols, ignore.case = TRUE)][1]
time_col <- all_cols[grep("permth|statsurv|followup|surv", all_cols, ignore.case = TRUE)][1]

# Si no se encuentran, crear placeholders para evitar errores en select()
if(is.na(mort_col)) {
  message("⚠️ Advertencia: No se detectó columna de mortalidad (mortstat/dies).")
  pheno$data$mort_placeholder <- NA
  pheno$data$time_placeholder <- NA
  mort_col <- "mort_placeholder"
  time_col <- "time_placeholder"
}

# 3b. Crear dataframe final uniendo todo por sampleID
data_final <- pheno$data %>%
  select(sampleID, age, gender, 
         mortality = all_of(mort_col), 
         survival_time = all_of(time_col), 
         phenoage_orig = phenoage, 
         phenoage_adv_orig = phenoage_advance) %>%
  left_join(hd$data %>% select(sampleID, hd, hd_log), by = "sampleID") %>%
  left_join(kdm$data %>% select(sampleID, kdm, kdm_advance), by = "sampleID")

# Unir sub-relojes
for (nm in names(subclocks)) {
  tmp <- subclocks[[nm]]$data %>%
    select(sampleID, phenoage) %>%
    rename(!!paste0("pheno_", nm) := phenoage)
  data_final <- left_join(data_final, tmp, by = "sampleID")
}

# ==========================================================
# 4. CONTROL DE CALIDAD (QC)
# ==========================================================
r_val <- cor(data_final$age, data_final$phenoage_orig, use = "complete.obs")

cat("\n--- RESULTADOS DE CONTROL DE CALIDAD ---\n")
cat("Correlación PhenoAge vs Edad Cronológica (Meta: ~0.94): ", round(r_val, 4), "\n")
cat("----------------------------------------\n\n")

# ==========================================================
# 5. ANÁLISIS DE SUPERVIVENCIA (HAZARD RATIOS)
# ==========================================================

# 5a. Estandarización (Z-Scores de los avances)
# Esto permite comparar la potencia de los relojes en la misma escala
data_final <- data_final %>%
  mutate(across(starts_with("pheno_"), ~ as.numeric(scale(. - age)), .names = "{.col}_z")) %>%
  mutate(pheno_orig_z = as.numeric(scale(phenoage_adv_orig)),
         kdm_z = as.numeric(scale(kdm_advance)),
         hd_z = as.numeric(scale(hd_log)))

# 5b. Ejecutar tabla de mortalidad (solo si hay datos de mortalidad reales)
vars_z <- c("pheno_orig_z", "kdm_z", "hd_z", grep("pheno_.*_z$", names(data_final), value = TRUE))
labels_z <- gsub("pheno_|_z", "", vars_z)

# Solo ejecutamos si la columna de mortalidad no es puras NA
if(!all(is.na(data_final$mortality))) {
  surv_table <- table_surv(data_final, vars_z, labels_z)
  print(surv_table)
  
  # Preparar datos para Forest Plot
  hr_plot_data <- surv_table %>%
    mutate(
      HR = as.numeric(str_extract(Mortality, "^[0-9.]*")),
      Lower = as.numeric(str_extract(Mortality, "(?<=\\()[0-9.]*")),
      Upper = as.numeric(str_extract(Mortality, "(?<=, )[0-9.]*"))
    )
  
  # Visualización de Hazard Ratios
  ggplot(hr_plot_data, aes(x = reorder(Label, HR), y = HR)) +
    geom_point(size = 4, color = "#2c3e50") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "#2c3e50") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
    coord_flip() +
    labs(title = "Hazard Ratios: Reloj Global vs Órgano-Específicos",
         subtitle = "Riesgo de mortalidad por 1 Desviación Estándar de avance",
         x = "Reloj Biológico", y = "Hazard Ratio (95% CI)") +
    theme_minimal()
} else {
  message(">> Saltando Forest Plot: No hay datos de mortalidad vinculados.")
}

# ==========================================================
# 6. VISUALIZACIÓN DE CORRELACIÓN (QC Gráfico)
# ==========================================================
plot_ba(data_final, "phenoage_orig", "PhenoAge Original (QC)")

library(survival)
library(corrplot)

# ==========================================================
# 6. QC / DIAGNOSTIC TABLE
# ==========================================================
# Primero aseguramos que todos los avances de subrelojes existan en 'data'
for (nm in names(subclocks)) {
  varname <- paste0("pheno_", nm)
  advname <- paste0(varname, "_advance")
  if(varname %in% names(data)) {
    data[[advname]] <- data[[varname]] - data$age
  }
}

all_clocks <- c(
  "kdm_advance", "phenoage_adv_orig", "hd_log",
  grep("pheno_.*_advance$", names(data), value = TRUE)
)

qc_table <- lapply(all_clocks, function(v){
  if(!v %in% names(data)) return(NULL)
  x <- data[[v]]
  data.frame(
    variable   = v,
    n_total    = nrow(data),
    n_nonNA    = sum(!is.na(x)),
    n_finite   = sum(is.finite(x), na.rm=TRUE),
    sd_raw     = sd(x, na.rm=TRUE),
    frac_valid = mean(is.finite(x), na.rm=TRUE)
  )
}) %>% bind_rows()

message(">> QC Table:")
print(qc_table[order(qc_table$frac_valid),], row.names=FALSE)

# Filtramos relojes con suficiente n para que el HR sea estable
valid_clocks <- qc_table %>%
  filter(n_finite > 2000, sd_raw > 0) %>%
  pull(variable)

# ==========================================================
# 7. RESIDUALIZATION & STANDARDIZATION (Z-score: SD = 1)
# ==========================================================
# Extraemos el residuo de la regresión contra edad cronológica
# y luego dividimos por su SD para que el aumento sea por "1 SD"
for(v in valid_clocks){
  f <- as.formula(paste(v, "~ age"))
  # Residualizar
  res <- resid(lm(f, data = data, na.action = na.exclude))
  # Estandarizar el residuo (Z-score: mean 0, sd 1)
  data[[paste0(v, "_resid_z")]] <- as.numeric(scale(res))
}

# ==========================================================
# 8. HAZARD RATIO / SURVIVAL ANALYSIS (REVISADO)
# ==========================================================
mort_col <- intersect(names(data), c("mortstat", "dies", "dead"))[1]
time_col <- intersect(names(data), c("permth_int", "statsurv", "followup"))[1]

if(!is.na(mort_col) & !is.na(time_col)){
  
  # Usamos map_dfr para asegurar que el resultado sea un dataframe limpio
  surv_results <- map_dfr(valid_clocks, function(v){
    dv <- paste0(v, "_resid_z")
    
    # Seleccionamos y limpiamos
    sub_data <- data %>% 
      select(all_of(c(mort_col, time_col, dv, "age", "gender"))) %>% 
      filter(is.finite(get(dv))) %>%
      na.omit()
    
    if(nrow(sub_data) < 2000) return(NULL)
    
    # Cox Model
    formula_cox <- as.formula(paste0("Surv(", time_col, ",", mort_col, ") ~ ", dv, " + age + gender"))
    cox <- coxph(formula_cox, data = sub_data)
    
    # Extraer métricas con nombres explícitos para evitar el error de ggplot
    data.frame(
      clock  = gsub("pheno_|age_adv_orig|_advance|_adv", "", v),
      N      = nrow(sub_data),
      HR     = exp(coef(cox)[[1]]),
      LCL    = exp(confint(cox)[1,1]),
      UCL    = exp(confint(cox)[1,2]),
      Cindex = summary(cox)$concordance[1]
    )
  })
  
  # Verificación de la tabla
  print("Tabla de Supervivencia generada:")
  print(surv_results)
  
} else {
  message(">> No se encontraron columnas de mortalidad.")
}

# ==========================================================
# 10. VISUALIZACIÓN MEJORADA
# ==========================================================
if(exists("surv_results") && nrow(surv_results) > 0) {
  
  # Forest Plot de Hazard Ratios
  p1 <- ggplot(surv_results, aes(x = reorder(clock, HR), y = HR)) +
    geom_point(size = 4, color = "steelblue") +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2, color = "steelblue") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = "Hazard Ratio por Sistema",
         subtitle = "Efecto de 1 SD de Edad Residual (Ajustado)",
         x = "Sistema Fisiológico", 
         y = "Hazard Ratio (Mortalidad)") +
    theme_minimal()
  
  print(p1)
  
  
  
  # Gráfico de Correlación de Residuos
  # Aseguramos que los nombres coincidan con los del HR para comparar
  resid_vars <- grep("_resid_z$", names(data), value = TRUE)
  cor_resid <- cor(data[, resid_vars], use = "pairwise.complete.obs")
  colnames(cor_resid) <- rownames(cor_resid) <- gsub("pheno_|age_adv_orig|_advance|_adv|_resid_z", "", colnames(cor_resid))
  
  corrplot(cor_resid, method="color", type="upper", 
           addCoef.col="black", tl.col="black", 
           title="\nIndependencia de Información (Residuos)",
           mar=c(0,0,1,0))
  
  
}