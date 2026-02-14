##########################################################
# BIOLOGICAL AGING – FULL OPTIMIZED PIPELINE (CORRECTED)
# QC: PhenoAge Original (r=0.94) + Subclocks Analysis
##########################################################

library(BioAge)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(stringr)
library(survival)

# ==========================================================
# 1. DEFINICIÓN DE PANELES
# ==========================================================
biomarkers_kdm_hd <- c("albumin","alp","lncrp","totchol","lncreat","hba1c","sbp","bun","uap","lymph","mcv","wbc")
biomarkers_pheno_orig <- c("albumin_gL","alp","lncrp","totchol","lncreat_umol","hba1c","sbp","bun","uap","lymph","mcv","wbc")

panels_sub <- list(
  hepatic_enzime  = c("alp","ggt","insulin"),
  hepatic_lipid   = c("trig","totchol"),
  hematologic     = c("rdw","mcv","rbc","wbc","lymph"),
  renal           = c("lncreat_umol","bun")
)

# ==========================================================
# 2. ENTRENAMIENTO DE RELOJES
# ==========================================================
message(">> Entrenando relojes...")
hd    <- hd_nhanes(biomarkers_kdm_hd)
kdm   <- kdm_nhanes(biomarkers_kdm_hd)
pheno <- phenoage_nhanes(biomarkers_pheno_orig)

subclocks <- map(panels_sub, ~ tryCatch(phenoage_nhanes(.x), error = function(e) NULL))
subclocks <- subclocks[!map_lgl(subclocks, is.null)]

# ==========================================================
# 3. ENSAMBLAJE DE DATOS (CORREGIDO PARA TIME/STATUS)
# ==========================================================
# Extraemos la data base de PhenoAge vinculando 'status' y 'time' correctamente
data_final <- pheno$data %>%
  select(sampleID, age, gender, 
         mortality = status,       # Mapeo clave: status es la mortalidad
         survival_time = time,     # Mapeo clave: time es el seguimiento
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
cat("Correlación PhenoAge vs Edad Cronológica: ", round(r_val, 4), "\n")
cat("----------------------------------------\n\n")

# ==========================================================
# 7. CÁLCULO DE EDAD RESIDUAL (Ω) Y ESTANDARIZACIÓN
# ==========================================================
message(">> Calculando residuos Ω estandarizados (Z-scores)...")

# Identificamos las columnas de edad biológica (original y sub-relojes)
ba_vars <- c("phenoage_orig", grep("^pheno_", names(data_final), value = TRUE))
ba_vars <- unique(ba_vars)

for(v in ba_vars) {
  if(v %in% names(data_final)){
    f_res <- as.formula(paste(v, "~ age"))
    # Calculamos residuo puro (independiente de la edad cronológica)
    residuos_brutos <- resid(lm(f_res, data = data_final, na.action = na.exclude))
    # Z-score para que todos los HR sean comparables (Riesgo por 1 SD de exceso de edad)
    data_final[[paste0(v, "_resid_z")]] <- as.numeric(scale(residuos_brutos))
  }
}

# ==========================================================
# 8. ANÁLISIS DE SUPERVIVENCIA (HAZARD RATIO DE RESIDUOS)
# ==========================================================
resid_vars <- grep("_resid_z$", names(data_final), value = TRUE)

# Verificamos que existan datos de mortalidad antes de correr
if(!all(is.na(data_final$mortality))) {
  
  results_hr <- map_dfr(resid_vars, function(v) {
    # Subset independiente por cada reloj para maximizar N
    df_temp <- data_final %>%
      select(mortality, survival_time, age, gender, all_of(v)) %>%
      na.omit()
    
    if(nrow(df_temp) < 200) return(NULL) 
    
    # Modelo de Cox: ¿Cuánto riesgo aporta cada SD de exceso biológico?
    formula_cox <- as.formula(paste0("Surv(survival_time, mortality) ~ ", v, " + age + gender"))
    fit <- survival::coxph(formula_cox, data = df_temp)
    
    res <- summary(fit)
    data.frame(
      Reloj = gsub("pheno_|_resid_z", "", v),
      N = res$n,
      Events = res$nevent,
      HR = res$coefficients[1, "exp(coef)"],
      LCL = exp(confint(fit)[1, 1]),
      UCL = exp(confint(fit)[1, 2]),
      p_val = res$coefficients[1, "Pr(>|z|)"],
      Cindex = res$concordance[1]
    )
  }) %>% arrange(desc(HR))
  
  print("--- TABLA FINAL DE HAZARD RATIOS (RESIDUOS Z) ---")
  print(results_hr)
  
  # ==========================================================
  # 9. VISUALIZACIÓN (FOREST PLOT)
  # ==========================================================
  
  
  ggplot(results_hr, aes(x = reorder(Reloj, HR), y = HR)) +
    geom_point(aes(size = N), color = "#b2182b") + 
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2, color = "#2166ac") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    coord_flip() +
    labs(title = "Impacto de la Edad Residual en la Mortalidad",
         subtitle = "Hazard Ratio por cada 1 SD de exceso biológico (Ajustado por Edad/Sexo)",
         x = "Sistema / Órgano", 
         y = "Hazard Ratio (95% CI)",
         size = "N de la muestra") +
    theme_minimal()
  
} else {
  message(">> Error: No se pudo vincular la mortalidad. Revisa los nombres de las columnas.")
}

# Gráfico de control PhenoAge Original
plot_ba(data_final, "phenoage_orig", "PhenoAge Original (QC)")

# ==========================================================
# 11. ANÁLISIS DE INDEPENDENCIA BIOLÓGICA (VALIDACIÓN ESTRUCTURAL)
# ==========================================================
library(corrplot)
library(tidyr)

message(">> Iniciando validación de independencia de sistemas...")

# 1. Identificar variables de residuos
resid_cols <- grep("_resid_z$", names(data_final), value = TRUE)

# ----------------------------------------------------------
# PRUEBA A: Análisis de Solapamiento (Overlap)
# ----------------------------------------------------------
# Esto nos dice cuántas personas tienen TODOS los relojes a la vez
data_overlap <- data_final %>% select(all_of(resid_cols))
n_total <- nrow(data_overlap)
n_completo <- nrow(na.omit(data_overlap))
pct_overlap <- round((n_completo / n_total) * 100, 1)

cat("\n--- DIAGNÓSTICO DE MUESTRA ---")
cat("\nSujetos totales en dataset:", n_total)
cat("\nSujetos con TODOS los sistemas calculados:", n_completo, "(", pct_overlap, "% )")
cat("\n------------------------------\n")

# ----------------------------------------------------------
# PRUEBA B: Correlación en Población Idéntica (Listwise)
# ----------------------------------------------------------
# Para asegurar que la independencia no es por poblaciones distintas,
# filtramos solo a los que tienen todos los datos.
data_strict <- data_final %>% 
  select(all_of(resid_cols)) %>% 
  na.omit()

# Limpiar nombres para visualización
colnames(data_strict) <- gsub("pheno_|_resid_z", "", colnames(data_strict))
colnames(data_strict) <- str_to_title(colnames(data_strict))

# Calcular matriz
cor_matrix_strict <- cor(data_strict, method = "pearson")

# ----------------------------------------------------------
# PRUEBA C: Visualización de Independencia (Clustering)
# ----------------------------------------------------------
# Si n_completo es muy bajo, avisamos, pero igual graficamos 
# la matriz de los sujetos compartidos.

if(n_completo < 100) {
  warning("⚠️ El solapamiento de sujetos es muy bajo para conclusiones firmes.")
}

cat("\n>> Generando Matriz de Correlación (Sujetos Compartidos)...\n")

corrplot(cor_matrix_strict, 
         method = "color", 
         type = "upper", 
         order = "hclust",      # Agrupa por similitud biológica
         addCoef.col = "black",  # Coeficientes visibles
         tl.col = "black", 
         tl.srt = 45, 
         diag = FALSE,
         title = paste0("\nIndependencia Biológica (N = ", n_completo, " sujetos compartidos)"),
         mar = c(0,0,3,0))

# ----------------------------------------------------------
# EXTRA: Comparación de Edad entre Grupos de Relojes
# ----------------------------------------------------------
# Verificamos si un reloj atrae a gente más vieja que otro
resumen_edad <- map_dfr(resid_cols, function(col) {
  reloj_nom <- gsub("pheno_|_resid_z", "", col)
  data_final %>% 
    filter(!is.na(!!sym(col))) %>% 
    summarise(Sistema = reloj_nom, 
              N = n(),
              Edad_Media = mean(age, na.rm=TRUE),
              Edad_SD = sd(age, na.rm=TRUE))
})

cat("\n--- COMPARATIVA DEMOGRÁFICA POR SISTEMA ---")
print(resumen_edad)
cat("-------------------------------------------\n")