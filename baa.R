############################################################
# RESIDUAL BIOLOGICAL AGE PIPELINE
# Clinical & Mortality-Oriented Framework
# NO BAA – ONLY RESIDUAL AGE
############################################################

# ======================
# 0. LIBRERÍAS
# ======================

library(dplyr)
library(survival)
library(ggplot2)
library(tidyr)
library(forcats)
library(purrr)

# ======================
# 1. CHEQUEO DE VARIABLES
# ======================

required_vars <- c(
  "age","time","status",
  "phenoage_global",
  "pheno_metabolic",
  "pheno_renal",
  "pheno_hema",
  "pheno_exposome",
  "kdm"
)

stopifnot(all(required_vars %in% names(data_all)))

# ======================
# 2. FUNCIÓN RESIDUAL AGE
# ======================
# Residual = desviación biológica pura
# na.exclude mantiene longitud del dataset

calc_residual <- function(var, data) {
  f <- as.formula(paste0(var, " ~ age"))
  m <- lm(f, data = data, na.action = na.exclude)
  resid(m)
}

# ======================
# 3. CALCULAR RESIDUAL AGE
# ======================

data_resid <- data_all %>%
  mutate(
    resid_pheno_global   = calc_residual("phenoage_global", .),
    resid_pheno_metab    = calc_residual("pheno_metabolic", .),
    resid_pheno_renal    = calc_residual("pheno_renal", .),
    resid_pheno_hema     = calc_residual("pheno_hema", .),
    resid_pheno_exposome = calc_residual("pheno_exposome", .),
    resid_kdm            = calc_residual("kdm", .)
  )

# ======================
# 4. SANITY CHECK
# ======================

stopifnot(nrow(data_resid) == length(data_resid$resid_pheno_global))

resid_vars <- grep("^resid_", names(data_resid), value = TRUE)

# ======================
# 5. MODELOS DE COX
# ======================
# Ajustados por edad cronológica

cox_models <- map(
  resid_vars,
  ~ coxph(
    as.formula(paste0("Surv(time, status) ~ age + ", .x)),
    data = data_resid
  )
)

names(cox_models) <- resid_vars

# ======================
# 6. TABLA DE HAZARD RATIOS
# ======================

hr_table <- map_df(
  cox_models,
  function(m) {
    s <- summary(m)
    data.frame(
      HR    = s$coefficients[2,"exp(coef)"],
      lower = s$conf.int[2,"lower .95"],
      upper = s$conf.int[2,"upper .95"],
      p     = s$coefficients[2,"Pr(>|z|)"]
    )
  },
  .id = "clock"
) %>%
  mutate(clock = fct_reorder(clock, HR)) %>%
  arrange(desc(HR))

print(hr_table)

# ======================
# 7. FOREST PLOT
# ======================

ggplot(
  hr_table,
  aes(x = clock, y = HR, ymin = lower, ymax = upper)
) +
  geom_pointrange(size = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Residual Biological Age and Mortality",
    subtitle = "Hazard Ratio per +1 residual year (adjusted for age)",
    x = "",
    y = "Hazard Ratio"
  )

# ======================
# 8. DISTRIBUCIÓN DE RESIDUAL AGE
# ======================

data_resid %>%
  select(all_of(resid_vars)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "clock",
    values_to = "residual_age"
  ) %>%
  ggplot(aes(x = residual_age)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~ clock, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of Residual Biological Age",
    x = "Residual age (years)",
    y = "Density"
  )

# ======================
# 9. ESTRATIFICACIÓN CLÍNICA
# ======================
# Terciles de riesgo

data_resid <- data_resid %>%
  mutate(
    resid_global_q = ntile(resid_pheno_global, 3)
  )

data_resid$resid_global_q <- factor(
  data_resid$resid_global_q,
  levels = c(1,2,3),
  labels = c(
    "Low biological risk",
    "Intermediate risk",
    "High biological risk"
  )
)

# ======================
# 10. KAPLAN–MEIER
# ======================

km_fit <- survfit(
  Surv(time, status) ~ resid_global_q,
  data = data_resid
)

plot(
  km_fit,
  col = c("darkgreen","gray40","darkred"),
  lwd = 2,
  xlab = "Follow-up time",
  ylab = "Survival probability",
  main = "Survival by Residual Biological Age"
)

legend(
  "bottomleft",
  legend = levels(data_resid$resid_global_q),
  col = c("darkgreen","gray40","darkred"),
  lwd = 2,
  bty = "n"
)

# ======================
# 11. COX CATEGÓRICO
# ======================

cox_cat <- coxph(
  Surv(time, status) ~ age + resid_global_q,
  data = data_resid
)

summary(cox_cat)

############################################################
# FIN
############################################################
