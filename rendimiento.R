############################################################
# BIOLOGICAL AGING – COMPARATIVE ANALYSIS (REDUCED SET)
# Correlations, Heatmap, Global vs Subsystems, Cox
############################################################

library(dplyr)
library(survival)
library(purrr)
library(tibble)
library(ggplot2)
library(reshape2)

# ==========================================================
# 1. DEFINICIÓN AUTOMÁTICA DE RELOJES
# ==========================================================

base_clocks <- c(
  "phenoage_orig",
  "phenoage_global",
  "kdm",
  "hd"
)

sub_clocks <- names(data_all) %>%
  grep("^pheno_", ., value = TRUE) %>%
  setdiff("phenoage_orig")

clock_vars <- c(base_clocks, sub_clocks)
clock_vars <- clock_vars[clock_vars %in% names(data_all)]

# ==========================================================
# 2. CORRELACIÓN CON EDAD CRONOLÓGICA
# ==========================================================

cor_age <- map_dfr(clock_vars, function(v) {
  tibble(
    clock = v,
    cor_age = cor(data_all[[v]], data_all$age, use = "complete.obs")
  )
})

# ==========================================================
# 3. CORRELACIÓN VS PHENOAGE ORIGINAL
# ==========================================================

cor_pheno_orig <- map_dfr(
  setdiff(clock_vars, "phenoage_orig"),
  function(v) {
    tibble(
      clock = v,
      cor_pheno_orig = cor(
        data_all[[v]],
        data_all$phenoage_orig,
        use = "complete.obs"
      )
    )
  }
)

# ==========================================================
# 4. GLOBAL VS SUB-RELOJES
# ==========================================================

cor_global_subs <- map_dfr(
  sub_clocks,
  function(v) {
    tibble(
      subsystem = v,
      cor_global = cor(
        data_all$phenoage_global,
        data_all[[v]],
        use = "complete.obs"
      )
    )
  }
)

# ==========================================================
# 5. COX – AJUSTADO POR EDAD
# ==========================================================

run_cox_simple <- function(var, data) {
  
  df <- data %>%
    select(time, status, age, all_of(var)) %>%
    na.omit()
  
  if (nrow(df) < 200) return(NULL)
  
  fit <- coxph(
    as.formula(paste0("Surv(time, status) ~ age + ", var)),
    data = df
  )
  
  s <- summary(fit)
  
  tibble(
    clock = var,
    n = nrow(df),
    HR = s$coefficients[2, "exp(coef)"],
    lower = s$conf.int[2, "lower .95"],
    upper = s$conf.int[2, "upper .95"],
    p = s$coefficients[2, "Pr(>|z|)"],
    cindex = s$concordance[1]
  )
}

hr_table <- map_dfr(clock_vars, run_cox_simple, data = data_all)

# ==========================================================
# 6. TABLA RESUMEN INTEGRADA
# ==========================================================

summary_table <- cor_age %>%
  left_join(cor_pheno_orig, by = "clock") %>%
  left_join(hr_table, by = "clock") %>%
  arrange(desc(HR))

print(summary_table, n = Inf)

# ==========================================================
# 7. HEATMAP DE CORRELACIONES ENTRE RELOJES
# ==========================================================

cor_matrix <- cor(
  data_all[, clock_vars],
  use = "pairwise.complete.obs"
)

cor_long <- melt(cor_matrix)

ggplot(cor_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(
    title = "Correlation heatmap – Biological aging clocks",
    subtitle = "Red = shared biology | Blue = orthogonal systems",
    fill = "r"
  )

# ==========================================================
# 8. DISCORDANCIA ENTRE SISTEMAS (ADAPTADO)
# ==========================================================

# Z-scores por reloj
z_clocks <- data_all %>%
  mutate(across(all_of(sub_clocks), scale)) %>%
  select(all_of(sub_clocks))

# Definición explícita de sistemas
discordance_table <- z_clocks %>%
  mutate(
    hepatic = rowMeans(
      select(., pheno_hepatic_enzime_insulin, pheno_hepatic_lipid),
      na.rm = TRUE
    ),
    renal = z_clocks$pheno_renal_A,
    hematologic = z_clocks$pheno_hema_integrated,
    micronutrients = z_clocks$pheno_micronutrient_methylation,
    
    hepatic_vs_renal = hepatic - renal,
    hepatic_vs_hema = hepatic - hematologic,
    hepatic_vs_micro = hepatic - micronutrients
  ) %>%
  select(starts_with("hepatic"))

summary(discordance_table)
