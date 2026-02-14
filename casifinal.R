##########################################################
# BIOLOGICAL AGING – FINAL CLEAN BioAge Pipeline
# NHANES III → NHANES IV
# Canonical + Subclocks + QC + Residuals + HR + Correlation
##########################################################

# =========================
# LIBRERÍAS
# =========================
library(BioAge)
library(dplyr)
library(purrr)
library(survival)
if(!require(corrplot)) install.packages("corrplot")
library(corrplot)

# =========================
# 1. BIOMARKER SETS
# =========================
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

panel_global <- biomarkers_pheno_orig

# =========================
# 2. TRAIN CANONICAL CLOCKS
# =========================
message(">> Training canonical BioAge clocks")
hd      <- hd_nhanes(biomarkers_kdm_hd)
kdm     <- kdm_nhanes(biomarkers_kdm_hd)
pheno   <- phenoage_nhanes(biomarkers_pheno_orig)
pheno_g <- phenoage_nhanes(panel_global)

# =========================
# 3. TRAIN SUBCLOCKS
# =========================
panels_sub <- list(
  hepatic_enzime_insulin = c("alp","ggt","insulin"),
  hepatic_lipid          = c("trig","totchol"),
  hema_integrated        = c("rdw","mcv","rbc","wbc","lymph"),
  micronutrient_methyl   = c("vitaminB12","hba1c","rdw"),
  renal_A                = c("lncreat_umol","bun")
)

message(">> Training sub-clocks")
subclocks <- map(
  panels_sub,
  ~ tryCatch(phenoage_nhanes(.x), error=function(e) NULL)
)
subclocks <- subclocks[!map_lgl(subclocks, is.null)]

# =========================
# 4. ASSEMBLE NHANES IV
# =========================
message(">> Assembling NHANES IV data")
data <- hd$data %>%
  left_join(kdm$data, by="sampleID", suffix=c("_hd","_kdm")) %>%
  left_join(pheno$data, by="sampleID") %>%
  left_join(pheno_g$data %>%
              select(sampleID, phenoage) %>%
              rename(phenoage_global = phenoage),
            by="sampleID")

# attach subclocks
for(nm in names(subclocks)){
  tmp <- subclocks[[nm]]$data %>%
    select(sampleID, phenoage) %>%
    rename(!!paste0("pheno_", nm) := phenoage)
  data <- left_join(data, tmp, by="sampleID")
}

# =========================
# 5. COMPUTE ADVANCEMENT (BA − CA)
# =========================
data <- data %>%
  mutate(
    phenoage_advance        = phenoage - age,
    phenoage_global_advance = phenoage_global - age,
    kdm_advance             = kdm - age,
    hd_log                  = log(hd)
  )

for(nm in names(subclocks)){
  v <- paste0("pheno_",nm)
  data[[paste0(v,"_advance")]] <- data[[v]] - data$age
}

# =========================
# 5b. CORRELATION WITH AGE & BETWEEN BA & BAA
# =========================

# BA vs Chronological Age
ba_vars <- c("kdm","phenoage","phenoage_global", paste0("pheno_", names(subclocks)))
ba_vars <- intersect(ba_vars, names(data))
cor_ba_age <- cor(data[, ba_vars], data$age, use="pairwise.complete.obs")
print("Correlación entre BA y edad:")
print(round(cor_ba_age,2))

# BA vs BAA
baa_vars <- c("kdm_advance","phenoage_advance","phenoage_global_advance",
              paste0("pheno_", names(subclocks), "_advance"))
baa_vars <- intersect(baa_vars, names(data))
cor_ba_baa <- cor(data[, ba_vars], data[, baa_vars], use="pairwise.complete.obs")
print("Correlación entre BA y BAA:")
print(round(cor_ba_baa,2))


# ==========================================================
# 6b. PLOT: BA vs BAA
# ==========================================================

ba_vars  <- agevar_ba[1:length(agevar_ba)]  # tus BA
baa_vars <- agevar_baa[1:length(agevar_baa)] # tus BAA correspondientes

# iterar y graficar
par(mfrow=c(ceiling(length(ba_vars)/2), 2), mar=c(4,4,2,1))  # layout limpio

for(i in seq_along(ba_vars)){
  ba  <- ba_vars[i]
  baa <- baa_vars[i]
  
  if(!(ba %in% names(data)) | !(baa %in% names(data))) next
  
  plot(data[[ba]], data[[baa]],
       xlab = label_ba[i],
       ylab = label_baa[baa],
       pch  = 19,
       col  = rgb(0.2,0.4,0.6,0.5),
       main = paste(label_ba[i], "vs", label_baa[baa]))
  
  abline(lm(data[[baa]] ~ data[[ba]]), col="red", lwd=2)  # regresión lineal
}


# =========================
# 6. QC / DIAGNOSTIC TABLE
# =========================
all_clocks <- c(
  "kdm_advance","phenoage_advance","phenoage_global_advance",
  "hd","hd_log",
  paste0("pheno_", names(subclocks), "_advance")
)

qc_table <- lapply(all_clocks, function(v){
  if(!v %in% names(data)) return(NULL)
  x <- data[[v]]
  data.frame(
    variable   = v,
    n_total    = nrow(data),
    n_nonNA    = sum(!is.na(x)),
    n_finite   = sum(is.finite(x), na.rm=TRUE),
    sd         = sd(x, na.rm=TRUE),
    frac_valid = mean(is.finite(x), na.rm=TRUE)
  )
}) %>% bind_rows()

print(qc_table[order(qc_table$frac_valid),], row.names=FALSE)

valid_clocks <- qc_table %>%
  filter(n_finite > 2000, sd > 0) %>%
  pull(variable)

# =========================
# 7. RESIDUALIZATION (Population Ωᵢ)
# =========================
for(v in valid_clocks){
  f <- as.formula(paste(v,"~ age"))
  data[[paste0(v,"_resid")]] <- resid(lm(f, data=data, na.action = na.exclude))
}

# =========================
# 8. HAZARD RATIO / SURVIVAL ANALYSIS
# =========================
if(all(c("mortstat","permth_int") %in% names(data))){
  
  surv_results <- lapply(valid_clocks, function(v){
    dv <- paste0(v,"_resid")
    sub <- data %>% select(mortstat, permth_int, all_of(dv)) %>% na.omit()
    if(nrow(sub) < 2000) return(NULL)
    cox <- coxph(Surv(permth_int, mortstat) ~ sub[[dv]], data=sub)
    data.frame(
      clock  = v,
      N      = nrow(sub),
      HR     = exp(coef(cox)),
      LCL    = exp(confint(cox)[1]),
      UCL    = exp(confint(cox)[2]),
      Cindex = summary(cox)$concordance[1],
      AIC    = AIC(cox)
    )
  }) %>% bind_rows() %>% arrange(desc(Cindex))
  
  print(surv_results)
  
} else {
  message(">> mortality columns not found, skipping HR analysis")
}

# =========================
# 9. CORRELATION BETWEEN SUBCLOCKS, BAA Y EDAD
# =========================
# Residuals
resid_vars <- paste0(valid_clocks, "_resid")
resid_vars <- intersect(resid_vars, names(data))
cor_resid <- cor(data[, resid_vars], use="pairwise.complete.obs")
print(round(cor_resid,2))

# BAA (advancement)
advance_vars <- names(data)[grepl("_advance$", names(data))]
cor_baa <- cor(data[, advance_vars], use="pairwise.complete.obs")
print(round(cor_baa,2))

# Correlación de BAA con edad
if("age" %in% names(data)){
  cor_age_baa <- cor(data$age, data[, advance_vars], use="pairwise.complete.obs")
  print(round(cor_age_baa,2))
}

# =========================
# 10. HEATMAP / VISUALIZACIÓN DE CORRELACIÓN
# =========================
plot.new()
if(dev.cur() > 1) dev.off()

# Residuals heatmap
short_names <- substr(resid_vars,1,20)
rownames(cor_resid) <- short_names
colnames(cor_resid) <- short_names

corrplot(cor_resid,
         method="color",
         type="upper",
         tl.col="black",
         tl.srt=45,
         addCoef.col="white",
         number.cex=0.6,
         tl.cex=0.7,
         col=colorRampPalette(c("blue","white","red"))(200),
         mar=c(0,0,1,0))
title("Correlación entre relojes residualizados", line=1)

# BAA heatmap
short_names_baa <- substr(advance_vars,1,20)
rownames(cor_baa) <- short_names_baa
colnames(cor_baa) <- short_names_baa

corrplot(cor_baa,
         method="color",
         type="upper",
         tl.col="black",
         tl.srt=45,
         addCoef.col="white",
         number.cex=0.6,
         tl.cex=0.7,
         col=colorRampPalette(c("blue","white","red"))(200),
         mar=c(0,0,1,0))
title("Correlación entre BAA (advancement)", line=1)
