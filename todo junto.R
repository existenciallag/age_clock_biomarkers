##########################################################
# BIOLOGICAL AGING – FULL PIPELINE
# NHANES III → NHANES IV
# Canonical clocks + Global PhenoAge + Subclocks + Residuals + HR
##########################################################

# =========================
# 0. LIBRERÍAS
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
pheno_g <- phenoage_nhanes(panel_global) # Global check

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
# 4. ASSEMBLE NHANES IV DATA
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
# 6. CORRELATION CHECKS
# =========================
# BA vs Chronological Age
ba_vars <- c("kdm","phenoage","phenoage_global", paste0("pheno_", names(subclocks)))
ba_vars <- intersect(ba_vars, names(data))
cor_ba_age <- cor(data[, ba_vars], data$age, use="pairwise.complete.obs")
print("Correlación BA vs Edad:")
print(round(cor_ba_age,2))

# BA vs BAA
baa_vars <- c("kdm_advance","phenoage_advance","phenoage_global_advance",
              paste0("pheno_", names(subclocks), "_advance"))
baa_vars <- intersect(baa_vars, names(data))
cor_ba_baa <- cor(data[, ba_vars], data[, baa_vars], use="pairwise.complete.obs")
print("Correlación BA vs BAA:")
print(round(cor_ba_baa,2))

# =========================
# 7. QC / DIAGNOSTIC TABLE
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
# 8. RESIDUALIZATION (BAA residual)
# =========================
for(v in valid_clocks){
  f <- as.formula(paste(v,"~ age"))
  data[[paste0(v,"_resid")]] <- resid(lm(f, data=data, na.action=na.exclude))
}

# =========================
# 9. HAZARD RATIO / SURVIVAL ANALYSIS (seguro)
# =========================
surv_results <- NULL

if(all(c("mortstat","permth_int") %in% names(data))){
  
  surv_results <- lapply(valid_clocks, function(v){
    dv <- paste0(v,"_resid")
    sub <- data %>% select(mortstat, permth_int, all_of(dv)) %>% na.omit()
    if(nrow(sub) < 2000) return(NULL)  # aseguramos suficiente tamaño
    
    # Cox model
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
  
  if(nrow(surv_results) == 0){
    message(">> HR analysis could not run: not enough data per clock")
  } else {
    print(surv_results)
  }
  
} else {
  message(">> mortality columns not found. HR analysis skipped. Returning empty table")
  
  # Crear placeholder vacío con columnas correctas
  surv_results <- data.frame(
    clock  = character(0),
    N      = integer(0),
    HR     = numeric(0),
    LCL    = numeric(0),
    UCL    = numeric(0),
    Cindex = numeric(0),
    AIC    = numeric(0)
  )
}

# =========================
# 10. CORRELATION MATRIX / HEATMAP
# =========================
resid_vars <- paste0(valid_clocks, "_resid")
resid_vars <- intersect(resid_vars, names(data))
cor_resid <- cor(data[, resid_vars], use="pairwise.complete.obs")
print(round(cor_resid,2))

advance_vars <- names(data)[grepl("_advance$", names(data))]
cor_baa <- cor(data[, advance_vars], use="pairwise.complete.obs")
print(round(cor_baa,2))

if("age" %in% names(data)){
  cor_age_baa <- cor(data$age, data[, advance_vars], use="pairwise.complete.obs")
  print(round(cor_age_baa,2))
}

# Residuals heatmap
rownames(cor_resid) <- colnames(cor_resid) <- substr(resid_vars,1,20)
corrplot(cor_resid, method="color", type="upper", tl.col="black",
         tl.srt=45, addCoef.col="white", number.cex=0.6,
         tl.cex=0.7, col=colorRampPalette(c("blue","white","red"))(200))
title("Correlación entre relojes residualizados", line=1)

# BAA heatmap
rownames(cor_baa) <- colnames(cor_baa) <- substr(advance_vars,1,20)
corrplot(cor_baa, method="color", type="upper", tl.col="black",
         tl.srt=45, addCoef.col="white", number.cex=0.6,
         tl.cex=0.7, col=colorRampPalette(c("blue","white","red"))(200))
title("Correlación entre BAA (advancement)", line=1)


