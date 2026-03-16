# Diagnóstico: ¿por qué hema_glucose falla en 93% de casos?

library(BioAge)
library(dplyr)

PROJECT_ROOT <- getwd()
source(file.path(PROJECT_ROOT, "R/config.R"))
source(file.path(PROJECT_ROOT, "R/export.R"))

# Load bundle
bundle <- load_bundle(file.path(PROJECT_ROOT, "bioage_deployment_bundle.rds"))
fit <- bundle$subclock_models$hema_glucose

# 1. Print model coefficients
cat("=== hema_glucose model coefficients ===\n")
print(fit$coef)
cat("\nm_n =", fit$m_n, "\n")
cat("m_d =", fit$m_d, "\n")
cat("BA_n =", fit$BA_n, "\n")
cat("BA_d =", fit$BA_d, "\n")
cat("BA_i =", fit$BA_i, "\n")

# Compare with hema_integrated
fit_hi <- bundle$subclock_models$hema_integrated
cat("\n=== hema_integrated model coefficients ===\n")
print(fit_hi$coef)
cat("\nm_n =", fit_hi$m_n, "\n")
cat("m_d =", fit_hi$m_d, "\n")

# 2. Load clinical data and compute xb manually
f1 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "phenoage_master_PARTE2.csv")
raw <- read.csv(f1, stringsAsFactors = FALSE)
if (file.exists(f2)) raw <- rbind(raw, read.csv(f2, stringsAsFactors = FALSE))

df <- raw %>% mutate(
  age = as.numeric(age),
  glucose = as.numeric(glucose),
  lymph = as.numeric(lymphocyte),
  mcv = as.numeric(mcv),
  rdw = as.numeric(rdw),
  wbc = ifelse(as.numeric(wbc) > 300, as.numeric(wbc)/1000, as.numeric(wbc)),
  rbc = as.numeric(rbc),
  rbc = ifelse(rbc > 100, rbc/1e6, rbc)
)
df <- df[!is.na(df$age) & df$age >= 1 & df$age <= 120, ]

bm <- c("rdw", "mcv", "rbc", "wbc", "lymph", "glucose")
bm_age <- c(bm, "age")
complete <- complete.cases(df[, bm_age])
dat <- df[complete, bm_age]
cat("\n=== Complete cases:", nrow(dat), "===\n")

# 3. Compute xb manually
coef <- fit$coef
cat("\n=== Glucose values in Argentine cohort ===\n")
cat("Range:", range(dat$glucose), "\n")
cat("Mean:", mean(dat$glucose), "\n")
cat("Median:", median(dat$glucose), "\n")
cat("SD:", sd(dat$glucose), "\n")
q <- quantile(dat$glucose, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
cat("Quantiles:\n")
print(q)

# NHANES glucose for comparison
nhanes <- phenoage_nhanes(BIOMARKERS_PHENO_LEVINE)$data
cat("\n=== Glucose values in NHANES IV ===\n")
cat("Range:", range(nhanes$glucose, na.rm=TRUE), "\n")
cat("Mean:", mean(nhanes$glucose, na.rm=TRUE), "\n")
cat("Median:", median(nhanes$glucose, na.rm=TRUE), "\n")
q2 <- quantile(nhanes$glucose, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm=TRUE)
cat("Quantiles:\n")
print(q2)

# 4. Manual xb computation
bm_dat <- t(dat[, bm_age])
n1 <- bm_dat
for (r in 1:nrow(n1)) {
  x <- rownames(n1)[r]
  n1[r,] <- bm_dat[x,] * coef[x, "coef"]
}
xb <- apply(n1, 2, sum) + coef[2, "coef"]  # row 2 = "rate" = intercept

cat("\n=== xb (linear predictor) distribution ===\n")
cat("Range:", range(xb, na.rm=TRUE), "\n")
cat("Mean:", mean(xb, na.rm=TRUE), "\n")
cat("Quantiles:\n")
print(quantile(xb, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm=TRUE))
cat("\nexp(xb) overflow (Inf):", sum(is.infinite(exp(xb))), "/", length(xb), "\n")
cat("exp(xb) > 1e300:", sum(exp(xb) > 1e300, na.rm=TRUE), "\n")

# 5. Check which row is "rate" vs "shape"
cat("\n=== Coefficient row names ===\n")
cat(rownames(fit$coef), sep="\n")
cat("\nRow 1 =", rownames(fit$coef)[1], "=", fit$coef[1, "coef"], "\n")
cat("Row 2 =", rownames(fit$coef)[2], "=", fit$coef[2, "coef"], "\n")

# 6. Compute m and phenoage
m <- 1 - exp((fit$m_n * exp(xb)) / fit$m_d)
pa <- (log(fit$BA_n * log(1 - m)) / fit$BA_d) + fit$BA_i
n_valid <- sum(!is.na(pa) & is.finite(pa))
cat("\n=== PhenoAge results ===\n")
cat("Valid:", n_valid, "/", length(pa), "\n")
cat("% valid:", round(100*n_valid/length(pa), 1), "\n")

# 7. What glucose values produce valid results?
valid_mask <- !is.na(pa) & is.finite(pa)
cat("\n=== Glucose in VALID vs INVALID results ===\n")
cat("Valid glucose - Mean:", mean(dat$glucose[valid_mask]), 
    " Range:", range(dat$glucose[valid_mask]), "\n")
cat("Invalid glucose - Mean:", mean(dat$glucose[!valid_mask]),
    " Range:", range(dat$glucose[!valid_mask]), "\n")

# 8. Contribution of each term to xb
cat("\n=== Mean contribution of each biomarker to xb ===\n")
for (v in bm_age) {
  contrib <- dat[[v]] * coef[v, "coef"]
  cat(sprintf("  %-10s coef=%.6f  mean_val=%.1f  mean_contrib=%.3f\n",
              v, coef[v,"coef"], mean(dat[[v]], na.rm=TRUE), mean(contrib, na.rm=TRUE)))
}
cat(sprintf("  %-10s (rate/intercept) = %.3f\n", "intercept", coef[2,"coef"]))
