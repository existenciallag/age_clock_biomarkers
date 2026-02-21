# Biological Aging Clock Analysis

Multi-clock pipeline for computing biological age from routine blood biomarkers using the [BioAge](https://github.com/dayoonkwon/BioAge) package. Trains PhenoAge, KDM, and Homeostatic Dysregulation clocks on NHANES III, projects to NHANES IV, and evaluates mortality prediction through Cox proportional hazards models.

Based on:
- **Levine ME (2018)** - PhenoAge: An aging measure based on blood biomarkers
- **Klemera & Doubal (2006)** - KDM biological age
- **Cohen AA et al. (2014)** - Homeostatic Dysregulation

## Quick Start

```r
# 1. Install dependencies
install.packages(c("dplyr", "purrr", "survival", "ggplot2", "tidyr",
                    "tibble", "gridExtra", "corrplot"))
devtools::install_github("dayoonkwon/BioAge")

# 2. Run the full pipeline
source("run_pipeline.R")

# 3. Generate PDF report with all statistics and plots
source("generate_report_pdf.R")

# 4. Score new patients using the trained models
source("run_clinical.R")
```

## Project Structure

```
age_clock_biomarkers/
  run_pipeline.R            Main orchestration (runs steps 1-8)
  run_clinical.R            Score new patients with trained models
  generate_report_pdf.R     PDF report: tables, plots, equations, betas
  phenoage.Rproj            RStudio project file
  R/
    config.R                Panel definitions, thresholds, constants
    train.R                 Train all biological age clocks
    assemble.R              Merge clock outputs into one dataset
    qc.R                    Quality control diagnostics
    residualize.R           Residualize against age, Z-score
    survival.R              Cox models, HR tables, forest plots
    correlation.R           Correlation matrices, heatmaps
    stratify.R              Risk stratification, KM curves
    export.R                Deployment bundle for clinical use
    clinical_score.R        Score new patients from bundle
```

## Pipeline Steps

| Step | Module | What it does |
|------|--------|-------------|
| 1 | `train.R` | Train PhenoAge (original + global), KDM, HD, and 5 sub-clocks on NHANES III |
| 2 | `assemble.R` | Merge all clock projections into a single NHANES IV dataset |
| 3 | `qc.R` | Validate clocks (N, SD, PhenoAge vs age correlation ~0.94) |
| 4 | `residualize.R` | Regress out chronological age, Z-score standardize |
| 5 | `survival.R` | Cox proportional hazards models, hazard ratio tables |
| 6 | `correlation.R` | Inter-clock correlations, independence analysis |
| 7 | `stratify.R` | Risk tertiles, Kaplan-Meier curves, categorical Cox |
| 8 | `export.R` | Save deployment bundle with all coefficients |

## Clocks Trained

### Canonical Clocks

| Clock | Method | Biomarkers | Output |
|-------|--------|------------|--------|
| **PhenoAge Original** | Gompertz mortality model | 12 (albumin, ALP, CRP, cholesterol, creatinine, HbA1c, SBP, BUN, uric acid, lymphocytes, MCV, WBC) | Biological age in years |
| **PhenoAge Global** | Same as original | Same 12 (configurable) | Biological age in years |
| **KDM** | Klemera-Doubal weighted regression | 12 shared biomarkers | Biological age in years |
| **HD** | Mahalanobis distance from young reference | 12 shared biomarkers | Dysregulation score (log-transformed) |

### Sub-Clocks (organ-system PhenoAge models)

| Sub-clock | System | Biomarkers |
|-----------|--------|------------|
| `hepatic_enzime_insulin` | Hepatic | ALP, GGT, insulin |
| `hepatic_lipid` | Hepatic | triglycerides, total cholesterol |
| `hema_integrated` | Hematologic | RDW, MCV, RBC, WBC, lymphocytes |
| `micronutrient_methylation` | Micronutrient | vitamin B12, HbA1c, RDW |
| `renal_A` | Renal | creatinine (log), BUN |

## Module Reference

### R/config.R - Configuration

All biomarker panel definitions and thresholds. Edit this file to customize panels.

| Constant | Description |
|----------|-------------|
| `BIOMARKERS_PHENO_ORIG` | 12-biomarker PhenoAge panel |
| `BIOMARKERS_KDM_HD` | Shared KDM/HD panel |
| `PANEL_GLOBAL` | Global PhenoAge panel (default = original) |
| `PANELS_SUB` | Named list of sub-clock biomarker vectors |
| `PANEL_SYSTEM` | Maps sub-clock names to physiological systems |
| `QC_MIN_N` | Minimum observations for valid clock (2000) |
| `QC_MIN_SD` | Minimum SD for valid clock (0) |
| `COX_MIN_N` | Minimum N for Cox model (200) |

### R/train.R - Clock Training

```r
train_all_clocks(panels_sub, train_kdm, train_hd, verbose)
```
Trains all biological age clocks using `BioAge::phenoage_nhanes()`, `BioAge::kdm_nhanes()`, and `BioAge::hd_nhanes()`. Returns a named list with `$pheno_orig`, `$pheno_global`, `$kdm`, `$hd`, `$subclocks`, `$failed`.

```r
train_subclock(biomarkers, name)
```
Train a single PhenoAge sub-clock. Returns the BioAge fit object or NULL on failure.

### R/assemble.R - Dataset Assembly

```r
assemble_data(clocks)
```
Merges all clock projections into a single data.frame keyed by `sampleID`. Auto-detects mortality columns (`mortstat`, `permth_int`). Computes advancement (BA - CA) for each clock.

```r
list_ba_columns(data)        # Raw biological age columns
list_advance_columns(data)   # Advancement (BA - CA) columns
list_clock_columns(data)     # All clock-related columns
```

### R/qc.R - Quality Control

```r
qc_table(data, clock_vars)
```
Returns a data.frame with N, SD, and fraction valid for each clock.

```r
valid_clocks(qc, min_n, min_sd)
```
Filters QC table to clocks meeting minimum thresholds. Returns character vector of valid column names.

```r
qc_pheno_age_cor(data)
```
Prints Pearson correlation between PhenoAge original and chronological age (expected ~0.94).

```r
qc_overlap(data, clock_vars)
```
Reports how many subjects have complete data across all clocks.

### R/residualize.R - Residualization and Z-Scoring

```r
residualize(y, age)     # lm(y ~ age) residuals
zscore(x)               # Standardize to mean=0, SD=1
add_residuals(data, clock_vars)  # Add _resid and _resid_z columns
add_zscores(data, clock_vars)    # Add z_ prefix columns
list_resid_columns(data)         # Find all _resid columns
list_resid_z_columns(data)       # Find all _resid_z columns
```

**Why residualize?** Raw advancement (BA - CA) retains residual correlation with chronological age, biasing Cox models. Residualization via `lm(advancement ~ age)` removes this confound. Z-scoring makes hazard ratios interpretable as "per +1 SD of biological deviation."

### R/survival.R - Cox Models and Forest Plots

```r
fit_cox_single(data, predictor, adjust_gender, min_n)
```
Fits `Surv(time, status) ~ age + predictor`. Returns a one-row data.frame with HR, 95% CI, p-value, C-index, and AIC.

```r
cox_raw(data, clock_vars, adjust_gender)
```
Batch Cox models on Z-scored BA columns (not residualized).

```r
cox_residual(data, clock_vars, adjust_gender)
```
Batch Cox models on residualized + Z-scored columns.

```r
hr_table(data, clock_vars, adjust_gender)
```
Combined table of raw + residual HR results for all clocks.

```r
forest_plot(hr_df, title, subtitle, colour)
```
ggplot2 forest plot of hazard ratios with 95% CI.

### R/correlation.R - Correlation Analysis

```r
cor_ba_vs_age(data, ba_vars)        # BA vs chronological age (calibration)
cor_matrix(data, vars)              # Pairwise correlation matrix
cor_global_vs_subs(data, ref_col)   # Global PhenoAge vs each sub-clock
discordance_analysis(data, sub_vars) # Pairwise Z-score differences
plot_heatmap_corrplot(mat, title)   # corrplot heatmap
plot_heatmap_gg(mat, title)         # ggplot2 tile heatmap fallback
```

### R/stratify.R - Risk Stratification

```r
add_risk_groups(data, resid_col, n_groups, labels)
```
Assigns subjects to risk tertiles based on residual biological age.

```r
plot_km(data, group_col, colours, title)
```
Kaplan-Meier survival curves by risk group.

```r
cox_categorical(data, group_col)
```
Cox model with risk group as categorical predictor (tests tertile 3 vs tertile 1).

```r
stratify_multi(data, resid_cols, n_groups)
```
Multi-clock stratification: runs categorical Cox for each residualized clock.

```r
plot_resid_density(data, resid_cols)
```
Faceted density plots of residual age distributions.

### R/export.R - Deployment Bundle

```r
build_deployment_bundle(clocks, data, clock_vars, hr_results, qc)
```
Creates a self-contained list with trained model objects, panel definitions, population reference stats, residualization coefficients, and HR results.

```r
save_bundle(bundle, path)   # Save to .rds
load_bundle(path)           # Load from .rds
extract_population_stats(data, clock_vars)   # Mean/SD for Z-scoring
extract_resid_coefficients(data, clock_vars) # lm intercept + slope
```

### R/clinical_score.R - Score New Patients

```r
score_patients(patient_data, bundle, name_map)
```
Applies all trained clocks to new patient lab data. Computes BA, advancement, residuals, and Z-scores using saved population statistics.

```r
score_phenoage(new_data, fit, biomarkers)
```
PhenoAge scoring from first principles (Levine 2018 Gompertz formula):
1. `xb = intercept + sum(beta_i * biomarker_i)`
2. `M = 1 - exp(-exp(xb) * (exp(120 * 0.090165) - 1) / 0.090165)`
3. `PhenoAge = 141.50225 + ln(-0.00553 * ln(1 - M)) / 0.090165`

```r
map_biomarker_names(patient_data, name_map)
```
Renames columns from clinical lab format to NHANES variable names.

```r
clinical_report_table(scored_data, bundle)
```
Per-patient summary with BA, advancement, residual, Z-score, and risk level for each clock.

## Output Files

| File | Contents |
|------|----------|
| `bioage_deployment_bundle.rds` | Trained models, coefficients, population stats |
| `bioage_full_report.pdf` | Complete statistical report (tables + plots) |

## PDF Report Contents

The report (`generate_report_pdf.R`) includes 10 sections:

1. **Biomarker panel definitions** - which biomarkers per clock
2. **Beta coefficients** - trained weights for every model
3. **Scoring equations** - PhenoAge, KDM, HD, Cox formulas
4. **QC summary** - N, SD, fraction valid per clock
5. **Population reference** - mean/SD for Z-scoring new patients
6. **Residualization coefficients** - intercept + slope for age regression
7. **Hazard ratios** - table + forest plots (raw + residual models)
8. **Diagnostic plots** - scatter, heatmaps, density, boxplots, KM curves
9. **Correlation tables** - BA vs age, global vs sub-clocks, discordance
10. **Clinical reference card** - step-by-step scoring guide

## Interpreting Results

### Z-Scores
| Z-score | Interpretation |
|---------|---------------|
| Z <= -1.0 | Biologically younger than expected (low risk) |
| -1.0 < Z < 1.0 | Average biological aging |
| Z >= 1.0 | Biologically older than expected (elevated risk) |

### Hazard Ratios
| HR | Interpretation |
|----|---------------|
| HR = 1.00 | No association with mortality |
| HR = 1.10 | 10% higher mortality risk per +1 SD |
| HR = 1.20 | 20% higher mortality risk per +1 SD |

### Sub-Clock Discordance
When sub-clocks show different Z-scores (e.g., liver old + kidney young), this suggests **organ-specific** accelerated aging and may guide targeted interventions.

## Dependencies

| Package | Purpose |
|---------|---------|
| [BioAge](https://github.com/dayoonkwon/BioAge) | Clock training (phenoage_nhanes, kdm_nhanes, hd_nhanes) |
| dplyr | Data manipulation |
| purrr | Functional programming |
| survival | Cox proportional hazards models |
| ggplot2 | Forest plots, scatter, density, heatmaps |
| tidyr | Data reshaping (pivot_longer) |
| tibble | Data frame utilities |
| gridExtra | Formatted PDF tables |
| corrplot | Correlation heatmaps (optional) |

## References

- Levine ME et al. (2018). "An epigenetic biomarker of aging for lifespan and healthspan." *Aging*, 10(4), 573-591.
- Klemera P, Doubal S (2006). "A new approach to the concept and computation of biological age." *Mechanisms of Ageing and Development*, 127(3), 240-248.
- Cohen AA et al. (2014). "Statistical distance as a measure of physiological dysregulation." *Experimental Gerontology*, 46(C), 1-10.
- Kwon D, Belsky DW (2021). "A toolkit for quantification of biological age from blood chemistry and organ function test data: BioAge." *GeroScience*, 43, 2795-2808.
