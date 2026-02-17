###############################################################################
# generate_report_pdf.R — Generate a comprehensive PDF report
#
# Produces a multi-page PDF with:
#   - All trained clock beta coefficients and equations
#   - QC summary table
#   - Hazard ratio results (raw + residual)
#   - Correlation matrices
#   - Population reference statistics
#   - Residualization coefficients
#   - Forest plot, heatmaps, KM curves, density plots
#
# Usage:
#   source("run_pipeline.R")          # Run the pipeline first
#   source("generate_report_pdf.R")   # Then generate the PDF
#
# Or standalone (loads a saved bundle):
#   source("generate_report_pdf.R")
#
# Output: bioage_full_report.pdf (in the project root)
###############################################################################

# ======================================================================
# 0. AUTO-DETECT PROJECT ROOT
# ======================================================================

REPORT_ROOT <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) getwd())
REPORT_ROOT <- normalizePath(REPORT_ROOT)

# Source modules if not already loaded
if (!exists("BIOMARKERS_PHENO_ORIG")) {
  source(file.path(REPORT_ROOT, "R", "config.R"), local = FALSE)
}

# ======================================================================
# 1. DETERMINE DATA SOURCE
# ======================================================================

# Option A: pipeline_results exists in memory (after run_pipeline.R)
# Option B: load from saved bundle
if (!exists("pipeline_results")) {
  bundle_path <- file.path(REPORT_ROOT, "bioage_deployment_bundle.rds")
  if (file.exists(bundle_path)) {
    message(">> Loading saved deployment bundle for report...")
    source(file.path(REPORT_ROOT, "R", "export.R"), local = FALSE)
    saved_bundle <- load_bundle(bundle_path)
    pipeline_results <- list(bundle = saved_bundle)
  } else {
    stop("No pipeline_results in memory and no saved bundle found.\n",
         "  Run the pipeline first: source('run_pipeline.R')")
  }
}

# ======================================================================
# 2. HELPER FUNCTIONS FOR PDF TEXT LAYOUT
# ======================================================================

suppressPackageStartupMessages({
  library(grDevices)
  library(graphics)
})

# Write a title page
write_title_page <- function() {
  plot.new()
  par(mar = c(0, 0, 0, 0))
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  text(0.5, 0.85, "BIOLOGICAL AGING CLOCK ANALYSIS",
       cex = 2.0, font = 2, col = "navy")
  text(0.5, 0.77, "Complete Statistical Report",
       cex = 1.5, font = 3, col = "gray30")

  text(0.5, 0.63, "Multi-Clock Pipeline (PhenoAge + KDM + HD)",
       cex = 1.1, col = "gray20")
  text(0.5, 0.57, "Based on BioAge (Kwon) & Levine (2018)",
       cex = 1.0, col = "gray40")

  text(0.5, 0.42, paste("Generated:", Sys.time()), cex = 0.9, col = "gray40")
  text(0.5, 0.37, paste("R version:", R.version.string), cex = 0.9, col = "gray40")

  if (!is.null(pipeline_results$bundle$meta$n_subjects)) {
    text(0.5, 0.30, paste("Training population: NHANES IV, N =",
                           pipeline_results$bundle$meta$n_subjects),
         cex = 1.0, col = "gray20")
  }

  # Clocks list
  if (!is.null(pipeline_results$bundle$meta$clocks_trained)) {
    ct <- pipeline_results$bundle$meta$clocks_trained
    text(0.5, 0.20, paste("Clocks trained:", length(ct)), cex = 1.0, col = "gray20")
    y0 <- 0.15
    for (i in seq_along(ct)) {
      text(0.5, y0 - (i - 1) * 0.03, ct[i], cex = 0.8, col = "gray40")
    }
  }
}

# Write a section header page
write_section_header <- function(section_num, title, subtitle = NULL) {
  plot.new()
  par(mar = c(0, 0, 0, 0))
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  rect(0, 0.7, 1, 0.9, col = "navy", border = NA)
  text(0.5, 0.80, paste("Section", section_num), cex = 1.8, font = 2, col = "white")
  text(0.5, 0.60, title, cex = 1.6, font = 2, col = "navy")
  if (!is.null(subtitle)) {
    text(0.5, 0.52, subtitle, cex = 1.0, font = 3, col = "gray40")
  }
}

# Write a data.frame as a text table on a plot
write_table_page <- function(df, title, subtitle = NULL, cex_text = 0.65,
                              max_rows = 40) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  # Split into pages if too many rows
  n_pages <- ceiling(nrow(df) / max_rows)

  for (page in seq_len(n_pages)) {
    start <- (page - 1) * max_rows + 1
    end   <- min(page * max_rows, nrow(df))
    chunk <- df[start:end, , drop = FALSE]

    plot.new()
    par(mar = c(1, 1, 3, 1))
    plot.window(xlim = c(0, 1), ylim = c(0, 1))

    page_title <- title
    if (n_pages > 1) page_title <- paste0(title, " (", page, "/", n_pages, ")")
    title(page_title, cex.main = 1.2, col.main = "navy", font.main = 2)

    if (!is.null(subtitle) && page == 1) {
      mtext(subtitle, side = 3, line = 0, cex = 0.8, col = "gray40", font = 3)
    }

    # Format numeric columns
    for (j in seq_len(ncol(chunk))) {
      if (is.numeric(chunk[[j]])) {
        chunk[[j]] <- ifelse(
          abs(chunk[[j]]) < 0.001 & chunk[[j]] != 0,
          formatC(chunk[[j]], format = "e", digits = 3),
          formatC(chunk[[j]], format = "f", digits = 4)
        )
      }
    }

    # Build text lines
    header <- paste(formatC(names(chunk), width = 18, flag = "-"), collapse = "")
    sep    <- paste(rep("-", nchar(header)), collapse = "")

    lines <- c(header, sep)
    for (i in seq_len(nrow(chunk))) {
      vals <- vapply(chunk[i, ], function(x) {
        formatC(as.character(x), width = 18, flag = "-")
      }, character(1))
      lines <- c(lines, paste(vals, collapse = ""))
    }

    # Write lines
    n_lines <- length(lines)
    y_positions <- seq(0.95, 0.05, length.out = max(n_lines, max_rows + 5))

    for (i in seq_along(lines)) {
      font_val <- if (i <= 2) 2 else 1
      col_val  <- if (i <= 2) "navy" else "black"
      text(0.02, y_positions[i], lines[i], adj = c(0, 1),
           family = "mono", cex = cex_text, font = font_val, col = col_val)
    }
  }
}

# Write text content page
write_text_page <- function(title, text_lines, cex_text = 0.75) {
  plot.new()
  par(mar = c(1, 2, 3, 1))
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  title(title, cex.main = 1.3, col.main = "navy", font.main = 2)

  y <- 0.92
  for (ln in text_lines) {
    # Check if this is a header line (starts with "==")
    if (grepl("^==", ln)) {
      text(0.02, y, ln, adj = c(0, 1), family = "mono",
           cex = cex_text, font = 2, col = "navy")
    } else if (grepl("^--", ln)) {
      text(0.02, y, ln, adj = c(0, 1), family = "mono",
           cex = cex_text, font = 2, col = "steelblue")
    } else {
      text(0.02, y, ln, adj = c(0, 1), family = "mono",
           cex = cex_text, col = "black")
    }
    y <- y - 0.028
    if (y < 0.02) break
  }
}

# ======================================================================
# 3. EXTRACT ALL DATA FOR THE REPORT
# ======================================================================

bundle <- pipeline_results$bundle

# ---- A. Beta coefficients from each clock model ----

extract_model_betas <- function(fit, clock_name) {
  if (is.null(fit)) return(NULL)

  # BioAge stores coefficients in different places depending on model type
  coefs <- NULL

  # Try $coefficients directly
  if (!is.null(fit$coefficients)) {
    coefs <- fit$coefficients
  }
  # Try $coef
  if (is.null(coefs) && !is.null(fit$coef)) {
    coefs <- fit$coef
  }
  # Try inside a glm/lm summary structure
  if (is.null(coefs) && !is.null(fit$results) && !is.null(fit$results$coefficients)) {
    coefs <- fit$results$coefficients
  }
  # Try survreg-style
  if (is.null(coefs) && is.list(fit)) {
    # Walk the list looking for a named numeric vector
    for (nm in names(fit)) {
      obj <- fit[[nm]]
      if (is.numeric(obj) && !is.null(names(obj)) && length(obj) > 1) {
        coefs <- obj
        break
      }
    }
  }

  if (is.null(coefs)) {
    return(data.frame(
      clock       = clock_name,
      variable    = "(no coefficients found)",
      beta        = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    clock    = clock_name,
    variable = names(coefs),
    beta     = as.numeric(coefs),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# Extract from all models
all_betas <- list()

if (!is.null(bundle$models$pheno_orig)) {
  all_betas$pheno_orig <- extract_model_betas(bundle$models$pheno_orig, "PhenoAge Original")
}
if (!is.null(bundle$models$pheno_global)) {
  all_betas$pheno_global <- extract_model_betas(bundle$models$pheno_global, "PhenoAge Global")
}
if (!is.null(bundle$models$kdm)) {
  all_betas$kdm <- extract_model_betas(bundle$models$kdm, "KDM BioAge")
}
if (!is.null(bundle$models$hd)) {
  all_betas$hd <- extract_model_betas(bundle$models$hd, "HD")
}
# Sub-clocks
if (!is.null(bundle$subclock_models)) {
  for (nm in names(bundle$subclock_models)) {
    label <- paste0("Sub: ", gsub("_", " ", nm))
    all_betas[[nm]] <- extract_model_betas(bundle$subclock_models[[nm]], label)
  }
}

betas_df <- do.call(rbind, all_betas)
rownames(betas_df) <- NULL

# ---- B. Build equation strings ----

build_equation_lines <- function() {
  lines <- character(0)

  lines <- c(lines,
    "== PHENOAGE SCORING EQUATION (Levine 2018) ==",
    "",
    "Step 1: Linear predictor (Gompertz mortality score)",
    "  xb = b0 + b1*albumin_gL + b2*alp + b3*lncrp + ...",
    "",
    "Step 2: Mortality probability",
    "  M = 1 - exp( -exp(xb) * (exp(120*gamma) - 1) / gamma )",
    "",
    "Step 3: Invert Gompertz to get PhenoAge",
    "  PhenoAge = 141.50225 + ln(-0.00553 * ln(1 - M)) / 0.090165",
    "",
    "  where gamma = 0.090165 (Gompertz aging rate)",
    ""
  )

  lines <- c(lines,
    "== ADVANCEMENT ==",
    "  Advancement = BiologicalAge - ChronologicalAge",
    "",
    "== RESIDUALIZATION ==",
    "  ResidualAge = residuals( lm(Advancement ~ Age) )",
    "  This removes any residual correlation with chronological age",
    "",
    "== Z-SCORING ==",
    "  Z = (ResidualAge - mean) / SD",
    "  Interpretation: Z > 0 = biologically older than expected",
    "                  Z < 0 = biologically younger than expected",
    "",
    "== KDM BIOLOGICAL AGE ==",
    "  Klemera-Doubal method: weighted average of age-regressed",
    "  biomarker estimates. Each biomarker contributes based on",
    "  the strength of its age correlation (weighted by 1/residual_SD).",
    "  KDM_BA = sum(wi * age_i) / sum(wi)",
    "  where age_i = (xi - intercept_i) / slope_i for biomarker i",
    "  and wi = (slope_i / residual_SD_i)^2",
    "",
    "== HOMEOSTATIC DYSREGULATION (HD) ==",
    "  Mahalanobis distance from the young-adult (20-30y) reference",
    "  distribution in biomarker space.",
    "  HD = sqrt( (x - mu_ref)' * Sigma_ref^{-1} * (x - mu_ref) )",
    "  hd_log = log(HD) is often used for normality.",
    "",
    "== COX PROPORTIONAL HAZARDS ==",
    "  h(t) = h0(t) * exp(b_age*Age + b_clock*Z_clock)",
    "  HR = exp(b_clock) = hazard ratio per +1 SD of residual age",
    ""
  )

  lines
}

# ---- C. Panel composition reference ----

build_panel_lines <- function() {
  lines <- character(0)

  lines <- c(lines, "== BIOMARKER PANELS ==", "")

  lines <- c(lines,
    "-- PhenoAge Original (12 biomarkers) --",
    paste("   ", BIOMARKERS_PHENO_ORIG, collapse = "\n"),
    ""
  )

  lines <- c(lines,
    "-- KDM / HD Panel --",
    paste("   ", BIOMARKERS_KDM_HD, collapse = "\n"),
    ""
  )

  if (exists("PANELS_SUB")) {
    for (nm in names(PANELS_SUB)) {
      lines <- c(lines,
        paste0("-- Sub-clock: ", nm, " --"),
        paste("   ", PANELS_SUB[[nm]], collapse = "\n"),
        ""
      )
    }
  }

  lines
}

# ======================================================================
# 4. GENERATE THE PDF
# ======================================================================

pdf_path <- file.path(REPORT_ROOT, "bioage_full_report.pdf")

message(">> Generating PDF report: ", pdf_path)

pdf(pdf_path, width = 11, height = 8.5, family = "Helvetica")

# ---- PAGE 1: Title ----
write_title_page()

# ---- SECTION 1: Panel Definitions ----
write_section_header(1, "Biomarker Panel Definitions",
                     "Which biomarkers go into each clock")

panel_lines <- build_panel_lines()
# Split into chunks that fit on a page (~30 lines each)
for (start in seq(1, length(panel_lines), by = 30)) {
  end <- min(start + 29, length(panel_lines))
  write_text_page("Biomarker Panels", panel_lines[start:end])
}

# ---- SECTION 2: Beta Coefficients ----
write_section_header(2, "Beta Coefficients",
                     "Trained model weights for each clock")

if (!is.null(betas_df) && nrow(betas_df) > 0) {
  # One page per clock
  for (clock_nm in unique(betas_df$clock)) {
    sub <- betas_df[betas_df$clock == clock_nm, c("variable", "beta"), drop = FALSE]
    write_table_page(sub,
                     title = paste("Beta Coefficients:", clock_nm),
                     subtitle = "These are the trained model weights (log-hazard scale for PhenoAge)")
  }
} else {
  write_text_page("Beta Coefficients",
                  c("No coefficient data available from the model fit objects.",
                    "This may happen if BioAge stores fits in a non-standard format.",
                    "The fit objects are preserved in the deployment bundle (.rds)."))
}

# ---- SECTION 3: Equations ----
write_section_header(3, "Scoring Equations",
                     "How to compute each biological age from raw biomarkers")

eq_lines <- build_equation_lines()
for (start in seq(1, length(eq_lines), by = 30)) {
  end <- min(start + 29, length(eq_lines))
  write_text_page("Scoring Equations", eq_lines[start:end])
}

# ---- SECTION 4: QC Summary ----
write_section_header(4, "Quality Control",
                     "QC diagnostics for all trained clocks")

if (!is.null(bundle$qc)) {
  write_table_page(bundle$qc, title = "QC Summary Table",
                   subtitle = "n_finite > 2000 and sd > 0 required for valid clock")
}

if (!is.null(pipeline_results$qc)) {
  write_table_page(pipeline_results$qc, title = "QC Summary Table",
                   subtitle = "n_finite > 2000 and sd > 0 required for valid clock")
}

# ---- SECTION 5: Population Statistics ----
write_section_header(5, "Population Reference Statistics",
                     "Mean, SD, and N for Z-scoring new patients")

if (!is.null(bundle$pop_stats)) {
  write_table_page(bundle$pop_stats,
                   title = "Population Reference (NHANES IV)",
                   subtitle = "Use these values to Z-score new patient data")
}

# ---- SECTION 6: Residualization Coefficients ----
write_section_header(6, "Residualization Coefficients",
                     "lm(advancement ~ age) intercept and slope")

if (!is.null(bundle$resid_coefs)) {
  write_table_page(bundle$resid_coefs,
                   title = "Age Regression Coefficients",
                   subtitle = "residual = advancement - (intercept + slope * age)")
}

# ---- SECTION 7: Hazard Ratios ----
write_section_header(7, "Survival Analysis: Hazard Ratios",
                     "Cox PH models — HR per +1 SD of biological age")

hr_data <- NULL
if (!is.null(bundle$hr_results)) {
  hr_data <- bundle$hr_results
} else if (!is.null(pipeline_results$hr)) {
  hr_data <- pipeline_results$hr
}

if (!is.null(hr_data) && nrow(hr_data) > 0) {
  write_table_page(hr_data,
                   title = "Hazard Ratio Table (All Models)",
                   subtitle = "HR > 1 = higher mortality per +1 SD of biological age")

  # Forest plot
  tryCatch({
    library(ggplot2)
    source(file.path(REPORT_ROOT, "R", "survival.R"), local = FALSE)
    p <- forest_plot(hr_data,
                     title = "Hazard Ratios: All Biological Age Clocks",
                     subtitle = "Per +1 SD (adjusted for chronological age)")
    print(p)
  }, error = function(e) {
    message("  ! Could not generate forest plot: ", e$message)
  })
} else {
  write_text_page("Hazard Ratio Results",
                  c("No hazard ratio results available.",
                    "This is expected if no mortality data was present.",
                    "HR results require NHANES IV mortality follow-up data."))
}

# ---- SECTION 8: Correlation Data ----
write_section_header(8, "Correlation Analysis",
                     "Inter-clock correlations and independence")

# If data is in memory, compute correlations
if (!is.null(pipeline_results$data)) {
  tryCatch({
    source(file.path(REPORT_ROOT, "R", "assemble.R"), local = FALSE)
    source(file.path(REPORT_ROOT, "R", "residualize.R"), local = FALSE)
    source(file.path(REPORT_ROOT, "R", "correlation.R"), local = FALSE)

    d <- pipeline_results$data

    # BA vs age
    ba_vars <- list_ba_columns(d)
    if (length(ba_vars) > 0) {
      cor_age_df <- cor_ba_vs_age(d, ba_vars)
      write_table_page(cor_age_df,
                       title = "Correlation: BA vs Chronological Age",
                       subtitle = "Pearson r (calibration check)")
    }

    # Residual correlation matrix
    resid_cols <- list_resid_columns(d)
    if (length(resid_cols) >= 2) {
      mat <- cor_matrix(d, resid_cols)

      # Write numeric matrix
      mat_df <- as.data.frame(round(mat, 3))
      mat_df <- cbind(clock = rownames(mat_df), mat_df)
      write_table_page(mat_df,
                       title = "Residual Sub-Clock Correlation Matrix",
                       subtitle = "Low correlations = independent aging systems")

      # Heatmap plot
      has_corrplot <- requireNamespace("corrplot", quietly = TRUE)
      if (has_corrplot) {
        plot_heatmap_corrplot(mat, title = "Residual Sub-Clock Independence")
      } else {
        p <- plot_heatmap_gg(mat, "Residual Sub-Clock Independence")
        print(p)
      }
    }

    # Advancement correlation matrix
    adv_cols <- list_advance_columns(d)
    if (length(adv_cols) >= 2) {
      mat_adv <- cor_matrix(d, adv_cols)
      mat_adv_df <- as.data.frame(round(mat_adv, 3))
      mat_adv_df <- cbind(clock = rownames(mat_adv_df), mat_adv_df)
      write_table_page(mat_adv_df,
                       title = "Advancement (BA - CA) Correlation Matrix",
                       subtitle = "System concordance before residualization")

      if (has_corrplot) {
        plot_heatmap_corrplot(mat_adv, title = "BA Advancement Concordance")
      }
    }

    # Global vs sub-clocks
    tryCatch({
      cor_gs <- cor_global_vs_subs(d)
      if (nrow(cor_gs) > 0) {
        write_table_page(cor_gs,
                         title = "Global PhenoAge vs Sub-Clocks",
                         subtitle = "How much each sub-clock tracks the global clock")
      }
    }, error = function(e) NULL)

    # Discordance
    tryCatch({
      disc <- discordance_analysis(d)
      if (!is.null(disc) && nrow(disc) > 0) {
        write_table_page(disc,
                         title = "Pairwise System Discordance",
                         subtitle = "Mean Z-score difference between system pairs")
      }
    }, error = function(e) NULL)

  }, error = function(e) {
    write_text_page("Correlation Analysis",
                    c("Could not compute correlations from pipeline data.",
                      paste("Error:", e$message)))
  })
} else {
  write_text_page("Correlation Analysis",
                  c("Pipeline data not in memory.",
                    "Run source('run_pipeline.R') first, then generate the report."))
}

# ---- SECTION 9: Risk Stratification ----
write_section_header(9, "Risk Stratification",
                     "Kaplan-Meier curves and categorical Cox models")

if (!is.null(pipeline_results$data)) {
  d <- pipeline_results$data
  has_mortality <- all(c("time", "status") %in% names(d)) &&
    !all(is.na(d$status))

  if (has_mortality) {
    tryCatch({
      source(file.path(REPORT_ROOT, "R", "stratify.R"), local = FALSE)

      global_resid <- intersect(
        c("phenoage_orig_advance_resid", "phenoage_global_advance_resid"),
        names(d)
      )

      if (length(global_resid) > 0 && "risk_group" %in% names(d)) {
        # KM plot
        plot_km(d, title = "Survival by Residual Biological Age (PhenoAge)")
      } else if (length(global_resid) > 0) {
        d <- add_risk_groups(d, global_resid[1])
        plot_km(d, title = "Survival by Residual Biological Age (PhenoAge)")
      }

      # Density plots
      resid_cols <- list_resid_columns(d)
      if (length(resid_cols) > 0) {
        p_dens <- plot_resid_density(d, resid_cols)
        print(p_dens)
      }
    }, error = function(e) {
      message("  ! Risk stratification plot error: ", e$message)
    })
  } else {
    write_text_page("Risk Stratification",
                    c("No mortality data available.",
                      "KM curves and categorical Cox models require mortality follow-up."))
  }
}

# ---- SECTION 10: Clinical Use Reference Card ----
write_section_header(10, "Clinical Reference Card",
                     "How to interpret and apply the clocks")

ref_lines <- c(
  "== HOW TO SCORE A NEW PATIENT ==",
  "",
  "1. Collect blood panel: CBC, BMP, LFTs, HbA1c, CRP, Uric Acid",
  "2. Compute PhenoAge using the beta coefficients in Section 2:",
  "     xb = sum(beta_i * biomarker_i)",
  "     M  = 1 - exp(-exp(xb) * (exp(120*0.090165) - 1) / 0.090165)",
  "     PhenoAge = 141.50225 + ln(-0.00553 * ln(1-M)) / 0.090165",
  "3. Advancement = PhenoAge - ChronologicalAge",
  "4. Residualize using coefficients in Section 6:",
  "     residual = advancement - (intercept + slope * age)",
  "5. Z-score using population stats in Section 5:",
  "     Z = (residual - pop_mean) / pop_sd",
  "",
  "== INTERPRETING Z-SCORES ==",
  "",
  "  Z <= -1.0    Biologically younger (low risk)",
  "  -1.0 < Z < 1.0  Average biological aging",
  "  Z >= 1.0     Biologically older (elevated risk)",
  "",
  "== INTERPRETING HAZARD RATIOS ==",
  "",
  "  HR = 1.00    No association with mortality",
  "  HR = 1.10    10% higher mortality risk per +1 SD",
  "  HR = 1.20    20% higher mortality risk per +1 SD",
  "",
  "== SUB-CLOCK SYSTEM INTERPRETATION ==",
  "",
  "  Each sub-clock captures a specific physiological domain:",
  "    hepatic_enzime_insulin  - Liver function / insulin resistance",
  "    hepatic_lipid           - Liver lipid metabolism",
  "    hema_integrated         - Hematologic / bone marrow health",
  "    micronutrient_methylation - B12 / methylation / epigenetic",
  "    renal_A                 - Kidney filtration function",
  "",
  "  Discordant sub-clocks (e.g., liver old + kidney young)",
  "  suggest organ-specific interventions.",
  "",
  "== DEPLOYMENT BUNDLE ==",
  "",
  "  The file 'bioage_deployment_bundle.rds' contains everything",
  "  needed to score new patients without retraining:",
  "    bundle <- readRDS('bioage_deployment_bundle.rds')",
  "    source('R/clinical_score.R')",
  "    results <- score_patients(patient_df, bundle)"
)

for (start in seq(1, length(ref_lines), by = 30)) {
  end <- min(start + 29, length(ref_lines))
  write_text_page("Clinical Reference Card", ref_lines[start:end])
}

# ---- Close PDF ----
dev.off()

message(">> PDF report saved to: ", pdf_path)
message(">> Total pages generated. Open with any PDF viewer.")
