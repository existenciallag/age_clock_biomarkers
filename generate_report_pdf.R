###############################################################################
# generate_report_pdf.R — Full statistical report: tables + all plots
#
# Usage:
#   source("run_pipeline.R")          # Run the pipeline first
#   source("generate_report_pdf.R")   # Then generate the PDF
#
# Output: bioage_full_report.pdf  (all tables + all plots)
###############################################################################

# ======================================================================
# 0. SETUP
# ======================================================================

REPORT_ROOT <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) getwd())
REPORT_ROOT <- normalizePath(REPORT_ROOT)

# Source ALL modules so functions are available
module_files <- c("config.R", "train.R", "assemble.R", "qc.R",
                  "residualize.R", "survival.R", "correlation.R",
                  "stratify.R", "export.R")
for (mf in module_files) {
  fp <- file.path(REPORT_ROOT, "R", mf)
  if (file.exists(fp)) source(fp, local = FALSE)
}

suppressPackageStartupMessages({
  library(grDevices)
  library(graphics)
  library(grid)
  library(ggplot2)
  library(survival)
  library(dplyr)
  library(tidyr)
})

# gridExtra for clean tables — install if missing
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  message(">> Installing gridExtra for formatted tables...")
  install.packages("gridExtra", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(gridExtra)

has_corrplot <- requireNamespace("corrplot", quietly = TRUE)
if (has_corrplot) library(corrplot)

# Verify pipeline_results exists
if (!exists("pipeline_results")) {
  bundle_path <- file.path(REPORT_ROOT, "bioage_deployment_bundle.rds")
  if (file.exists(bundle_path)) {
    message(">> Loading saved deployment bundle...")
    saved_bundle <- load_bundle(bundle_path)
    pipeline_results <- list(bundle = saved_bundle)
  } else {
    stop("Run source('run_pipeline.R') first, then source this script.")
  }
}

bundle <- pipeline_results$bundle
d      <- pipeline_results$data   # May be NULL if loaded from bundle only

# ======================================================================
# HELPER: draw a data.frame as a clean gridExtra table
# ======================================================================

draw_table <- function(df, title, subtitle = NULL, fontsize = 8) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  # Safely convert all columns
  for (j in seq_len(ncol(df))) {
    col <- df[[j]]
    if (is.list(col)) {
      df[[j]] <- vapply(col, function(x) {
        if (is.null(x)) ""
        else if (is.numeric(x) && length(x) == 1) sprintf("%.4g", x)
        else paste(as.character(x), collapse = ",")
      }, character(1))
    } else if (is.numeric(col)) {
      df[[j]] <- ifelse(is.na(col), "",
                        ifelse(abs(col) < 0.0001 & col != 0,
                               sprintf("%.3e", col), sprintf("%.4f", col)))
    } else if (is.factor(col)) {
      df[[j]] <- as.character(col)
    }
  }

  # Theme for gridExtra table
  tt <- ttheme_minimal(
    core    = list(fg_params = list(cex = fontsize / 10, hjust = 1, x = 0.95),
                   bg_params = list(fill = c("gray95", "white"))),
    colhead = list(fg_params = list(cex = fontsize / 10, fontface = "bold",
                                     col = "white"),
                   bg_params = list(fill = "navy")),
    rowhead = list(fg_params = list(cex = fontsize / 10, fontface = "italic"))
  )

  grobs <- list()

  # Title grob
  grobs[[length(grobs) + 1]] <- textGrob(title, gp = gpar(fontsize = 14,
                                          fontface = "bold", col = "navy"))
  if (!is.null(subtitle)) {
    grobs[[length(grobs) + 1]] <- textGrob(subtitle,
                                            gp = gpar(fontsize = 10,
                                                      fontface = "italic",
                                                      col = "gray40"))
  }

  # Split large tables across pages
  max_rows <- 35
  n_pages <- ceiling(nrow(df) / max_rows)

  for (pg in seq_len(n_pages)) {
    start_r <- (pg - 1) * max_rows + 1
    end_r   <- min(pg * max_rows, nrow(df))
    chunk   <- df[start_r:end_r, , drop = FALSE]

    tbl_grob <- tableGrob(chunk, rows = NULL, theme = tt)

    if (pg == 1) {
      grid.newpage()
      all_grobs <- c(grobs, list(tbl_grob))
      heights <- unit(c(rep(0.6, length(grobs)), 1), c(rep("cm", length(grobs)), "null"))
      grid.draw(arrangeGrob(grobs = all_grobs, heights = heights))
    } else {
      grid.newpage()
      pg_title <- textGrob(paste0(title, " (cont. ", pg, "/", n_pages, ")"),
                           gp = gpar(fontsize = 12, fontface = "bold", col = "navy"))
      grid.draw(arrangeGrob(pg_title, tbl_grob,
                             heights = unit(c(0.6, 1), c("cm", "null"))))
    }
  }
}

# ======================================================================
# HELPER: section header
# ======================================================================

draw_section <- function(num, title, subtitle = NULL) {
  grid.newpage()
  pushViewport(viewport(width = 1, height = 1))
  grid.rect(x = 0.5, y = 0.80, width = 1, height = 0.15,
            gp = gpar(fill = "navy", col = NA))
  grid.text(paste("Section", num), x = 0.5, y = 0.80,
            gp = gpar(fontsize = 24, fontface = "bold", col = "white"))
  grid.text(title, x = 0.5, y = 0.60,
            gp = gpar(fontsize = 20, fontface = "bold", col = "navy"))
  if (!is.null(subtitle)) {
    grid.text(subtitle, x = 0.5, y = 0.52,
              gp = gpar(fontsize = 13, fontface = "italic", col = "gray40"))
  }
  popViewport()
}

# ======================================================================
# HELPER: title page
# ======================================================================

draw_title_page <- function() {
  grid.newpage()
  pushViewport(viewport(width = 1, height = 1))
  grid.text("BIOLOGICAL AGING CLOCK ANALYSIS", x = 0.5, y = 0.85,
            gp = gpar(fontsize = 28, fontface = "bold", col = "navy"))
  grid.text("Complete Statistical Report", x = 0.5, y = 0.77,
            gp = gpar(fontsize = 18, fontface = "italic", col = "gray30"))
  grid.text("Multi-Clock Pipeline: PhenoAge + KDM + HD", x = 0.5, y = 0.65,
            gp = gpar(fontsize = 14, col = "gray20"))
  grid.text("Based on BioAge (Kwon) & Levine (2018)", x = 0.5, y = 0.60,
            gp = gpar(fontsize = 12, col = "gray40"))
  grid.text(paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
            x = 0.5, y = 0.48, gp = gpar(fontsize = 11, col = "gray40"))

  n_sub <- if (!is.null(bundle$meta$n_subjects)) bundle$meta$n_subjects else "?"
  grid.text(paste("Training population: NHANES IV, N =", n_sub),
            x = 0.5, y = 0.40, gp = gpar(fontsize = 12, col = "gray20"))

  ct <- bundle$meta$clocks_trained
  if (!is.null(ct)) {
    grid.text(paste("Clocks:", paste(ct, collapse = "  |  ")),
              x = 0.5, y = 0.33, gp = gpar(fontsize = 10, col = "gray40"))
  }
  popViewport()
}

# ======================================================================
# HELPER: text page
# ======================================================================

draw_text_page <- function(title, text_lines) {
  grid.newpage()
  pushViewport(viewport(width = 0.9, height = 0.9))
  grid.text(title, x = 0.5, y = 0.97,
            gp = gpar(fontsize = 16, fontface = "bold", col = "navy"))

  y <- 0.92
  for (ln in text_lines) {
    if (y < 0.02) break
    fp <- if (grepl("^==", ln)) {
      gpar(fontsize = 11, fontface = "bold", col = "navy", fontfamily = "mono")
    } else if (grepl("^--", ln)) {
      gpar(fontsize = 10, fontface = "bold", col = "steelblue", fontfamily = "mono")
    } else {
      gpar(fontsize = 9, col = "black", fontfamily = "mono")
    }
    grid.text(ln, x = 0.02, y = y, just = "left", gp = fp)
    y <- y - 0.028
  }
  popViewport()
}

# ======================================================================
# EXTRACT BETA COEFFICIENTS
# ======================================================================

extract_model_betas <- function(fit, clock_name) {
  if (is.null(fit)) return(NULL)

  collect_numeric <- function(obj, prefix = "") {
    rows <- list()
    if (is.null(obj)) return(rows)
    if (is.numeric(obj) && !is.null(names(obj)) && !is.matrix(obj)) {
      for (i in seq_along(obj)) {
        nm <- if (nzchar(prefix)) paste0(prefix, ".", names(obj)[i]) else names(obj)[i]
        rows[[length(rows) + 1]] <- data.frame(clock = clock_name, variable = nm,
          beta = obj[i], stringsAsFactors = FALSE, row.names = NULL)
      }
      return(rows)
    }
    if (is.numeric(obj) && length(obj) == 1 && is.null(names(obj)) && nzchar(prefix)) {
      rows[[1]] <- data.frame(clock = clock_name, variable = prefix,
        beta = obj, stringsAsFactors = FALSE, row.names = NULL)
      return(rows)
    }
    if (is.matrix(obj) || is.data.frame(obj)) {
      if (nrow(obj) > 50) return(rows)
      col_pick <- intersect(c("coef", "Estimate"), colnames(obj))
      vals <- if (length(col_pick) > 0) as.numeric(obj[, col_pick[1]])
              else as.numeric(obj[, 1])
      nms <- rownames(obj)
      if (is.null(nms)) nms <- paste0("V", seq_along(vals))
      for (i in seq_along(vals)) {
        lab <- if (nzchar(prefix)) paste0(prefix, ".", nms[i]) else nms[i]
        rows[[length(rows) + 1]] <- data.frame(clock = clock_name, variable = lab,
          beta = vals[i], stringsAsFactors = FALSE, row.names = NULL)
      }
      return(rows)
    }
    if (is.list(obj) && !is.null(names(obj))) {
      for (nm in names(obj)) {
        child <- obj[[nm]]
        if (is.function(child) || inherits(child, "formula") ||
            inherits(child, "call")) next
        if (is.data.frame(child) && nrow(child) > 50) next
        child_prefix <- if (nzchar(prefix)) paste0(prefix, ".", nm) else nm
        rows <- c(rows, collect_numeric(child, child_prefix))
      }
    }
    rows
  }

  all_rows <- collect_numeric(fit)
  if (length(all_rows) == 0) {
    return(data.frame(clock = clock_name,
                      variable = paste("(structure:", paste(names(fit), collapse = ", "), ")"),
                      beta = NA_real_, stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, all_rows)
  rownames(out) <- NULL
  out
}

all_betas <- list()
if (!is.null(bundle$models$pheno_orig))
  all_betas$po <- extract_model_betas(bundle$models$pheno_orig, "PhenoAge Original")
if (!is.null(bundle$models$pheno_global))
  all_betas$pg <- extract_model_betas(bundle$models$pheno_global, "PhenoAge Global")
if (!is.null(bundle$models$kdm))
  all_betas$kdm <- extract_model_betas(bundle$models$kdm, "KDM BioAge")
if (!is.null(bundle$models$hd))
  all_betas$hd <- extract_model_betas(bundle$models$hd, "HD")
if (!is.null(bundle$subclock_models)) {
  for (nm in names(bundle$subclock_models)) {
    all_betas[[nm]] <- extract_model_betas(bundle$subclock_models[[nm]],
                                            paste0("Sub: ", gsub("_", " ", nm)))
  }
}
betas_df <- do.call(rbind, all_betas)
if (!is.null(betas_df)) rownames(betas_df) <- NULL

# ======================================================================
# GET / RECOMPUTE HR DATA
# ======================================================================

hr_data <- NULL

# Try stored sources first
tryCatch({
  tmp <- bundle$hr_results
  if (!is.null(tmp) && is.data.frame(tmp) && nrow(tmp) > 0) hr_data <- tmp
}, error = function(e) NULL)

if (is.null(hr_data)) tryCatch({
  tmp <- pipeline_results$hr
  if (!is.null(tmp) && is.data.frame(tmp) && nrow(tmp) > 0) hr_data <- tmp
}, error = function(e) NULL)

# Recompute from data as final fallback
if (is.null(hr_data) && !is.null(d)) {
  has_mort <- all(c("time", "status") %in% names(d)) &&
    any(!is.na(d$status) & d$status == 1)
  if (has_mort) {
    message(">> Recomputing HR from data...")
    tryCatch({
      ba_cols_r  <- list_ba_columns(d)
      rz_avail   <- grep("_resid_z$", names(d), value = TRUE)
      valid_r    <- sub("_resid_z$", "", rz_avail)
      hr_raw_r   <- cox_raw(d, ba_cols_r, adjust_gender = FALSE)
      hr_res_r   <- cox_residual(d, valid_r, adjust_gender = FALSE)
      hr_data    <- rbind(hr_raw_r, hr_res_r)
    }, error = function(e) message("  ! HR recompute: ", e$message))
  }
}

message(">> HR data: ", if (!is.null(hr_data)) paste(nrow(hr_data), "models") else "none")

# ======================================================================
# BUILD ALL GGPLOT OBJECTS (for viewer + PDF)
# ======================================================================

plots <- list()

# --- Plot 1: Forest plot of HR ---
if (!is.null(hr_data) && nrow(hr_data) > 0) {
  plots$forest <- forest_plot(hr_data,
    title = "Hazard Ratios: All Biological Age Clocks",
    subtitle = "Per +1 SD (adjusted for chronological age)")

  # Separate raw vs residual if both present
  if ("model" %in% names(hr_data)) {
    for (mdl in unique(hr_data$model)) {
      sub_hr <- hr_data[hr_data$model == mdl, , drop = FALSE]
      if (nrow(sub_hr) > 0) {
        nm <- paste0("forest_", gsub(" ", "_", tolower(mdl)))
        plots[[nm]] <- forest_plot(sub_hr, title = paste("HR:", mdl),
                                    subtitle = "Per +1 SD of biological age")
      }
    }
  }
}

# --- Plot 2: BA vs Age scatter ---
if (!is.null(d) && "phenoage_orig" %in% names(d)) {
  plots$ba_vs_age <- ggplot(d, aes(x = age, y = phenoage_orig)) +
    geom_point(alpha = 0.15, size = 0.8, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red", linewidth = 0.8) +
    geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
    theme_minimal(base_size = 13) +
    labs(title = "PhenoAge vs Chronological Age",
         subtitle = paste0("r = ", round(cor(d$phenoage_orig, d$age, use = "complete.obs"), 3)),
         x = "Chronological Age", y = "PhenoAge (Biological Age)")
}

# --- Plot 3: Density of residual ages ---
if (!is.null(d)) {
  resid_cols <- list_resid_columns(d)
  if (length(resid_cols) > 0) {
    long_resid <- tidyr::pivot_longer(d[, resid_cols, drop = FALSE],
      cols = everything(), names_to = "clock", values_to = "residual_age")
    long_resid$clock <- gsub("pheno_|phenoage_|_advance|_resid|_orig", "", long_resid$clock)
    plots$density <- ggplot(long_resid, aes(x = residual_age)) +
      geom_density(fill = "steelblue", alpha = 0.4, colour = "navy") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
      facet_wrap(~ clock, scales = "free") +
      theme_minimal(base_size = 12) +
      labs(title = "Distribution of Residual Biological Age",
           subtitle = "Values > 0 = biologically older than expected",
           x = "Residual Age (years)", y = "Density")
  }
}

# --- Plot 4: Correlation heatmap (ggplot) ---
if (!is.null(d)) {
  resid_cols <- list_resid_columns(d)
  if (length(resid_cols) >= 2) {
    mat <- cor_matrix(d, resid_cols)
    short <- gsub("pheno_|phenoage_|_advance|_resid_z|_resid|_orig", "", colnames(mat))
    rownames(mat) <- colnames(mat) <- short
    long_cor <- data.frame(
      Var1 = rep(rownames(mat), ncol(mat)),
      Var2 = rep(colnames(mat), each = nrow(mat)),
      value = as.vector(mat), stringsAsFactors = FALSE)
    plots$heatmap_resid <- ggplot(long_cor, aes(Var1, Var2, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3.5) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           midpoint = 0, limits = c(-1, 1)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank()) +
      labs(title = "Residual Sub-Clock Correlation", fill = "r",
           subtitle = "Low values = independent aging systems")
  }

  adv_cols <- list_advance_columns(d)
  if (length(adv_cols) >= 2) {
    mat_a <- cor_matrix(d, adv_cols)
    short_a <- gsub("pheno_|phenoage_|_advance|_orig", "", colnames(mat_a))
    rownames(mat_a) <- colnames(mat_a) <- short_a
    long_a <- data.frame(
      Var1 = rep(rownames(mat_a), ncol(mat_a)),
      Var2 = rep(colnames(mat_a), each = nrow(mat_a)),
      value = as.vector(mat_a), stringsAsFactors = FALSE)
    plots$heatmap_adv <- ggplot(long_a, aes(Var1, Var2, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3.5) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           midpoint = 0, limits = c(-1, 1)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank()) +
      labs(title = "BA Advancement Correlation", fill = "r",
           subtitle = "System concordance (BA - CA)")
  }
}

# --- Plot 5: KM curve ---
if (!is.null(d)) {
  has_mortality <- all(c("time", "status") %in% names(d)) &&
    any(!is.na(d$status) & d$status == 1)
  if (has_mortality) {
    global_resid <- intersect(
      c("phenoage_orig_advance_resid", "phenoage_global_advance_resid"), names(d))
    if (length(global_resid) > 0) {
      if (!"risk_group" %in% names(d))
        d <- add_risk_groups(d, global_resid[1])

      # Store KM as a base-R plot function (can't be a ggplot)
      plots$km <- "base_km"
    }
  }
}

# --- Plot 6: Advancement boxplots by age decade ---
if (!is.null(d)) {
  adv_cols <- list_advance_columns(d)
  if (length(adv_cols) > 0 && "age" %in% names(d)) {
    d$age_decade <- cut(d$age, breaks = seq(20, 90, 10),
                        labels = paste0(seq(20, 80, 10), "s"),
                        include.lowest = TRUE, right = FALSE)
    long_adv <- tryCatch({
      tidyr::pivot_longer(d[, c("age_decade", adv_cols), drop = FALSE],
        cols = all_of(adv_cols), names_to = "clock", values_to = "advancement")
    }, error = function(e) NULL)
    if (!is.null(long_adv)) {
      long_adv$clock <- gsub("pheno_|phenoage_|_advance|_orig", "", long_adv$clock)
      plots$boxplot_adv <- ggplot(long_adv, aes(x = age_decade, y = advancement)) +
        geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
        facet_wrap(~ clock, scales = "free_y") +
        theme_minimal(base_size = 12) +
        labs(title = "Biological Age Advancement by Age Decade",
             subtitle = "Values > 0 = biologically older than chronological age",
             x = "Age Decade", y = "Advancement (years)")
    }
  }
}

# ======================================================================
# SHOW ALL PLOTS IN VIEWER (before opening PDF)
# ======================================================================

message(">> Displaying plots in viewer...")
for (nm in names(plots)) {
  if (nm == "km") next  # base plot, handle separately
  tryCatch({
    print(plots[[nm]])
  }, error = function(e) message("  ! Plot '", nm, "' viewer error: ", e$message))
}

# KM in viewer
if ("km" %in% names(plots) && "risk_group" %in% names(d)) {
  tryCatch({
    plot_km(d, title = "Survival by Residual Biological Age")
  }, error = function(e) NULL)
}

# ======================================================================
# GENERATE PDF
# ======================================================================

pdf_path <- file.path(REPORT_ROOT, "bioage_full_report.pdf")
message(">> Writing PDF: ", pdf_path)

pdf(pdf_path, width = 11, height = 8.5)

# ---- Title page ----
draw_title_page()

# ---- Section 1: Panel Definitions ----
draw_section(1, "Biomarker Panel Definitions",
             "Which biomarkers go into each clock")

panel_lines <- c(
  "== PhenoAge Original (12 biomarkers) ==",
  paste(" ", BIOMARKERS_PHENO_ORIG),
  "",
  "== KDM / HD Panel ==",
  paste(" ", BIOMARKERS_KDM_HD),
  ""
)
for (nm in names(PANELS_SUB)) {
  panel_lines <- c(panel_lines,
    paste0("== Sub-clock: ", nm, " =="),
    paste(" ", PANELS_SUB[[nm]]), "")
}
draw_text_page("Biomarker Panels", panel_lines)

# ---- Section 2: Beta Coefficients ----
draw_section(2, "Beta Coefficients",
             "Trained model weights for each clock")

if (!is.null(betas_df) && nrow(betas_df) > 0) {
  for (clk in unique(betas_df$clock)) {
    sub <- betas_df[betas_df$clock == clk, c("variable", "beta"), drop = FALSE]
    sub$beta <- sprintf("%.6g", sub$beta)
    draw_table(sub, title = paste("Coefficients:", clk))
  }
} else {
  draw_text_page("Beta Coefficients",
    c("Coefficients stored in deployment bundle.",
      "Load with: readRDS('bioage_deployment_bundle.rds')"))
}

# ---- Section 3: Equations ----
draw_section(3, "Scoring Equations",
             "How to compute biological age from biomarkers")

eq_lines <- c(
  "== PHENOAGE (Levine 2018) ==", "",
  "Step 1: Linear predictor",
  "  xb = b0 + b1*albumin_gL + b2*alp + b3*lncrp + ...",
  "",
  "Step 2: Mortality probability (Gompertz)",
  "  M = 1 - exp(-exp(xb) * (exp(120*g)-1)/g)",
  "  where g = 0.090165",
  "",
  "Step 3: Invert to PhenoAge",
  "  PhenoAge = 141.50225 + ln(-0.00553*ln(1-M)) / 0.090165",
  "",
  "== KDM BIOLOGICAL AGE ==",
  "  KDM = sum(wi * age_i) / sum(wi)",
  "  age_i = (xi - intercept_i) / slope_i",
  "  wi = (slope_i / residual_SD_i)^2",
  "",
  "== HOMEOSTATIC DYSREGULATION (HD) ==",
  "  HD = sqrt((x-mu)' * Sigma^-1 * (x-mu))",
  "  Reference: young adults 20-30y",
  "",
  "== ADVANCEMENT ==",
  "  Advancement = BiologicalAge - ChronologicalAge",
  "",
  "== RESIDUALIZATION ==",
  "  Residual = residuals(lm(Advancement ~ Age))",
  "",
  "== Z-SCORING ==",
  "  Z = (Residual - mean) / SD",
  "  Z > 0 = biologically older",
  "",
  "== COX MODEL ==",
  "  h(t) = h0(t) * exp(b_age*Age + b_clock*Z)",
  "  HR = exp(b_clock) per +1 SD"
)
draw_text_page("Scoring Equations", eq_lines)

# ---- Section 4: QC ----
draw_section(4, "Quality Control", "Diagnostics for trained clocks")

qc_data <- NULL
if (!is.null(pipeline_results$qc) && is.data.frame(pipeline_results$qc))
  qc_data <- pipeline_results$qc
if (is.null(qc_data) && !is.null(bundle$qc) && is.data.frame(bundle$qc))
  qc_data <- bundle$qc
if (!is.null(qc_data))
  draw_table(qc_data, "QC Summary", "n_finite > 2000 and sd > 0 for valid clock")

# ---- Section 5: Population Stats ----
draw_section(5, "Population Reference Statistics",
             "Mean and SD for Z-scoring new patients")

if (!is.null(bundle$pop_stats))
  draw_table(bundle$pop_stats, "Population Reference (NHANES IV)",
             "Use these to Z-score new patient data")

# ---- Section 6: Residualization ----
draw_section(6, "Residualization Coefficients",
             "lm(advancement ~ age): intercept + slope")

if (!is.null(bundle$resid_coefs))
  draw_table(bundle$resid_coefs, "Age Regression Coefficients",
             "residual = advancement - (intercept + slope * age)")

# ---- Section 7: Hazard Ratios ----
draw_section(7, "Survival Analysis",
             "Cox PH models: Hazard Ratios per +1 SD")

if (!is.null(hr_data) && is.data.frame(hr_data) && nrow(hr_data) > 0) {
  draw_table(hr_data, "Hazard Ratio Table",
             "HR > 1 = increased mortality risk per +1 SD")
  # Forest plots
  for (nm in grep("^forest", names(plots), value = TRUE)) {
    tryCatch(print(plots[[nm]]), error = function(e) NULL)
  }
} else {
  draw_text_page("Hazard Ratios", c("No HR data available."))
}

# ---- Section 8: Plots ----
draw_section(8, "Diagnostic Plots",
             "Visual analysis of all clocks")

# BA vs Age scatter
if (!is.null(plots$ba_vs_age)) {
  tryCatch(print(plots$ba_vs_age), error = function(e) NULL)
}

# Heatmaps
if (!is.null(plots$heatmap_resid)) {
  tryCatch(print(plots$heatmap_resid), error = function(e) NULL)
}
if (!is.null(plots$heatmap_adv)) {
  tryCatch(print(plots$heatmap_adv), error = function(e) NULL)
}

# corrplot heatmaps (base R)
if (has_corrplot && !is.null(d)) {
  resid_cols <- list_resid_columns(d)
  if (length(resid_cols) >= 2) {
    tryCatch({
      mat <- cor_matrix(d, resid_cols)
      plot_heatmap_corrplot(mat, title = "Residual Sub-Clock Independence")
    }, error = function(e) NULL)
  }
  adv_cols <- list_advance_columns(d)
  if (length(adv_cols) >= 2) {
    tryCatch({
      mat_a <- cor_matrix(d, adv_cols)
      plot_heatmap_corrplot(mat_a, title = "BA Advancement Concordance")
    }, error = function(e) NULL)
  }
}

# Density
if (!is.null(plots$density)) {
  tryCatch(print(plots$density), error = function(e) NULL)
}

# Boxplot by decade
if (!is.null(plots$boxplot_adv)) {
  tryCatch(print(plots$boxplot_adv), error = function(e) NULL)
}

# KM curve
if ("km" %in% names(plots) && "risk_group" %in% names(d)) {
  tryCatch(plot_km(d, title = "Survival by Residual Biological Age"),
           error = function(e) NULL)
}

# ---- Section 9: Correlation Tables ----
draw_section(9, "Correlation Tables",
             "Numeric correlation data")

if (!is.null(d)) {
  tryCatch({
    ba_vars <- list_ba_columns(d)
    if (length(ba_vars) > 0) {
      cor_age_df <- cor_ba_vs_age(d, ba_vars)
      draw_table(cor_age_df, "BA vs Chronological Age", "Pearson r")
    }
  }, error = function(e) NULL)

  tryCatch({
    cor_gs <- cor_global_vs_subs(d)
    if (!is.null(cor_gs) && nrow(cor_gs) > 0)
      draw_table(cor_gs, "Global PhenoAge vs Sub-Clocks")
  }, error = function(e) NULL)

  tryCatch({
    disc <- discordance_analysis(d)
    if (!is.null(disc) && nrow(disc) > 0)
      draw_table(disc, "Pairwise System Discordance",
                 "Mean Z-score difference between systems")
  }, error = function(e) NULL)
}

# ---- Section 10: Clinical Reference ----
draw_section(10, "Clinical Reference Card",
             "How to interpret and apply the clocks")

ref <- c(
  "== HOW TO SCORE A NEW PATIENT ==", "",
  "1. Collect: CBC, BMP, LFTs, HbA1c, CRP, Uric Acid",
  "2. xb = sum(beta_i * biomarker_i)",
  "3. M = 1 - exp(-exp(xb)*(exp(120*0.090165)-1)/0.090165)",
  "4. PhenoAge = 141.50225 + ln(-0.00553*ln(1-M))/0.090165",
  "5. Advancement = PhenoAge - Age",
  "6. Residual = adv - (intercept + slope*age)  [Sec 6]",
  "7. Z = (residual - mean) / SD               [Sec 5]",
  "",
  "== Z-SCORE INTERPRETATION ==", "",
  "  Z <= -1.0   Biologically younger (low risk)",
  "  -1 < Z < 1  Average biological aging",
  "  Z >= 1.0    Biologically older (elevated risk)",
  "",
  "== HAZARD RATIO INTERPRETATION ==", "",
  "  HR = 1.00   No mortality association",
  "  HR = 1.10   10% higher risk per +1 SD",
  "  HR = 1.20   20% higher risk per +1 SD",
  "",
  "== SUB-CLOCK SYSTEMS ==", "",
  "  hepatic_enzime_insulin   Liver / insulin resistance",
  "  hepatic_lipid            Liver lipid metabolism",
  "  hema_integrated          Hematologic / bone marrow",
  "  micronutrient_methylation B12 / methylation",
  "  renal_A                  Kidney filtration",
  "",
  "== DEPLOYMENT ==", "",
  "  bundle <- readRDS('bioage_deployment_bundle.rds')",
  "  source('R/clinical_score.R')",
  "  results <- score_patients(patient_df, bundle)"
)
draw_text_page("Clinical Reference Card", ref)

# ---- Close PDF ----
dev.off()

cat("\n")
message("========================================")
message("  PDF REPORT COMPLETE")
message("  File: ", pdf_path)
message("  Plots shown in viewer: ", length(plots))
message("========================================")
