###############################################################################
# assemble.R — Build a unified dataset from trained clock objects
#
# Takes the output of train_all_clocks() and merges every clock's NHANES IV
# projection into a single data.frame keyed by sampleID.
#
# Columns produced:
#   sampleID, age, gender                — demographics
#   time, status                         — mortality follow-up (auto-detected)
#   phenoage_orig, phenoage_global       — canonical PhenoAge BA
#   kdm, hd, hd_log                      — KDM & HD BA
#   pheno_<name>                          — sub-clock BA values
#   *_advance                             — BA − CA (biological age advancement)
###############################################################################

# ── Auto-detect mortality columns ──────────────────────────────────────────

detect_column <- function(df, patterns, label) {
  for (pat in patterns) {
    hit <- grep(pat, names(df), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  message(">> Warning: could not detect '", label, "' column")
  NA_character_
}

# ── Main assembly function ─────────────────────────────────────────────────

#' Assemble a unified biological-age dataset
#'
#' @param clocks  The list returned by train_all_clocks().
#' @return A data.frame with one row per NHANES IV subject, containing
#'         all clock BA values, advancement scores, and mortality columns.
assemble_data <- function(clocks) {

  # --- Start from PhenoAge original (guaranteed to exist) ------------------
  base <- clocks$pheno_orig$data

  # Auto-detect mortality / follow-up columns
  mort_col <- detect_column(base, MORT_STATUS_PATTERNS, "mortality status")
  time_col <- detect_column(base, MORT_TIME_PATTERNS,   "follow-up time")

  # Core columns always present
  keep <- c("sampleID", "age", "gender")
  if (!is.na(mort_col)) keep <- c(keep, mort_col)
  if (!is.na(time_col)) keep <- c(keep, time_col)
  keep <- intersect(keep, names(base))

  data <- base[, c(keep, "phenoage", "phenoage_advance"), drop = FALSE]
  names(data)[names(data) == "phenoage"]         <- "phenoage_orig"
  names(data)[names(data) == "phenoage_advance"] <- "phenoage_orig_advance"

  # Standardise mortality column names for downstream code
  if (!is.na(mort_col) && mort_col != "status") {
    data$status <- data[[mort_col]]
  }
  if (!is.na(time_col) && time_col != "time") {
    data$time <- data[[time_col]]
  }

  # --- PhenoAge global -----------------------------------------------------
  if (!is.null(clocks$pheno_global)) {
    tmp <- clocks$pheno_global$data[, c("sampleID", "phenoage"), drop = FALSE]
    names(tmp)[2] <- "phenoage_global"
    data <- merge(data, tmp, by = "sampleID", all.x = TRUE)
    data$phenoage_global_advance <- data$phenoage_global - data$age
  }

  # --- KDM -----------------------------------------------------------------
  if (!is.null(clocks$kdm)) {
    cols <- intersect(c("sampleID", "kdm", "kdm_advance"), names(clocks$kdm$data))
    tmp  <- clocks$kdm$data[, cols, drop = FALSE]
    data <- merge(data, tmp, by = "sampleID", all.x = TRUE)
  }

  # --- HD ------------------------------------------------------------------
  if (!is.null(clocks$hd)) {
    cols <- intersect(c("sampleID", "hd", "hd_log"), names(clocks$hd$data))
    tmp  <- clocks$hd$data[, cols, drop = FALSE]
    data <- merge(data, tmp, by = "sampleID", all.x = TRUE)
    if (!"hd_log" %in% names(data) && "hd" %in% names(data)) {
      data$hd_log <- log(data$hd)
    }
  }

  # --- Sub-clocks ----------------------------------------------------------
  for (nm in names(clocks$subclocks)) {
    col_ba  <- paste0("pheno_", nm)
    col_adv <- paste0("pheno_", nm, "_advance")

    tmp <- clocks$subclocks[[nm]]$data[, c("sampleID", "phenoage"), drop = FALSE]
    names(tmp)[2] <- col_ba
    data <- merge(data, tmp, by = "sampleID", all.x = TRUE)
    data[[col_adv]] <- data[[col_ba]] - data$age
  }

  # --- Advancement for canonical clocks (if not already present) -----------
  if (!"kdm_advance" %in% names(data) && "kdm" %in% names(data)) {
    data$kdm_advance <- data$kdm - data$age
  }

  message(">> Assembly complete: ",
          nrow(data), " subjects, ",
          ncol(data), " columns")

  data
}

# ── Helper: list all clock column names in the assembled dataset ───────────

#' Return the names of all BA / advancement / HD columns
#' @param data The assembled data.frame.
#' @return Character vector of clock-related column names.
list_clock_columns <- function(data) {
  grep("^pheno_|^phenoage_|^kdm|^hd", names(data), value = TRUE)
}

#' Return only the "raw BA" columns (no advancement / residuals)
list_ba_columns <- function(data) {
  ba <- c("phenoage_orig", "phenoage_global", "kdm", "hd", "hd_log")
  sub <- grep("^pheno_[a-z]", names(data), value = TRUE)
  sub <- sub[!grepl("_advance$|_resid|_z$", sub)]
  intersect(c(ba, sub), names(data))
}

#' Return only advancement columns
list_advance_columns <- function(data) {
  grep("_advance$", names(data), value = TRUE)
}
