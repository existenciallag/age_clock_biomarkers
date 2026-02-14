###############################################################################
# correlation.R — Correlation matrices and heatmaps between clocks
#
# Three types of correlation analysis:
#   1. BA vs chronological age  — calibration check
#   2. Between residualized sub-clocks — independence / shared biology
#   3. Between advancement scores — system-level concordance
#
# Heatmaps use corrplot when available; ggplot2 tile fallback otherwise.
###############################################################################

# ── Correlation: BA vs chronological age ──────────────────────────────────

#' Pearson correlation of each BA column with chronological age
#'
#' @param data    Assembled data.frame.
#' @param ba_vars Character vector of BA column names.
#' @return A data.frame with columns: clock, cor_age.
cor_ba_vs_age <- function(data, ba_vars = NULL) {
  if (is.null(ba_vars)) ba_vars <- list_ba_columns(data)
  ba_vars <- intersect(ba_vars, names(data))

  out <- data.frame(
    clock   = ba_vars,
    cor_age = vapply(ba_vars, function(v) {
      cor(data[[v]], data$age, use = "pairwise.complete.obs")
    }, numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out[order(abs(out$cor_age), decreasing = TRUE), ]
}

# ── Correlation: clock-vs-clock matrix ────────────────────────────────────

#' Pairwise Pearson correlation matrix between a set of columns
#'
#' @param data data.frame.
#' @param vars Character vector of column names.
#' @return A numeric correlation matrix.
cor_matrix <- function(data, vars) {
  vars <- intersect(vars, names(data))
  cor(data[, vars, drop = FALSE], use = "pairwise.complete.obs")
}

# ── Correlation: global vs sub-clocks ─────────────────────────────────────

#' Correlation of PhenoAge global (or any reference) with each sub-clock
#'
#' @param data      data.frame.
#' @param ref_col   Name of the reference clock column.
#' @param sub_vars  Character vector of sub-clock columns.
#' @return A data.frame: subsystem, cor_ref.
cor_global_vs_subs <- function(data,
                               ref_col  = "phenoage_global",
                               sub_vars = NULL) {
  if (is.null(sub_vars)) {
    sub_vars <- grep("^pheno_[a-z]", names(data), value = TRUE)
    sub_vars <- sub_vars[!grepl("_advance|_resid|_z", sub_vars)]
  }
  sub_vars <- intersect(sub_vars, names(data))

  out <- data.frame(
    subsystem = sub_vars,
    cor_ref   = vapply(sub_vars, function(v) {
      cor(data[[ref_col]], data[[v]], use = "pairwise.complete.obs")
    }, numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out[order(abs(out$cor_ref), decreasing = TRUE), ]
}

# ── Heatmap: corrplot version ─────────────────────────────────────────────

#' Draw a correlation heatmap using corrplot
#'
#' @param mat      Correlation matrix (from cor_matrix()).
#' @param title    Plot title.
#' @param labels   Optional short labels for rows/columns.
#' @param palette  Colour ramp vector (low, mid, high).
plot_heatmap_corrplot <- function(mat,
                                  title   = "Clock Correlation",
                                  labels  = NULL,
                                  palette = HEATMAP_PALETTE) {
  if (!requireNamespace("corrplot", quietly = TRUE)) {
    message("corrplot not available; falling back to ggplot heatmap")
    return(plot_heatmap_gg(mat, title))
  }

  if (!is.null(labels)) {
    rownames(mat) <- colnames(mat) <- labels
  } else {
    short <- gsub("pheno_|phenoage_|_advance|_resid_z|_resid|_orig", "",
                  colnames(mat))
    rownames(mat) <- colnames(mat) <- short
  }

  corrplot::corrplot(
    mat,
    method       = "color",
    type         = "upper",
    tl.col       = "black",
    tl.srt       = 45,
    addCoef.col  = "white",
    number.cex   = 0.6,
    tl.cex       = 0.7,
    col          = grDevices::colorRampPalette(palette)(200),
    mar          = c(0, 0, 2, 0)
  )
  title(title, line = 1)
}

# ── Heatmap: ggplot2 fallback ─────────────────────────────────────────────

#' Draw a correlation heatmap using ggplot2
#'
#' @param mat   Correlation matrix.
#' @param title Plot title.
#' @return A ggplot object.
plot_heatmap_gg <- function(mat, title = "Clock Correlation") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Neither corrplot nor ggplot2 available for heatmap")
  }

  short <- gsub("pheno_|phenoage_|_advance|_resid_z|_resid|_orig", "",
                colnames(mat))
  rownames(mat) <- colnames(mat) <- short

  # Melt to long form
  long <- data.frame(
    Var1  = rep(rownames(mat), ncol(mat)),
    Var2  = rep(colnames(mat), each = nrow(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(long, ggplot2::aes(Var1, Var2, fill = value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "firebrick",
      midpoint = 0, limits = c(-1, 1)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title  = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = title, fill = "r")
}

# ── Discordance analysis ──────────────────────────────────────────────────

#' Compute pairwise system discordance from Z-scored sub-clocks
#'
#' For each pair of sub-clock systems, computes the per-subject difference
#' in Z-scored BA.  Large absolute differences indicate discordant aging.
#'
#' @param data      data.frame.
#' @param sub_vars  Character vector of sub-clock BA column names.
#' @return A data.frame with pairwise difference statistics.
discordance_analysis <- function(data, sub_vars = NULL) {
  if (is.null(sub_vars)) {
    sub_vars <- grep("^pheno_[a-z]", names(data), value = TRUE)
    sub_vars <- sub_vars[!grepl("_advance|_resid|_z", sub_vars)]
  }
  sub_vars <- intersect(sub_vars, names(data))
  if (length(sub_vars) < 2) {
    message(">> Need at least 2 sub-clocks for discordance analysis")
    return(NULL)
  }

  # Z-score each sub-clock
  z_mat <- as.data.frame(lapply(sub_vars, function(v) zscore(data[[v]])))
  names(z_mat) <- sub_vars

  pairs <- combn(sub_vars, 2, simplify = FALSE)

  out <- lapply(pairs, function(p) {
    diff <- z_mat[[p[1]]] - z_mat[[p[2]]]
    data.frame(
      clock_a   = p[1],
      clock_b   = p[2],
      mean_diff = mean(diff, na.rm = TRUE),
      sd_diff   = sd(diff, na.rm = TRUE),
      median_diff = median(diff, na.rm = TRUE),
      n         = sum(!is.na(diff)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}
