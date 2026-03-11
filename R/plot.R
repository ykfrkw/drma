#' Plot a dose-response curve
#'
#' Draws the estimated dose-response curve with optional 95% confidence band,
#' bubble overlay (study-level data), and rug marks for tested doses.
#'
#' @param x         A `drma` object.
#' @param ref_dose  Reference dose for the y-axis (default: minimum observed
#'   dose).
#' @param xlim      x-axis limits.  Default: `c(0, max(dose))`.
#' @param ylim      y-axis limits.  Auto-computed if `NULL`.
#' @param xlab      x-axis label (default `"Dose"`).
#' @param ylab      y-axis label (default: inferred from `measure`).
#' @param n_pred    Number of prediction points (default 300).
#' @param ci        Logical; draw the 95% confidence band (default `TRUE`).
#' @param ci_lty    Line type for the confidence band (default `"dashed"`).
#' @param ci_col    Colour for the confidence band (default `1`).
#' @param add_abline Logical; draw a horizontal reference line at 1 (ratios)
#'   or 0 (differences) (default `TRUE`).
#' @param bubble    Logical; overlay a bubble plot showing the observed
#'   log-effect per arm, with bubble area proportional to arm sample size
#'   (default `FALSE`).
#' @param bubble_scale Numeric scale factor for bubble sizes (default `1`).
#' @param bubble_col   Colour for bubbles (default `"steelblue"`).
#' @param bubble_alpha Transparency for bubbles, passed as `col` alpha
#'   (0–1, default `0.5`).
#' @param rug        Logical; draw rug marks along the x-axis for each
#'   evaluated dose (default `FALSE`).
#' @param rug_col    Colour for rug marks (default `"grey40"`).
#' @param add        Logical; if `TRUE`, add the curve to an existing plot
#'   without drawing axes (for overlaying multiple models, default `FALSE`).
#' @param ...       Graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#' @export
plot.drma <- function(
    x,
    ref_dose     = NULL,
    xlim         = NULL,
    ylim         = NULL,
    xlab         = "Dose",
    ylab         = NULL,
    n_pred       = 300,
    ci           = TRUE,
    ci_lty       = "dashed",
    ci_col       = 1,
    add_abline   = TRUE,
    bubble       = FALSE,
    bubble_scale = 1,
    bubble_col   = "steelblue",
    bubble_alpha = 0.5,
    rug          = FALSE,
    rug_col      = "grey40",
    add          = FALSE,
    ...
) {
  is_ratio <- x$measure %in% c("OR", "RR")

  if (is.null(ylab))     ylab     <- .effect_label(x$measure)
  if (is.null(ref_dose)) ref_dose <- min(x$data$.dose, na.rm = TRUE)
  if (is.null(xlim))     xlim     <- c(0, max(x$data$.dose, na.rm = TRUE))

  nd   <- data.frame(.dose = seq(xlim[1], xlim[2], length.out = n_pred))
  pred <- predict(x$model, nd, xref = ref_dose, exp = is_ratio)
  xs   <- pred[[.rcs_col(pred)]]

  if (is.null(ylim)) {
    ys   <- c(pred$pred, pred$ci.lb, pred$ci.ub)
    ylim <- range(ys[is.finite(ys)], na.rm = TRUE)
    if (is_ratio) ylim[1] <- max(ylim[1], 0.01)
  }

  log_arg <- if (is_ratio) "y" else ""

  if (!add) {
    graphics::plot(xs, pred$pred,
                   type = "l", xlim = xlim, ylim = ylim,
                   xlab = xlab, ylab = ylab,
                   log  = log_arg, bty  = "l", las  = 1, ...)
    if (add_abline)
      graphics::abline(h = if (is_ratio) 1 else 0, lty = 3, col = "grey60")
  } else {
    graphics::lines(xs, pred$pred, ...)
  }

  if (ci)
    graphics::matlines(xs, cbind(pred$ci.ub, pred$ci.lb),
                       col = ci_col, lty = ci_lty)

  # ── Bubble plot ────────────────────────────────────────────────────────────
  if (bubble) {
    d_bub <- x$data_fit[x$data_fit$.yi != 0, ]  # exclude reference arms
    if (is_ratio) {
      y_bub <- exp(d_bub$.yi + pred$pred[1] - pred$pred[1])  # offset to ref_dose
      # Re-predict at each study dose relative to ref_dose
      nd_bub <- data.frame(.dose = d_bub$.dose)
      pr_bub <- predict(x$model, nd_bub, xref = ref_dose, exp = TRUE)
      y_obs  <- exp(d_bub$.yi)  # observed OR/RR vs own reference arm
    } else {
      y_obs <- d_bub$.yi
    }
    # Bubble size proportional to sqrt(n) for area ~ n
    sz <- if (!is.null(d_bub$.n)) {
      bubble_scale * sqrt(d_bub$.n) / max(sqrt(d_bub$.n), na.rm = TRUE) * 3
    } else {
      rep(bubble_scale, nrow(d_bub))
    }
    col_rgb <- grDevices::col2rgb(bubble_col) / 255
    bub_col <- grDevices::rgb(col_rgb[1], col_rgb[2], col_rgb[3],
                               alpha = bubble_alpha)
    graphics::symbols(d_bub$.dose, y_obs,
                      circles  = sz,
                      inches   = FALSE,
                      add      = TRUE,
                      bg       = bub_col,
                      fg       = bubble_col)
  }

  # ── Rug marks ─────────────────────────────────────────────────────────────
  if (rug) {
    tested <- unique(x$data$.dose)
    graphics::rug(tested, col = rug_col, ticksize = 0.03)
  }

  invisible(x)
}


#' Add a dose-response curve to an existing plot
#'
#' Convenience wrapper around `plot.drma(..., add = TRUE)` for overlaying
#' multiple dose-response curves.
#'
#' @param x       A `drma` object.
#' @param ref_dose Reference dose (default: minimum observed dose).
#' @param n_pred  Number of prediction points (default 300).
#' @param ci      Logical; draw the confidence band (default `FALSE`).
#' @param xlim    x range used for the prediction grid.  If `NULL`, uses the
#'   observed dose range of `x`.
#' @param ...     Graphical parameters (e.g. `col`, `lty`, `lwd`).
#'
#' @return Invisibly returns `x`.
#' @export
lines.drma <- function(
    x,
    ref_dose = NULL,
    n_pred   = 300,
    ci       = FALSE,
    xlim     = NULL,
    ...
) {
  is_ratio <- x$measure %in% c("OR", "RR")
  if (is.null(ref_dose)) ref_dose <- min(x$data$.dose, na.rm = TRUE)
  if (is.null(xlim))     xlim     <- c(0, max(x$data$.dose, na.rm = TRUE))

  nd   <- data.frame(.dose = seq(xlim[1], xlim[2], length.out = n_pred))
  pred <- predict(x$model, nd, xref = ref_dose, exp = is_ratio)
  xs   <- pred[[.rcs_col(pred)]]

  graphics::lines(xs, pred$pred, ...)
  if (ci)
    graphics::matlines(xs, cbind(pred$ci.ub, pred$ci.lb), lty = "dashed", ...)

  invisible(x)
}


#' VPC (Variance Partitioning Coefficient) plot
#'
#' @param x   A `drma` object.
#' @param ... Graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns the VPC values.
#' @export
vpc_plot <- function(x, ...) {
  vcps <- dosresmeta::vpc(x$model)
  d    <- x$data_fit
  graphics::plot(d$.dose, vcps,
                 xlab = "Dose", ylab = "VPC",
                 pch  = 16, col = "steelblue", ...)
  graphics::lines(stats::lowess(d$.dose, vcps), col = "red")
  invisible(vcps)
}
