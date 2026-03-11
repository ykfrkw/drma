#' Plot a dose-response curve
#'
#' Draws the estimated dose-response curve with optional 95% confidence band,
#' bubble overlay (study-level data points), and rug marks for tested doses.
#'
#' @param x          A `drma` object.
#' @param ref_dose   Reference dose for the y-axis (default: minimum observed
#'   dose).
#' @param xlim       x-axis limits.  Default: `c(0, max(dose))`.
#' @param ylim       y-axis limits.  Auto-computed if `NULL`.
#' @param xlab       x-axis label (default `"Dose"`).
#' @param ylab       y-axis label (default: inferred from `sm`).
#' @param n_pred     Number of prediction points for the curve (default 300).
#' @param ci         Logical; draw the 95% confidence band (default `TRUE`).
#' @param ci_lty     Line type for the confidence band (default `"dashed"`).
#' @param ci_col     Colour for the confidence band (default `1`).
#' @param add_abline Logical; draw a reference line at 1 (ratios) or 0
#'   (differences) (default `TRUE`).
#' @param bubble     Logical; overlay observed arm-level log-effects as
#'   bubbles, with area proportional to arm sample size (default `FALSE`).
#' @param bubble_scale Numeric scale factor for bubble sizes.  The largest
#'   bubble has diameter `bubble_scale * 0.2` inches on the printed figure
#'   (default `1`, i.e. 0.2 inches max).
#' @param bubble_col   Fill colour for bubbles (default `"steelblue"`).
#' @param bubble_alpha Transparency for bubbles (0–1, default `0.5`).
#' @param rug        Logical; draw rug tick marks along the x-axis for each
#'   evaluated dose (default `FALSE`).
#' @param rug_col    Colour for rug marks (default `"grey20"`).
#' @param rug_ticksize Relative length of rug tick marks as a fraction of
#'   the plot height (default `0.05`).
#' @param add        Logical; add the curve to an existing plot without
#'   redrawing axes — useful for overlaying multiple models (default `FALSE`).
#' @param ...        Graphical parameters passed to [graphics::plot()] or
#'   [graphics::lines()].
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \dontrun{
#' # Basic plot
#' plot(res, ylab = "Response (OR)", ylim = c(0.5, 3))
#'
#' # With bubbles and rug marks
#' plot(res, bubble = TRUE, rug = TRUE, ref_dose = 0)
#'
#' # Overlay two models
#' plot(res_primary, col = "black", ylim = c(0.5, 3))
#' lines(res_s1, col = "red",  lty = 2)
#' lines(res_s2, col = "blue", lty = 3)
#' }
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
    rug_col      = "grey20",
    rug_ticksize = 0.05,
    add          = FALSE,
    ...
) {
  is_ratio <- x$sm %in% c("OR", "RR")

  if (is.null(ylab))     ylab     <- .effect_label(x$sm)
  if (is.null(ref_dose)) ref_dose <- min(x$data$.dose, na.rm = TRUE)
  if (is.null(xlim))     xlim     <- c(0, max(x$data$.dose, na.rm = TRUE))

  nd   <- data.frame(.dose = seq(xlim[1], xlim[2], length.out = n_pred))
  pred <- predict(x$model, nd, xref = ref_dose, exp = is_ratio)
  xs   <- nd$.dose   # use directly so all curve types work

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
    y_obs <- if (is_ratio) exp(d_bub$.yi) else d_bub$.yi

    # Bubble radius proportional to sqrt(n) so area ~ n.
    # sz is normalised to [0, 1]; the largest bubble = bubble_scale * 0.2 inches.
    sz <- if (!is.null(d_bub$.n)) {
      sqrt(d_bub$.n) / max(sqrt(d_bub$.n), na.rm = TRUE)
    } else {
      rep(1, nrow(d_bub))
    }
    col_rgb <- grDevices::col2rgb(bubble_col) / 255
    bub_col <- grDevices::rgb(col_rgb[1], col_rgb[2], col_rgb[3],
                               alpha = bubble_alpha)
    graphics::symbols(d_bub$.dose, y_obs,
                      circles  = sz,
                      inches   = bubble_scale * 0.2,
                      add      = TRUE,
                      bg       = bub_col,
                      fg       = bubble_col)
  }

  # ── Rug marks ─────────────────────────────────────────────────────────────
  if (rug) {
    tested <- sort(unique(x$data$.dose))
    graphics::rug(tested, col = rug_col, ticksize = rug_ticksize, lwd = 1.5)
  }

  invisible(x)
}


#' Add a dose-response curve to an existing plot
#'
#' Convenience function for overlaying multiple dose-response curves on a
#' single plot.  Equivalent to `plot.drma(..., add = TRUE)`.
#'
#' @param x        A `drma` object.
#' @param ref_dose Reference dose (default: minimum observed dose in `x`).
#' @param n_pred   Number of prediction points (default 300).
#' @param ci       Logical; draw the confidence band (default `FALSE`).
#' @param xlim     x range for the prediction grid.  If `NULL`, uses the
#'   observed dose range of `x`.
#' @param ...      Graphical parameters, e.g. `col`, `lty`, `lwd`.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \dontrun{
#' plot(res_primary, col = "black", ylim = c(0.5, 3))
#' lines(res_s1, col = "red",  lty = 2, lwd = 1.5)
#' lines(res_s2, col = "blue", lty = 3, lwd = 1.5)
#' legend("topright",
#'        legend = c("Primary", "Sensitivity 1", "Sensitivity 2"),
#'        col    = c("black", "red", "blue"),
#'        lty    = 1:3, lwd = 1.5, bty = "n")
#' }
#' @export
lines.drma <- function(
    x,
    ref_dose = NULL,
    n_pred   = 300,
    ci       = FALSE,
    xlim     = NULL,
    ...
) {
  is_ratio <- x$sm %in% c("OR", "RR")
  if (is.null(ref_dose)) ref_dose <- min(x$data$.dose, na.rm = TRUE)
  if (is.null(xlim))     xlim     <- c(0, max(x$data$.dose, na.rm = TRUE))

  nd   <- data.frame(.dose = seq(xlim[1], xlim[2], length.out = n_pred))
  pred <- predict(x$model, nd, xref = ref_dose, exp = is_ratio)
  xs   <- nd$.dose

  graphics::lines(xs, pred$pred, ...)
  if (ci)
    graphics::matlines(xs, cbind(pred$ci.ub, pred$ci.lb), lty = "dashed", ...)

  invisible(x)
}


#' VPC (Variance Partitioning Coefficient) plot
#'
#' Visualises how well the dose-response model accounts for the between-study
#' heterogeneity.
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
