#' Plot a dose-response curve
#'
#' Draws the estimated dose-response curve with optional shaded 95% confidence
#' band, bubble overlay (study-level data points), and rug marks for tested
#' doses.
#'
#' @param x          A `drma` object.
#' @param ref_dose   Reference dose for the y-axis (default: minimum observed
#'   dose).  Any value in the dose range is accepted; bubbles are automatically
#'   re-centred so they align with the shifted curve.
#' @param xlim       x-axis limits.  Default: `c(0, max(dose))`.
#' @param ylim       y-axis limits.  Auto-computed if `NULL`.
#' @param xlab       x-axis label (default `"Dose"`).
#' @param ylab       y-axis label (default: inferred from `sm`).
#' @param n_pred     Number of prediction points for the curve (default 300).
#' @param ci         Logical; draw the 95% confidence band as a shaded polygon
#'   (default `TRUE`).
#' @param ci_col     Colour for the confidence band fill (default `"grey50"`).
#' @param ci_alpha   Transparency of the CI fill (0–1, default `0.25`).
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
#' # Non-placebo reference (e.g. compare all doses to 1 mg)
#' plot(res, ref_dose = 1, bubble = TRUE)
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
  ci_col       = "grey50",
  ci_alpha     = 0.25,
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
  # is_ratio: TRUE for OR/RR and for precomputed log-ratio data.
  # Falls back to sm-based detection for drma objects built before is_ratio
  # was stored (backward compatibility).
  is_ratio <- if (!is.null(x$is_ratio)) x$is_ratio else x$sm %in% c("OR", "RR")

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
                   type = "n", xlim = xlim, ylim = ylim,
                   xlab = xlab, ylab = ylab,
                   log  = log_arg, bty  = "l", las  = 1, ...)
    if (add_abline)
      graphics::abline(h = if (is_ratio) 1 else 0, lty = 3, col = "grey60")
  }

  # ── CI shaded polygon ────────────────────────────────────────────────────
  if (ci) {
    ok <- is.finite(pred$ci.lb) & is.finite(pred$ci.ub)
    if (any(ok)) {
      col_rgb  <- grDevices::col2rgb(ci_col) / 255
      ci_fill  <- grDevices::rgb(col_rgb[1], col_rgb[2], col_rgb[3],
                               alpha = ci_alpha)
      poly_x   <- c(xs[ok], rev(xs[ok]))
      poly_y   <- c(pred$ci.ub[ok], rev(pred$ci.lb[ok]))
      graphics::polygon(poly_x, poly_y, col = ci_fill, border = NA)
    }
  }

  # Draw the curve on top of the CI band
  if (!add) {
    graphics::lines(xs, pred$pred, ...)
  } else {
    graphics::lines(xs, pred$pred, ...)
  }

  # ── Bubble plot ────────────────────────────────────────────────────────────
  if (bubble) {
    d_bub <- x$data_fit[x$data_fit$.yi != 0, ]  # exclude reference arms

    # Re-centre bubble y-values for the chosen ref_dose.
    # data_fit$.yi are contrasts vs. the fitting reference (min dose per study).
    # When ref_dose differs from the fitting reference, subtract the model's
    # predicted log-effect at ref_dose (relative to the fitting reference) so
    # bubbles align with the shifted curve.
    fitting_ref <- min(x$data$.dose, na.rm = TRUE)
    if (!isTRUE(all.equal(ref_dose, fitting_ref))) {
      adj <- predict(x$model,
                     data.frame(.dose = ref_dose),
                     xref = fitting_ref, exp = FALSE)$pred
    } else {
      adj <- 0
    }
    y_obs <- if (is_ratio) exp(d_bub$.yi - adj) else d_bub$.yi - adj  # nolint

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
#' @param ci       Logical; draw the confidence band as a shaded polygon
#'   (default `FALSE`).
#' @param ci_col   Fill colour for the CI band (default `"grey50"`).
#' @param ci_alpha Transparency of the CI fill (0–1, default `0.25`).
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
  ci_col   = "grey50",
  ci_alpha = 0.25,
  xlim     = NULL,
  ...
) {
  is_ratio <- if (!is.null(x$is_ratio)) x$is_ratio else x$sm %in% c("OR", "RR")
  if (is.null(ref_dose)) ref_dose <- min(x$data$.dose, na.rm = TRUE)
  if (is.null(xlim))     xlim     <- c(0, max(x$data$.dose, na.rm = TRUE))

  nd   <- data.frame(.dose = seq(xlim[1], xlim[2], length.out = n_pred))
  pred <- predict(x$model, nd, xref = ref_dose, exp = is_ratio)
  xs   <- nd$.dose

  if (ci) {
    ok <- is.finite(pred$ci.lb) & is.finite(pred$ci.ub)
    if (any(ok)) {
      col_rgb <- grDevices::col2rgb(ci_col) / 255
      ci_fill <- grDevices::rgb(col_rgb[1], col_rgb[2], col_rgb[3],
                                alpha = ci_alpha)
      graphics::polygon(c(xs[ok], rev(xs[ok])),
                        c(pred$ci.ub[ok], rev(pred$ci.lb[ok])),
                        col = ci_fill, border = NA)
    }
  }
  graphics::lines(xs, pred$pred, ...)

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
