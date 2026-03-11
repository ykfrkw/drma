#' Plot a dose-response curve
#'
#' Draws the estimated dose-response curve with optional shaded 95% confidence
#' band, bubble overlay (study-level data points), and rug marks for tested
#' doses.  Returns a `ggplot` object that can be further modified with `+`.
#'
#' @param x          A `drma` object.
#' @param ref_dose   Reference dose for the y-axis (default: minimum observed
#'   dose).  Any value in the dose range is accepted; bubbles are automatically
#'   re-centred so they align with the shifted curve.
#' @param xlim       x-axis limits.  Default: `c(0, max(dose))`.
#' @param ylim       y-axis limits in the **display scale**.
#'   For OR/RR use the exponentiated scale (e.g. `c(0.75, 2)`, not
#'   `c(log(0.75), log(2))`).  Auto-computed if `NULL`.
#' @param xlab       x-axis label (default `"Dose"`).
#' @param ylab       y-axis label (default: inferred from `sm`).
#' @param n_pred     Number of prediction points for the curve (default 300).
#' @param ci         Logical; draw the 95% confidence band (default `TRUE`).
#' @param ci_col     Colour for the confidence band fill (default `"grey50"`).
#' @param ci_alpha   Transparency of the CI fill (0-1, default `0.25`).
#' @param add_abline Logical; draw a reference line at 1 (ratios) or 0
#'   (differences) (default `TRUE`).
#' @param bubble     Logical; overlay observed arm-level effects as bubbles,
#'   with area proportional to arm sample size (default `FALSE`).
#' @param bubble_scale `max_size` passed to [ggplot2::scale_size_area()]
#'   (default `5`).
#' @param bubble_col   Fill colour for bubbles (default `"steelblue"`).
#' @param bubble_alpha Transparency for bubbles (0-1, default `0.5`).
#' @param rug        Logical; draw rug tick marks along the x-axis for each
#'   evaluated dose (default `FALSE`).
#' @param rug_col    Colour for rug marks (default `"grey20"`).
#' @param col        Curve colour (default `"black"`).
#' @param lty        Line type — same codes as base R: 1=solid, 2=dashed,
#'   3=dotted, etc. (default `1`).
#' @param lwd        Line width in mm (default `0.5`).
#' @param ...        Currently unused; reserved for future parameters.
#'
#' @return A `ggplot` object.  Sensitivity-analysis overlays can be added with
#'   `+` and [lines.drma()]:
#'   ```r
#'   plot(res_primary, col = "black") +
#'     lines(res_s1, col = "tomato",    lty = 2) +
#'     lines(res_s2, col = "steelblue", lty = 3)
#'   ```
#'
#' @examples
#' \dontrun{
#' # Basic plot
#' plot(res, ylab = "Response (OR)", ylim = c(0.5, 3))
#'
#' # With bubbles and rug marks
#' plot(res, bubble = TRUE, rug = TRUE, ref_dose = 0)
#'
#' # Non-placebo reference (compare all doses to 1 mg)
#' plot(res, ref_dose = 1, bubble = TRUE)
#'
#' # Overlay sensitivity analyses
#' plot(res_primary, col = "black", ylim = c(0.5, 3)) +
#'   lines(res_s1, col = "tomato",    lty = 2) +
#'   lines(res_s2, col = "steelblue", lty = 3)
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
  bubble_scale = 5,
  bubble_col   = "steelblue",
  bubble_alpha = 0.5,
  rug          = FALSE,
  rug_col      = "grey20",
  col          = "black",
  lty          = 1,
  lwd          = 0.5,
  ...
) {
  is_ratio <- if (!is.null(x$is_ratio)) x$is_ratio else x$sm %in% c("OR", "RR")
  if (is.null(ylab))     ylab     <- .effect_label(x$sm)
  if (is.null(ref_dose)) ref_dose <- min(x$data$.dose, na.rm = TRUE)
  if (is.null(xlim))     xlim     <- c(0, max(x$data$.dose, na.rm = TRUE))

  # Prediction grid — ref_dose must be in grid for predict.dosresmeta(xref=)
  dose_seq <- seq(xlim[1], xlim[2], length.out = n_pred)
  if (!any(dose_seq == ref_dose)) dose_seq <- sort(c(dose_seq, ref_dose))
  nd   <- data.frame(.dose = dose_seq)
  pred <- predict(x$model, nd, xref = ref_dose, exp = is_ratio)
  df_pred <- data.frame(
    dose  = nd$.dose,
    pred  = pred$pred,
    ci_lb = pred$ci.lb,
    ci_ub = pred$ci.ub
  )

  p <- ggplot2::ggplot(df_pred, ggplot2::aes(x = .data$dose, y = .data$pred))

  # Reference line
  if (add_abline)
    p <- p + ggplot2::geom_hline(
      yintercept = if (is_ratio) 1 else 0,
      linetype = "dotted", color = "grey60"
    )

  # CI ribbon
  if (ci)
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lb, ymax = .data$ci_ub),
      fill = ci_col, alpha = ci_alpha
    )

  # Curve
  p <- p + ggplot2::geom_line(color = col, linetype = lty, linewidth = lwd)

  # Bubbles
  if (bubble) {
    d_bub <- x$data_fit[x$data_fit$.yi != 0, ]

    # Re-centre bubbles when ref_dose differs from the fitting reference
    # (which is always the minimum dose per study).
    fitting_ref <- min(x$data$.dose, na.rm = TRUE)
    adj <- if (!isTRUE(all.equal(ref_dose, fitting_ref))) {
      predict(x$model,
              data.frame(.dose = ref_dose),
              xref = fitting_ref, exp = FALSE)$pred
    } else {
      0
    }
    y_obs <- if (is_ratio) exp(d_bub$.yi - adj) else d_bub$.yi - adj

    df_bub <- data.frame(
      dose = d_bub$.dose,
      y    = y_obs,
      n    = if (!is.null(d_bub$.n)) d_bub$.n else rep(1, nrow(d_bub))
    )
    p <- p +
      ggplot2::geom_point(
        data        = df_bub,
        ggplot2::aes(x = .data$dose, y = .data$y, size = .data$n),
        color       = bubble_col,
        fill        = bubble_col,
        alpha       = bubble_alpha,
        shape       = 21,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_size_area(max_size = bubble_scale, guide = "none")
  }

  # Rug
  if (rug)
    p <- p + ggplot2::geom_rug(
      data        = data.frame(dose = sort(unique(x$data$.dose))),
      ggplot2::aes(x = .data$dose),
      color       = rug_col,
      sides       = "b",
      inherit.aes = FALSE
    )

  # Axis labels
  p <- p + ggplot2::labs(x = xlab, y = ylab)

  # y-axis scale.  ylim is always in display units (OR/RR scale for ratios).
  # scale_y_log10(limits=) takes original-scale values directly — no log conversion needed.
  if (is_ratio) {
    p <- p + ggplot2::scale_y_log10(limits = ylim, oob = scales::squish)
  } else if (!is.null(ylim)) {
    p <- p + ggplot2::scale_y_continuous(limits = ylim, oob = scales::squish)
  }
  p <- p + ggplot2::coord_cartesian(xlim = xlim)

  p + ggplot2::theme_classic()
}


#' Add a dose-response curve to an existing ggplot
#'
#' Returns a list of [ggplot2] layers (geom_line + optional geom_ribbon) that
#' can be added to a plot produced by [plot.drma()] using `+`.
#'
#' @param x        A `drma` object.
#' @param ref_dose Reference dose (default: minimum observed dose in `x`).
#' @param n_pred   Number of prediction points (default 300).
#' @param ci       Logical; draw the confidence band (default `FALSE`).
#' @param ci_col   Fill colour for the CI band (default `"grey50"`).
#' @param ci_alpha Transparency of the CI fill (0-1, default `0.25`).
#' @param xlim     x range for the prediction grid.  If `NULL`, uses the
#'   observed dose range of `x`.
#' @param col      Line colour (default `"black"`).
#' @param lty      Line type (default `1`).
#' @param lwd      Line width in mm (default `0.5`).
#' @param ...      Currently unused.
#'
#' @return A list of ggplot2 layers.  Add to a `plot.drma()` result with `+`.
#'
#' @examples
#' \dontrun{
#' plot(res_primary, col = "black", ylim = c(0.5, 3)) +
#'   lines(res_s1, col = "red",  lty = 2, lwd = 0.8) +
#'   lines(res_s2, col = "blue", lty = 3, lwd = 0.8)
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
  col      = "black",
  lty      = 1,
  lwd      = 0.5,
  ...
) {
  is_ratio <- if (!is.null(x$is_ratio)) x$is_ratio else x$sm %in% c("OR", "RR")
  if (is.null(ref_dose)) ref_dose <- min(x$data$.dose, na.rm = TRUE)
  if (is.null(xlim))     xlim     <- c(0, max(x$data$.dose, na.rm = TRUE))

  dose_seq <- seq(xlim[1], xlim[2], length.out = n_pred)
  if (!any(dose_seq == ref_dose)) dose_seq <- sort(c(dose_seq, ref_dose))
  nd   <- data.frame(.dose = dose_seq)
  pred <- predict(x$model, nd, xref = ref_dose, exp = is_ratio)
  df_pred <- data.frame(
    dose  = nd$.dose,
    pred  = pred$pred,
    ci_lb = pred$ci.lb,
    ci_ub = pred$ci.ub
  )

  layers <- list()
  if (ci)
    layers[[length(layers) + 1]] <- ggplot2::geom_ribbon(
      data        = df_pred,
      ggplot2::aes(x = .data$dose, ymin = .data$ci_lb, ymax = .data$ci_ub),
      fill        = ci_col,
      alpha       = ci_alpha,
      inherit.aes = FALSE
    )
  layers[[length(layers) + 1]] <- ggplot2::geom_line(
    data        = df_pred,
    ggplot2::aes(x = .data$dose, y = .data$pred),
    color       = col,
    linetype    = lty,
    linewidth   = lwd,
    inherit.aes = FALSE
  )
  layers
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
