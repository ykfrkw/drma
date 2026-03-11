#' Predictions from a drma model
#'
#' @param object  A `drma` object.
#' @param doses   Doses at which to predict.  If `NULL`, 100 evenly-spaced
#'   values from `ref_dose` to the observed maximum are used.
#' @param ref_dose Reference dose (default: minimum observed dose).
#' @param expo    Logical; exponentiate predictions?  Default `TRUE` for
#'   OR/RR, `FALSE` for MD/SMD.
#' @param ...     Additional arguments passed to
#'   [dosresmeta::predict.dosresmeta()].
#'
#' @return The output of [dosresmeta::predict.dosresmeta()].
#' @export
predict.drma <- function(
    object,
    doses    = NULL,
    ref_dose = NULL,
    expo     = NULL,
    ...
) {
  if (is.null(ref_dose)) ref_dose <- min(object$data$.dose, na.rm = TRUE)
  if (is.null(expo))     expo     <- object$measure %in% c("OR", "RR")
  if (is.null(doses))
    doses <- seq(ref_dose,
                 max(object$data$.dose, na.rm = TRUE),
                 length.out = 100)
  nd <- data.frame(.dose = doses)
  predict(object$model, newdata = nd, xref = ref_dose, exp = expo, ...)
}


#' Predicted effects at specified doses
#'
#' Returns a tidy table of predicted effects (with 95% CI) at user-specified
#' doses.  Optionally converts OR/RR to absolute risk using a baseline
#' proportion.
#'
#' @param x             A `drma` object.
#' @param doses         Numeric vector of doses.
#' @param ref_dose      Reference dose (default: minimum observed dose).
#' @param baseline_prop Baseline proportion used to convert OR/RR to absolute
#'   risk (e.g. `0.30` for 30%).  If `NULL`, absolute risk columns are omitted.
#'
#' @return A `data.frame` with columns `dose`, the effect estimate and 95% CI
#'   (column names depend on `measure`), and optionally `abs_risk`,
#'   `abs_risk.lb`, `abs_risk.ub` (%).
#' @export
predict_table <- function(x, doses, ref_dose = NULL, baseline_prop = NULL) {
  if (is.null(ref_dose)) ref_dose <- min(x$data$.dose, na.rm = TRUE)
  is_ratio <- x$measure %in% c("OR", "RR")

  nd   <- data.frame(.dose = doses)
  pred <- predict(x$model, newdata = nd, xref = ref_dose, exp = is_ratio)

  eff_names <- switch(x$measure,
    OR  = c("OR",  "OR.lb",  "OR.ub"),
    RR  = c("RR",  "RR.lb",  "RR.ub"),
    MD  = c("MD",  "MD.lb",  "MD.ub"),
    SMD = c("SMD", "SMD.lb", "SMD.ub"),
    c("effect", "ci.lb", "ci.ub")
  )

  res <- data.frame(
    dose  = doses,
    e     = round(pred$pred,  3),
    ci.lb = round(pred$ci.lb, 3),
    ci.ub = round(pred$ci.ub, 3)
  )
  names(res)[2:4] <- eff_names

  if (!is.null(baseline_prop) && is_ratio) {
    bl  <- baseline_prop / (1 - baseline_prop)
    eff <- res[[2]]; lb <- res[[3]]; ub <- res[[4]]
    res$abs_risk    <- round(100 * bl * eff / (1 + bl * eff), 1)
    res$abs_risk.lb <- round(100 * bl * lb  / (1 + bl * lb),  1)
    res$abs_risk.ub <- round(100 * bl * ub  / (1 + bl * ub),  1)
  }
  res
}


#' Estimate target doses (ED50, ED95, etc.)
#'
#' Finds the dose corresponding to a given fraction of the maximum predicted
#' effect (analogous to ED50, ED95, ED100 in pharmacology).
#'
#' @param x      A `drma` object.
#' @param p      Numeric vector of fractions of the maximum effect
#'   (e.g. `c(0.5, 0.95, 1)`).
#' @param n_pred Number of prediction points for the search grid (default
#'   1 000).
#' @param trunc  If `TRUE`, return `NA` when the plateau is not reached within
#'   the observed dose range (default `FALSE`).
#'
#' @return A `data.frame` with columns `p`, `dose`, `log_effect`, and
#'   `effect`.
#' @export
target_dose <- function(x, p = c(0.5, 0.95, 1), n_pred = 1000, trunc = FALSE) {
  dose_range <- range(x$data$.dose, na.rm = TRUE)
  doses      <- seq(dose_range[1], dose_range[2], length.out = n_pred)
  nd         <- data.frame(.dose = doses)
  pred_log   <- predict(x$model, nd)$pred  # log scale

  max_eff <- max(pred_log, na.rm = TRUE)
  ed_max  <- doses[which.max(pred_log)]

  if (trunc && ed_max >= max(doses) * 0.99) {
    return(data.frame(p = p, dose = NA_real_,
                      log_effect = NA_real_, effect = NA_real_))
  }

  ED <- vapply(p, function(pi) {
    target <- pi * max_eff
    idx    <- doses < ed_max
    if (!any(idx)) return(NA_real_)
    doses[idx][which.min(abs(pred_log[idx] - target))]
  }, numeric(1))

  is_ratio <- x$measure %in% c("OR", "RR")
  data.frame(
    p          = p,
    dose       = round(ED, 4),
    log_effect = round(p * max_eff, 4),
    effect     = if (is_ratio) round(exp(p * max_eff), 4)
                 else          round(p * max_eff, 4)
  )
}


#' Pairwise league table at specified doses
#'
#' For each pair of doses (row dose vs column dose), computes the predicted
#' relative effect and 95% CI — analogous to a network meta-analysis league
#' table.
#'
#' Multiple `drma` objects can be passed to compare results across models
#' (e.g. primary analysis vs sensitivity analyses).
#'
#' @param ...    One or more `drma` objects (separated by commas).
#' @param doses  Numeric vector of doses forming the rows/columns of the table.
#' @param labels Optional character vector of labels for each model (used when
#'   multiple models are passed).
#' @param digits Number of decimal places (default 2).
#' @param sep    Separator between point estimate and CI in cells (default
#'   `" "` — newline-style labels are also fine).
#'
#' @return
#'   * **Single model**: a character matrix where cell `[i, j]` gives the
#'     effect of `doses[i]` relative to `doses[j]` (upper triangle) and its
#'     95% CI (lower triangle).  Diagonal cells show the dose value.
#'   * **Multiple models**: a list of such matrices, one per model.
#'
#' @examples
#' \dontrun{
#' lt <- league_table(res_primary, res_sensitivity,
#'                    doses  = c(0, 1, 2, 3),
#'                    labels = c("Primary", "Sensitivity"))
#' lt$Primary
#' }
#' @export
league_table <- function(..., doses, labels = NULL, digits = 2, sep = " ") {
  models <- list(...)

  # Validate that all inputs are drma objects
  is_drma <- vapply(models, inherits, logical(1), what = "drma")
  if (!all(is_drma))
    stop("All positional arguments must be 'drma' objects.")

  if (is.null(labels))
    labels <- paste0("Model", seq_along(models))

  make_table <- function(mod) {
    is_ratio <- mod$measure %in% c("OR", "RR")

    # Predict log-effects at each dose (relative to dose[1] as common base)
    base  <- doses[1]
    nd    <- data.frame(.dose = doses)
    pred  <- predict(mod$model, newdata = nd, xref = base, exp = FALSE)

    rcs_v <- .rcs_col(pred)
    log_e <- pred$pred    # log-effect vs base
    log_lo <- pred$ci.lb  # lower bound (already on log scale from dosresmeta)
    log_hi <- pred$ci.ub

    k   <- length(doses)
    mat <- matrix("", nrow = k, ncol = k,
                  dimnames = list(as.character(doses), as.character(doses)))

    for (i in seq_len(k)) {
      for (j in seq_len(k)) {
        if (i == j) {
          mat[i, j] <- as.character(doses[i])
          next
        }
        # Effect of dose i vs dose j  (= log_e[i] - log_e[j])
        le   <- log_e[i]  - log_e[j]
        # Delta-method-based CI: variance = var_i + var_j - 2*cov_ij
        # We approximate CI width from the individual half-widths (conservative)
        hw_i <- (log_hi[i] - log_lo[i]) / (2 * 1.96)
        hw_j <- (log_hi[j] - log_lo[j]) / (2 * 1.96)
        se_ij <- sqrt(hw_i^2 + hw_j^2)
        lo   <- le - 1.96 * se_ij
        hi   <- le + 1.96 * se_ij

        if (is_ratio) {
          fmt_val <- function(v) format(round(exp(v), digits), nsmall = digits)
        } else {
          fmt_val <- function(v) format(round(v, digits), nsmall = digits)
        }

        mat[i, j] <- paste0(fmt_val(le), sep,
                             "(", fmt_val(lo), ", ", fmt_val(hi), ")")
      }
    }
    mat
  }

  if (length(models) == 1) {
    return(make_table(models[[1]]))
  }

  result <- lapply(models, make_table)
  names(result) <- labels
  result
}
