#' Dose-Response Meta-Analysis
#'
#' Fits a dose-response meta-analysis model using the `dosresmeta` engine.
#' Accepts arm-level data in long format and automatically computes OR, RR,
#' MD, or SMD.  The interface is modelled on `meta`/`netmeta`.
#'
#' @param data     A `data.frame` in **long format** — one row per arm.
#' @param studlab  Column name for the study identifier.
#'   Accepts bare (unquoted) names or quoted strings.
#' @param dose     Column name for the dose variable.
#'   Accepts bare (unquoted) names or quoted strings.
#' @param sm       Summary measure: `"OR"`, `"RR"`, `"MD"`, `"SMD"`,
#'   or `"precomputed"`.
#' @param event    Column name for event count — required for
#'   `sm = "OR"` or `"RR"`.
#' @param n        Column name for arm sample size.
#' @param mean     Column name for arm mean — required for
#'   `sm = "MD"` or `"SMD"`.
#' @param sd       Column name for arm SD — required for
#'   `sm = "MD"` or `"SMD"`.
#' @param yi       Column name for pre-computed log-effect
#'   (`sm = "precomputed"`).
#' @param sei      Column name for pre-computed SE
#'   (`sm = "precomputed"`).
#' @param precomputed_scale Scale of pre-computed effects:
#'   `"ratio"` (default) for log-OR / log-RR — y-axis is exponentiated and
#'   plotted on a log scale; `"difference"` for MD / SMD — y-axis is linear.
#' @param ref      Reference dose.
#'   * `NULL` (default): minimum dose within each study.
#'   * Single number: same reference for all studies.
#'   * Named numeric vector (names = study IDs): per-study references.
#' @param curve    Functional form of the dose-response curve:
#'   * `"rcs"` (default): restricted cubic spline.
#'   * `"linear"`: straight line, `effect ~ dose`.
#'   * `"log"`: log-linear, `effect ~ log(dose + log_shift)`.
#'   * `"quadratic"`: `effect ~ dose + dose^2`.
#' @param knots    Knot positions for `curve = "rcs"`.
#'   * **String** `"p1-p2-p3"` — percentile probabilities, e.g.
#'     `"0.1-0.5-0.9"` or `"0.25-0.50-0.75"`.
#'   * **Numeric vector** — actual dose values, e.g. `c(1, 2, 3)`.
#'   * **Single integer** — e.g. `3L`: that many knots placed
#'     at evenly-spaced quantiles (10th–90th percentile).
#' @param log_shift Shift added before log-transformation when
#'   `curve = "log"`: `log(dose + log_shift)`.  Default `1`.
#' @param method   Estimation method for `dosresmeta` (default `"ml"`).
#' @param proc     Procedure: `"1stage"` (default) or `"2stage"`.
#' @param type     Overrides the `type` argument of `dosresmeta`.
#'   If `NULL` (default) it is inferred from `sm`.
#' @param zero_add Continuity correction for zero cells in binary
#'   outcomes (default `0.5`).
#' @param ...      Additional arguments passed to
#'   [dosresmeta::dosresmeta()].
#'
#' @return An object of class `"drma"` with components:
#'   \item{model}{The fitted `dosresmeta` object.}
#'   \item{data}{Processed data including computed effect sizes.}
#'   \item{data_fit}{Subset of `data` used for model fitting.}
#'   \item{sm}{Summary measure.}
#'   \item{curve}{Curve type.}
#'   \item{knots}{Resolved knot positions (for `curve = "rcs"`).}
#'   \item{log_shift}{Log shift (for `curve = "log"`).}
#'   \item{type}{`type` passed to `dosresmeta`.}
#'   \item{is_ratio}{`TRUE` when the effect measure is on a ratio / log scale
#'     (OR, RR, or precomputed with `precomputed_scale = "ratio"`).}
#'   \item{call}{The matched call.}
#'
#' @seealso [plot.drma()], [predict_table()], [target_dose()],
#'   [league_table()]
#'
#' @examples
#' \dontrun{
#' data(brexpiprazole)
#'
#' # ── Binary outcome (OR) — arm-level counts ────────────────────────
#' res <- drma(
#'   data    = brexpiprazole,
#'   studlab = study_id,
#'   dose    = dose,
#'   sm      = "OR",
#'   event   = n_responders,
#'   n       = n_arm,
#'   knots   = c(1, 2, 3)
#' )
#'
#' # Knots as percentile string "p1-p2-p3"
#' res <- drma(data = brexpiprazole, studlab = study_id, dose = dose,
#'             sm = "OR", event = n_responders, n = n_arm,
#'             knots = "0.1-0.5-0.9")
#'
#' # ── Precomputed log-effects ───────────────────────────────────────
#' res_t <- drma(
#'   data    = brexpiprazole,
#'   studlab = study_id,
#'   dose    = dose,
#'   sm      = "precomputed",
#'   yi      = tolerability_logor,
#'   sei     = tolerability_se,
#'   event   = n_dropout_ae,
#'   n       = n_arm,
#'   knots   = c(1, 2, 3)
#' )
#'
#' # ── Standard workflow ─────────────────────────────────────────────
#' print(res)
#' plot(res, ylab = "Response (OR)", ylim = c(0.75, 2), ref_dose = 0,
#'      bubble = TRUE, rug = TRUE)
#' target_dose(res, p = c(0.5, 0.95, 1))
#' predict_table(res, doses = c(0, 1, 2, 3), baseline_prop = 0.183)
#' league_table(res, doses = c(0, 1, 2, 3))
#'
#' # ── Overlay sensitivity analyses ──────────────────────────────────
#' res_s1 <- drma(data = brexpiprazole, studlab = study_id, dose = dose,
#'                sm = "OR", event = n_responders, n = n_arm,
#'                knots = "0.1-0.5-0.9")
#' res_s2 <- drma(data = brexpiprazole, studlab = study_id, dose = dose,
#'                sm = "OR", event = n_responders, n = n_arm,
#'                knots = "0.25-0.50-0.75")
#' plot(res,  col = "black", ylim = c(0.75, 2), ref_dose = 0)
#' lines(res_s1, col = "tomato",    lty = 2)
#' lines(res_s2, col = "steelblue", lty = 3)
#' legend("topright",
#'   legend = c("Primary c(1,2,3)", "S1 10/50/90%", "S2 25/50/75%"),
#'   col = c("black", "tomato", "steelblue"), lty = 1:3, bty = "n")
#' }
#' @export
drma <- function(
  data,
  studlab   = NULL,
  dose      = NULL,
  sm        = c("OR", "RR", "MD", "SMD", "precomputed"),
  event     = NULL,
  n         = NULL,
  mean      = NULL,
  sd        = NULL,
  yi        = NULL,
  sei       = NULL,
  ref       = NULL,
  curve     = c("rcs", "linear", "log", "quadratic"),
  knots     = "0.25-0.50-0.75",
  log_shift = 1,
  method    = "ml",
  proc      = "1stage",
  type      = NULL,
  zero_add          = 0.5,   # continuity correction for zero cells; change only if needed
  precomputed_scale = c("ratio", "difference"),
  ...
) {
  sm                <- match.arg(sm)
  curve             <- match.arg(curve)
  precomputed_scale <- match.arg(precomputed_scale)
  cl                <- match.call()

  # ── 1. Resolve column names (quoted "x" or bare x both accepted) ───
  studlab <- .as_col(substitute(studlab))
  dose    <- .as_col(substitute(dose))
  event   <- .as_col(substitute(event))
  n       <- .as_col(substitute(n))
  mean    <- .as_col(substitute(mean))
  sd      <- .as_col(substitute(sd))
  yi      <- .as_col(substitute(yi))
  sei     <- .as_col(substitute(sei))

  d <- as.data.frame(data)
  if (is.null(studlab))
    stop("'studlab' must be specified (e.g. studlab = study_id).")
  if (!studlab %in% names(d))
    stop("Column '", studlab, "' not found in data.")
  if (is.null(dose))
    stop("'dose' must be specified (e.g. dose = dose).")
  if (!dose %in% names(d))
    stop("Column '", dose, "' not found in data.")
  d$.id   <- d[[studlab]]
  d$.dose <- as.numeric(d[[dose]])

  # ── 2. Compute effect sizes ────────────────────────────────────────
  if (sm == "precomputed") {
    if (is.null(yi) || is.null(sei))
      stop("'yi' and 'sei' must be provided when sm = 'precomputed'.")
    d$.yi    <- as.numeric(d[[yi]])
    d$.sei   <- as.numeric(d[[sei]])
    # Reference arms in published tables often have yi = 0 and sei = NA.
    # dosresmeta treats variance = 0 as a reference arm, so keep these rows
    # and set sei = 0 so they are not dropped by the NA filter below.
    is_ref_arm <- !is.na(d$.yi) & d$.yi == 0 & is.na(d$.sei)
    d$.sei[is_ref_arm] <- 0
    d$.cases <- if (!is.null(event)) as.numeric(d[[event]]) else NULL
    d$.n     <- if (!is.null(n))     as.numeric(d[[n]])     else NULL
    auto_type <- "cc"

  } else if (sm %in% c("OR", "RR")) {
    if (is.null(event) || is.null(n))
      stop("'event' and 'n' required for sm = '", sm, "'.")
    d$.cases  <- as.numeric(d[[event]])
    d$.n      <- as.numeric(d[[n]])
    auto_type <- if (sm == "OR") "cc" else "ci"
    d         <- .compute_binary(d, ref, sm, zero_add)

  } else {  # MD / SMD
    if (is.null(mean) || is.null(sd) || is.null(n))
      stop("'mean', 'sd', and 'n' required for sm = '", sm, "'.")
    d$.mean   <- as.numeric(d[[mean]])
    d$.sd     <- as.numeric(d[[sd]])
    d$.n      <- as.numeric(d[[n]])
    auto_type <- "d"
    d         <- .compute_continuous(d, ref, sm)
  }

  type_use <- if (!is.null(type)) type else auto_type

  # ── 3. Resolve knots (rcs only) ────────────────────────────────────
  dose_vals      <- d$.dose[!is.na(d$.yi)]
  resolved_knots <- if (curve == "rcs") {
    .resolve_knots(knots, dose_vals)
  } else {
    NULL
  }

  # ── 4. Build formula ───────────────────────────────────────────────
  # The formula environment (= current function env) keeps
  # resolved_knots and log_shift accessible during model.frame().
  fm <- switch(curve,
    rcs       = as.formula(".yi ~ rcs(.dose, resolved_knots)"),
    linear    = as.formula(".yi ~ .dose"),
    log       = as.formula(
      paste0(".yi ~ log(.dose + ", log_shift, ")")
    ),
    quadratic = as.formula(".yi ~ .dose + I(.dose^2)")
  )

  # ── 5. Subset to valid rows ────────────────────────────────────────
  # Keep rows with non-NA yi and sei, then drop studies that have only a
  # reference arm (yi = 0, sei = eps) with no non-reference arm — such
  # studies occur when all event counts are NA and would make dosresmeta's
  # within-study covariance computation fail.
  d_fit <- d[!is.na(d$.yi) & !is.na(d$.sei), ]
  # Reference arms have sei = 0; non-reference arms have sei > 0.
  # Drop studies where only the reference arm survived NA removal.
  arms_per_study <- tapply(d_fit$.sei > 0, d_fit$.id, sum)
  valid_ids <- names(arms_per_study)[arms_per_study > 0]
  excluded  <- setdiff(unique(d_fit$.id), valid_ids)
  if (length(excluded) > 0)
    message(
      "Excluded ", length(excluded), " study/studies with no non-reference ",
      "arm after NA removal: ",
      paste(excluded, collapse = ", ")
    )
  d_fit <- d_fit[d_fit$.id %in% valid_ids, ]

  # ── 6. Call dosresmeta ─────────────────────────────────────────────
  args <- list(
    formula = fm,
    type    = type_use,
    id      = d_fit$.id,
    se      = d_fit$.sei,
    data    = d_fit,
    method  = method,
    proc    = proc,
    ...
  )
  if (!is.null(d_fit$.cases)) args$cases <- d_fit$.cases
  if (!is.null(d_fit$.n))     args$n     <- d_fit$.n

  mod <- tryCatch(
    do.call(dosresmeta::dosresmeta, args),
    error = function(e) stop("dosresmeta error: ", conditionMessage(e))
  )

  # ── 7. Return ──────────────────────────────────────────────────────
  is_ratio_val <- sm %in% c("OR", "RR") ||
    (sm == "precomputed" && precomputed_scale == "ratio")

  structure(
    list(
      model     = mod,
      data      = d,
      data_fit  = d_fit,
      sm        = sm,
      is_ratio  = is_ratio_val,
      curve     = curve,
      knots     = resolved_knots,
      log_shift = log_shift,
      type      = type_use,
      call      = cl
    ),
    class = "drma"
  )
}


# ── S3: print ─────────────────────────────────────────────────────────────────

#' @export
print.drma <- function(x, ...) {
  n_studies <- length(unique(x$data$.id))
  n_arms    <- nrow(x$data_fit)

  # Header
  knot_str <- if (x$curve == "rcs" && !is.null(x$knots))
    paste0("  knots: ", paste(round(x$knots, 3), collapse = ", "))
  else if (x$curve == "log")
    paste0("  log_shift: ", x$log_shift)
  else ""

  cat(sprintf(
    "drma  sm=%s  curve=%s%s  method=%s  studies=%d  arms=%d\n",
    x$sm, x$curve, knot_str, x$model$method, n_studies, n_arms
  ))

  # Coefficient table with 95% CI
  b  <- stats::coef(x$model)
  se <- sqrt(diag(stats::vcov(x$model)))
  lo <- b - 1.96 * se
  hi <- b + 1.96 * se
  tbl <- data.frame(
    est   = round(b,  3),
    lower = round(lo, 3),
    upper = round(hi, 3),
    row.names = names(b)
  )
  names(tbl) <- c("est", "2.5%", "97.5%")
  print(tbl)

  # Fit stats
  ll  <- round(as.numeric(stats::logLik(x$model)), 1)
  aic <- round(stats::AIC(x$model), 1)
  cat(sprintf("logLik: %s  AIC: %s\n", ll, aic))

  invisible(x)
}


# ── S3: summary ───────────────────────────────────────────────────────────────

#' @export
summary.drma <- function(object, ...) {
  cat("=== Dose-Response Meta-Analysis ===\n")
  cat(sprintf("Measure : %s\n", object$sm))
  cat(sprintf("Curve   : %s\n", object$curve))
  if (object$curve == "rcs" && !is.null(object$knots))
    cat(sprintf(
      "Knots   : %s\n",
      paste(round(object$knots, 3), collapse = ", ")
    ))
  cat(sprintf("Studies : %d\n\n", length(unique(object$data$.id))))
  summary(object$model, ...)
}
