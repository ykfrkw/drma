#' Dose-Response Meta-Analysis
#'
#' Fits a dose-response meta-analysis model using the `dosresmeta` engine.
#' Accepts arm-level data in long format and automatically computes OR, RR,
#' MD, or SMD.  The interface is modelled on `meta`/`netmeta`.
#'
#' @param data     A `data.frame` in **long format** — one row per arm.
#' @param studlab  Column name for the study identifier
#'   (default `"studyID"`).
#' @param dose     Column name for the dose variable (default `"dose"`).
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
#'   \item{call}{The matched call.}
#'
#' @seealso [plot.drma()], [predict_table()], [target_dose()],
#'   [league_table()]
#'
#' @examples
#' \dontrun{
#' # ── Binary outcome (OR) ───────────────────────────────────────────
#' res <- drma(
#'   data    = mydata,
#'   studlab = "studyID",
#'   dose    = "dose",
#'   sm      = "OR",
#'   event   = "n_events",
#'   n       = "n_total",
#'   knots   = c(1, 2, 3)           # actual dose values
#' )
#'
#' # Knots as percentiles — string "p1-p2-p3"
#' res <- drma(..., knots = "0.1-0.5-0.9")
#' res <- drma(..., knots = "0.25-0.50-0.75")
#'
#' # Auto-place 3 knots (integer)
#' res <- drma(..., knots = 3L)
#'
#' # ── Continuous outcome (MD) ───────────────────────────────────────
#' res_md <- drma(
#'   data    = mydata,
#'   studlab = "studyID",
#'   dose    = "dose",
#'   sm      = "MD",
#'   mean    = "mean_arm",
#'   sd      = "sd_arm",
#'   n       = "n_arm",
#'   knots   = "0.25-0.50-0.75"
#' )
#'
#' # ── Log-linear curve ──────────────────────────────────────────────
#' res_log <- drma(..., curve = "log", log_shift = 1)
#'
#' # ── Pre-computed log-effects ──────────────────────────────────────
#' res_pre <- drma(
#'   data    = df,
#'   studlab = "id",
#'   dose    = "dose",
#'   sm      = "precomputed",
#'   yi      = "logor",
#'   sei     = "se",
#'   event   = "cases",
#'   n       = "n",
#'   knots   = c(1.5, 3, 6)
#' )
#'
#' # ── Standard workflow ─────────────────────────────────────────────
#' print(res)
#' plot(res, ylab = "Response (OR)", ylim = c(0.5, 3))
#' target_dose(res, p = c(0.5, 0.95, 1))
#' predict_table(res, doses = c(0, 1, 2, 3), baseline_prop = 0.30)
#' league_table(res, doses = c(0, 1, 2, 3))
#'
#' # ── Overlay sensitivity analyses ──────────────────────────────────
#' res_s1 <- drma(..., knots = "0.1-0.5-0.9")
#' res_s2 <- drma(..., knots = "0.25-0.50-0.75")
#' plot(res,  col = "black", ylim = c(0.5, 3))
#' lines(res_s1, col = "red",  lty = 2)
#' lines(res_s2, col = "blue", lty = 3)
#' legend("topright",
#'   legend = c("Primary c(1,2,3)", "S1 10/50/90%", "S2 25/50/75%"),
#'   col = c("black", "red", "blue"), lty = 1:3)
#' }
#' @export
drma <- function(
  data,
  studlab   = "studyID",
  dose      = "dose",
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
  zero_add  = 0.5,   # continuity correction for zero cells; change only if needed
  ...
) {
  sm    <- match.arg(sm)
  curve <- match.arg(curve)
  cl    <- match.call()

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
  if (!studlab %in% names(d))
    stop("Column '", studlab, "' not found in data.")
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
  d_fit <- d[!is.na(d$.yi) & !is.na(d$.sei), ]

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
  structure(
    list(
      model     = mod,
      data      = d,
      data_fit  = d_fit,
      sm        = sm,
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
  cat("Dose-Response Meta-Analysis\n")
  cat(strrep("-", 38), "\n")
  cat(sprintf("  Measure : %s\n", x$sm))
  cat(sprintf("  Curve   : %s\n", x$curve))
  if (x$curve == "rcs" && !is.null(x$knots))
    cat(sprintf(
      "  Knots   : %s\n",
      paste(round(x$knots, 3), collapse = ", ")
    ))
  if (x$curve == "log")
    cat(sprintf("  LogShift: %s\n", x$log_shift))
  cat(sprintf("  Method  : %s\n",   x$model$method))
  cat(sprintf("  Studies : %d\n\n", length(unique(x$data$.id))))
  print(x$model)
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
