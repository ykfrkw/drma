#' Dose-Response Meta-Analysis
#'
#' Fits a dose-response meta-analysis model using restricted cubic splines.
#' Accepts arm-level data in long format and automatically computes effect
#' sizes (OR, RR, MD, SMD) before passing to [dosresmeta::dosresmeta()].
#'
#' @param data      A `data.frame` in long format — one row per arm per study.
#' @param study     Column name for the study identifier.
#' @param dose      Column name for the dose variable.
#' @param measure   Effect measure: `"OR"`, `"RR"`, `"MD"`, `"SMD"`,
#'   or `"precomputed"`.
#' @param events    Column name for event count (required for OR / RR).
#' @param n         Column name for arm sample size.
#' @param mean      Column name for arm mean (required for MD / SMD).
#' @param sd        Column name for arm SD (required for MD / SMD).
#' @param yi        Column name for pre-computed log-effect (`"precomputed"`).
#' @param sei       Column name for pre-computed SE (`"precomputed"`).
#' @param ref       Reference dose.  `NULL` (default) uses the minimum dose in
#'   each study.  A single number applies to all studies.  A named numeric
#'   vector (names = study IDs) sets per-study references.
#' @param knots     Knot positions for the restricted cubic spline.
#'   * Numeric vector with all values in (0, 1): treated as quantile
#'     probabilities of the observed doses (e.g. `c(0.25, 0.5, 0.75)`).
#'   * Numeric vector with any value >= 1: treated as actual dose values
#'     (e.g. `c(1, 2, 3)`).
#'   * Single integer >= 2: automatically places that many knots at
#'     evenly-spaced quantiles (e.g. `3L`).
#' @param method    Estimation method passed to `dosresmeta` (default `"ml"`).
#' @param proc      Procedure passed to `dosresmeta`: `"1stage"` (default) or
#'   `"2stage"`.
#' @param type      Override the `type` argument of `dosresmeta`.  If `NULL`
#'   (default) it is inferred from `measure` (`"cc"` for OR, `"ci"` for RR,
#'   `"d"` for MD/SMD).
#' @param zero_add  Continuity correction added to zero cells in binary
#'   outcomes (default `0.5`).
#' @param ...       Additional arguments passed to [dosresmeta::dosresmeta()].
#'
#' @return An object of class `"drma"` with components:
#'   \item{model}{The fitted `dosresmeta` object.}
#'   \item{data}{Processed long-format data including computed effect sizes.}
#'   \item{data_fit}{Subset of `data` used for model fitting.}
#'   \item{measure}{Effect measure used.}
#'   \item{knots}{Resolved knot positions.}
#'   \item{type}{`type` argument passed to `dosresmeta`.}
#'   \item{call}{The matched call.}
#'
#' @examples
#' \dontrun{
#' res <- drma(
#'   data    = mydata,
#'   study   = "studyID",
#'   dose    = "dose",
#'   measure = "OR",
#'   events  = "n_events",
#'   n       = "n_total",
#'   knots   = c(1, 2, 3)
#' )
#' plot(res)
#' }
#' @export
drma <- function(
    data,
    study    = "studyID",
    dose     = "dose",
    measure  = c("OR", "RR", "MD", "SMD", "precomputed"),
    events   = NULL,
    n        = NULL,
    mean     = NULL,
    sd       = NULL,
    yi       = NULL,
    sei      = NULL,
    ref      = NULL,
    knots    = c(0.25, 0.5, 0.75),
    method   = "ml",
    proc     = "1stage",
    type     = NULL,
    zero_add = 0.5,
    ...
) {
  measure <- match.arg(measure)
  cl      <- match.call()

  # ── 1. Standardise column names ──────────────────────────────────────────
  d <- as.data.frame(data)
  if (!study %in% names(d)) stop("Column '", study, "' not found in data.")
  if (!dose  %in% names(d)) stop("Column '", dose,  "' not found in data.")
  d$.id   <- d[[study]]
  d$.dose <- as.numeric(d[[dose]])

  # ── 2. Compute effect sizes ───────────────────────────────────────────────
  if (measure == "precomputed") {
    if (is.null(yi) || is.null(sei))
      stop("'yi' and 'sei' must be provided when measure = 'precomputed'.")
    d$.yi    <- as.numeric(d[[yi]])
    d$.sei   <- as.numeric(d[[sei]])
    d$.cases <- if (!is.null(events)) as.numeric(d[[events]]) else NULL
    d$.n     <- if (!is.null(n))      as.numeric(d[[n]])      else NULL
    auto_type <- "cc"

  } else if (measure %in% c("OR", "RR")) {
    if (is.null(events) || is.null(n))
      stop("'events' and 'n' required for measure = '", measure, "'.")
    d$.cases <- as.numeric(d[[events]])
    d$.n     <- as.numeric(d[[n]])
    auto_type <- if (measure == "OR") "cc" else "ci"
    d         <- .compute_binary(d, ref, measure, zero_add)

  } else {  # MD / SMD
    if (is.null(mean) || is.null(sd) || is.null(n))
      stop("'mean', 'sd', and 'n' required for measure = '", measure, "'.")
    d$.mean <- as.numeric(d[[mean]])
    d$.sd   <- as.numeric(d[[sd]])
    d$.n    <- as.numeric(d[[n]])
    auto_type <- "d"
    d         <- .compute_continuous(d, ref, measure)
  }

  type_use <- if (!is.null(type)) type else auto_type

  # ── 3. Resolve knots ──────────────────────────────────────────────────────
  dose_vals      <- d$.dose[!is.na(d$.yi)]
  resolved_knots <- .resolve_knots(knots, dose_vals)

  # ── 4. Build formula (environment captures resolved_knots) ───────────────
  fm <- as.formula(".yi ~ rcs(.dose, resolved_knots)")

  # ── 5. Subset to valid rows ───────────────────────────────────────────────
  d_fit <- d[!is.na(d$.yi) & !is.na(d$.sei), ]

  # ── 6. Call dosresmeta ────────────────────────────────────────────────────
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

  # ── 7. Return ─────────────────────────────────────────────────────────────
  structure(
    list(
      model    = mod,
      data     = d,
      data_fit = d_fit,
      measure  = measure,
      knots    = resolved_knots,
      type     = type_use,
      call     = cl
    ),
    class = "drma"
  )
}


# ── S3: print ────────────────────────────────────────────────────────────────

#' @export
print.drma <- function(x, ...) {
  cat("Dose-Response Meta-Analysis\n")
  cat(strrep("-", 38), "\n")
  cat(sprintf("  Measure : %s\n",   x$measure))
  cat(sprintf("  Knots   : %s\n",   paste(round(x$knots, 3), collapse = ", ")))
  cat(sprintf("  Method  : %s\n",   x$model$method))
  cat(sprintf("  Studies : %d\n\n", length(unique(x$data$.id))))
  print(x$model)
  invisible(x)
}


# ── S3: summary ──────────────────────────────────────────────────────────────

#' @export
summary.drma <- function(object, ...) {
  cat("=== Dose-Response Meta-Analysis ===\n")
  cat(sprintf("Measure : %s\n",   object$measure))
  cat(sprintf("Knots   : %s\n",   paste(round(object$knots, 3), collapse = ", ")))
  cat(sprintf("Studies : %d\n\n", length(unique(object$data$.id))))
  summary(object$model, ...)
}
