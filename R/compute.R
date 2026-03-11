# Internal functions for computing effect sizes from arm-level data

# ── Binary outcomes (OR / RR) ─────────────────────────────────────────────────
.compute_binary <- function(d, ref, measure, zero_add) {
  d$.events_ref <- NA_real_
  d$.n_ref      <- NA_real_
  d$.ref_dose   <- NA_real_
  d$.yi         <- NA_real_
  d$.sei        <- NA_real_

  for (s in unique(d$.id)) {
    idx   <- which(d$.id == s)
    ds    <- d[idx, ]
    rdose <- .get_ref_dose(ds$.dose, ref, s)
    d$.ref_dose[idx] <- rdose

    ri <- which(ds$.dose == rdose)[1]
    if (length(ri) == 0)
      stop("Reference dose not found in study '", s, "'.")

    # Continuity correction for zero cells
    if (any(ds$.cases == 0, na.rm = TRUE)) {
      ds$.cases <- ds$.cases + zero_add
      ds$.n     <- ds$.n     + zero_add
      d$.cases[idx] <- ds$.cases
      d$.n[idx]     <- ds$.n
      message("Study '", s, "': zero cell — added ", zero_add, " to events and n.")
    }

    ev_ref <- ds$.cases[ri]
    n_ref  <- ds$.n[ri]
    d$.events_ref[idx] <- ev_ref
    d$.n_ref[idx]      <- n_ref

    if (measure == "OR") {
      d$.yi[idx]  <- log(ds$.cases / (ds$.n - ds$.cases)) -
                     log(ev_ref   / (n_ref  - ev_ref))
      d$.sei[idx] <- sqrt(1 / ds$.cases + 1 / (ds$.n - ds$.cases) +
                          1 / ev_ref    + 1 / (n_ref  - ev_ref))
    } else {  # RR
      d$.yi[idx]  <- log(ds$.cases / ds$.n) - log(ev_ref / n_ref)
      d$.sei[idx] <- sqrt(1 / ds$.cases - 1 / ds$.n +
                          1 / ev_ref    - 1 / n_ref)
    }

    # Reference arm: contrast = 0
    d$.yi [idx[ri]] <- 0
    d$.sei[idx[ri]] <- .Machine$double.eps
  }
  d
}


# ── Continuous outcomes (MD / SMD) ────────────────────────────────────────────
.compute_continuous <- function(d, ref, measure) {
  d$.mean_ref <- NA_real_
  d$.sd_ref   <- NA_real_
  d$.n_ref    <- NA_real_
  d$.ref_dose <- NA_real_
  d$.yi       <- NA_real_
  d$.sei      <- NA_real_

  for (s in unique(d$.id)) {
    idx   <- which(d$.id == s)
    ds    <- d[idx, ]
    rdose <- .get_ref_dose(ds$.dose, ref, s)
    d$.ref_dose[idx] <- rdose

    ri    <- which(ds$.dose == rdose)[1]
    if (length(ri) == 0)
      stop("Reference dose not found in study '", s, "'.")

    m_ref <- ds$.mean[ri]
    s_ref <- ds$.sd[ri]
    n_ref <- ds$.n[ri]
    d$.mean_ref[idx] <- m_ref
    d$.sd_ref[idx]   <- s_ref
    d$.n_ref[idx]    <- n_ref

    if (measure == "MD") {
      d$.yi[idx]  <- ds$.mean - m_ref
      d$.sei[idx] <- sqrt(ds$.sd^2 / ds$.n + s_ref^2 / n_ref)

    } else {  # SMD (Hedges' g)
      df_g <- ds$.n + n_ref - 2
      sp   <- sqrt(((ds$.n - 1) * ds$.sd^2 + (n_ref - 1) * s_ref^2) / df_g)
      g    <- (ds$.mean - m_ref) / sp
      J    <- 1 - 3 / (4 * df_g - 1)   # Hedges correction factor
      d$.yi[idx]  <- J * g
      d$.sei[idx] <- sqrt(ds$.n / n_ref / (ds$.n + n_ref) +
                          g^2 / (2 * df_g))
    }

    d$.yi [idx[ri]] <- 0
    d$.sei[idx[ri]] <- .Machine$double.eps
  }
  d
}


# ── Helpers ───────────────────────────────────────────────────────────────────

# Get reference dose for a study
.get_ref_dose <- function(doses, ref, study_id) {
  if (is.null(ref))                             return(min(doses, na.rm = TRUE))
  if (is.numeric(ref) && length(ref) == 1)      return(ref)
  if (!is.null(names(ref))) {
    key <- as.character(study_id)
    if (key %in% names(ref)) return(unname(ref[[key]]))
  }
  min(doses, na.rm = TRUE)  # fallback
}


# Resolve knot specification to actual dose values
#
# knots_type:
#   "auto"     – detect from values:
#                  single integer >= 2  -> equally spaced quantiles
#                  all in (0, 1)        -> quantile probabilities
#                  otherwise            -> actual dose values
#   "quantile" – always treat as quantile probabilities
#   "values"   – always treat as actual dose values
.resolve_knots <- function(knots, dose_vals,
                           knots_type = c("auto", "quantile", "values")) {
  knots_type <- match.arg(knots_type)
  dose_vals  <- dose_vals[!is.na(dose_vals)]

  if (knots_type == "quantile") {
    if (any(knots <= 0 | knots >= 1))
      stop("'knots' must be in (0, 1) when knots_type = 'quantile'.")
    return(unname(stats::quantile(dose_vals, probs = knots)))
  }

  if (knots_type == "values") {
    return(as.numeric(knots))
  }

  # ── auto detection ──────────────────────────────────────────────────────
  if (length(knots) == 1L && knots == round(knots) && knots >= 2L) {
    # Single integer -> equally spaced quantiles (10th to 90th percentile)
    probs <- seq(0.1, 0.9, length.out = as.integer(knots))
    return(unname(stats::quantile(dose_vals, probs = probs)))
  }

  if (all(knots > 0 & knots < 1)) {
    # All in (0, 1) -> interpret as quantile probabilities
    return(unname(stats::quantile(dose_vals, probs = knots)))
  }

  as.numeric(knots)  # actual dose values
}


# Extract the rcs() column name from a predict() result
.rcs_col <- function(pred) {
  grep("^rcs", names(pred), value = TRUE)[1]
}


# Human-readable effect label
.effect_label <- function(measure) {
  switch(measure,
    OR  = "Odds Ratio",
    RR  = "Risk Ratio",
    MD  = "Mean Difference",
    SMD = "Standardised Mean Difference",
    "Effect"
  )
}
