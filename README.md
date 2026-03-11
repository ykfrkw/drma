# drma <img src="https://img.shields.io/badge/version-0.1.0-blue" align="right"/>

**Dose-Response Meta-Analysis** — a `meta`/`netmeta`-style wrapper around
[`dosresmeta`](https://cran.r-project.org/package=dosresmeta) for
**one-stage dose-response meta-analysis**.

Pass arm-level data in long format and `drma` automatically computes OR,
RR, MD, or SMD before fitting a dose-response curve with restricted cubic
splines (or linear / log-linear / quadratic alternatives).

> **Statistical engine**: all modelling is performed by `dosresmeta`
> (Crippa & Orsini 2016).  `drma` provides the data-preparation,
> knot-resolution, plot, and utility layer on top of it.

> **Version 0.1.0** — first public release.  API may change in future versions.

## Installation

```r
# requires remotes
remotes::install_github("ykfrkw/drma")
```

---

## Complete example (copy-paste and run top to bottom)

Uses the `brexpiprazole` dataset bundled with the package
(Furukawa et al. 2022, *Psychiatry Clin Neurosci* 76:416–422).

```r
library(drma)
library(ggplot2)

# ── 0. Data ───────────────────────────────────────────────────────────────────
data(brexpiprazole)
brexpiprazole[, c("study_id", "dose", "n_arm", "n_responders")]

# ── 1. Fit (efficacy: arm-level counts → log-OR computed automatically) ───────
res <- drma(
  data    = brexpiprazole,
  studlab = study_id,
  dose    = dose,
  sm      = "OR",
  event   = n_responders,
  n       = n_arm,
  knots   = c(1, 2, 3)        # RCS knots at 1 / 2 / 3 mg
)

# ── 2. Inspect ────────────────────────────────────────────────────────────────
print(res)     # coefficients + 95% CI + logLik/AIC
summary(res)   # full dosresmeta summary

# ── 3. Plot ───────────────────────────────────────────────────────────────────
# ylim is always in display scale (OR/RR units, not log)
plot(res,
     ylab     = "Response (OR)",
     ylim     = c(0.75, 2.5),
     xlab     = "Brexpiprazole (mg)",
     ref_dose = 0,
     bubble   = TRUE,   # circles sized by arm N
     rug      = TRUE)   # tick marks for observed doses

# ── 4. Non-placebo reference ──────────────────────────────────────────────────
# ref_dose accepts any dose value — bubbles are re-centred automatically
plot(res,
     ylab     = "Response vs 1 mg (OR)",
     ylim     = c(0.5, 2),
     xlab     = "Brexpiprazole (mg)",
     ref_dose = 1,
     bubble   = TRUE)

# ── 5. Predictions at clinically relevant doses ───────────────────────────────
predict_table(
  res,
  doses         = c(0, 0.5, 1, 1.5, 2, 3),
  ref_dose      = 0,
  baseline_prop = 0.183   # converts OR → absolute response rate
)

# ── 6. Target doses (ED50 / ED95 / ED100) ────────────────────────────────────
target_dose(res, p = c(0.5, 0.95, 1))

# ── 7. Pairwise league table ──────────────────────────────────────────────────
league_table(res, doses = c(0, 1, 2, 3))

# ── 8. Sensitivity analyses — overlay multiple curves ────────────────────────
# plot() returns a ggplot object; add curves with + lines()
res_s1 <- drma(data = brexpiprazole, studlab = study_id, dose = dose,
               sm = "OR", event = n_responders, n = n_arm,
               knots = c(0.5, 1.5, 2.5))
res_s2 <- drma(data = brexpiprazole, studlab = study_id, dose = dose,
               sm = "OR", event = n_responders, n = n_arm,
               knots = "0.25-0.50-0.75")
res_s3 <- drma(data = brexpiprazole, studlab = study_id, dose = dose,
               sm = "OR", event = n_responders, n = n_arm,
               curve = "log")

plot(res,   col = "black", lwd = 1.5,
     ylim = c(0.75, 2.5), ylab = "Response (OR)",
     xlab = "Brexpiprazole (mg)", ref_dose = 0) +
  lines(res_s1, col = "tomato",      lty = 2, lwd = 1) +
  lines(res_s2, col = "steelblue",   lty = 3, lwd = 1) +
  lines(res_s3, col = "forestgreen", lty = 4, lwd = 1) +
  annotate("text", x = 3, y = c(1.55, 1.32, 1.22, 1.44),
           hjust = 1, size = 3,
           label = c("Primary c(1,2,3)", "S1 c(0.5,1.5,2.5)",
                     "S2 25/50/75th pct", "S3 log-linear"),
           color = c("black", "tomato", "steelblue", "forestgreen"))

# ── 9. Precomputed log-OR (tolerability / acceptability) ─────────────────────
#    Use sm = "precomputed" when yi and sei are already available.
res_t <- drma(
  data    = brexpiprazole,
  studlab = study_id,
  dose    = dose,
  sm      = "precomputed",
  yi      = tolerability_logor,
  sei     = tolerability_se,
  event   = n_dropout_ae,
  n       = n_arm,
  knots   = c(1, 2, 3)
)

plot(res_t,
     ylab     = "Dropout for AEs (OR)",
     ylim     = c(0.2, 10),
     xlab     = "Brexpiprazole (mg)",
     ref_dose = 0,
     bubble   = TRUE,
     rug      = TRUE)

# ── 10. VPC plot (Variance Partitioning Coefficient) ─────────────────────────
vpc_plot(res)
```

---

## Argument overview

### `drma()`

| Argument | Role |
|---|---|
| `studlab` | Column name for study identifier |
| `dose` | Column name for dose |
| `sm` | Summary measure: `"OR"` / `"RR"` / `"MD"` / `"SMD"` / `"precomputed"` |
| `event` / `n` | Event count and sample size (OR / RR) |
| `mean` / `sd` / `n` | Arm mean, SD, and N (MD / SMD) |
| `yi` / `sei` | Pre-computed log-effect and SE (`sm = "precomputed"`) |
| `ref` | Reference dose (default: minimum per study) |
| `curve` | Curve shape (see below) |
| `knots` | RCS knot specification (see below) |

### `plot.drma()`

| Argument | Role |
|---|---|
| `ref_dose` | Reference dose for the y-axis (default: minimum observed dose; **any value accepted**) |
| `ylim` | y-axis limits in **display scale** — OR/RR in ratio units (e.g. `c(0.75, 2)`), MD/SMD in difference units |
| `xlim` | x-axis limits (default: `c(0, max(dose))`) |
| `bubble` | Overlay observed effects as bubbles sized by arm N (`FALSE`) |
| `rug` | Draw rug tick marks for observed doses (`FALSE`) |
| `ci` | Draw 95% confidence band (`TRUE`) |
| `col` / `lty` / `lwd` | Curve colour, line type, line width |

---

## Curve types (`curve =`)

| Value | Formula | Notes |
|---|---|---|
| `"rcs"` (default) | restricted cubic spline | knots required |
| `"linear"` | `effect ~ dose` | |
| `"log"` | `effect ~ log(dose + log_shift)` | `log_shift = 1` by default |
| `"quadratic"` | `effect ~ dose + dose²` | |

---

## Knot specification for RCS

| `knots =` | Interpretation |
|---|---|
| `"0.1-0.5-0.9"` | 10th / 50th / 90th percentile of observed doses |
| `"0.25-0.50-0.75"` | 25th / 50th / 75th percentile (default) |
| `c(1, 2, 3)` | Actual dose values |
| `3L` | 3 knots auto-placed at evenly-spaced quantiles |

---

## Notes

**`ylim` scale**: always specify in display units.  For OR/RR use the
exponentiated scale (e.g. `ylim = c(0.75, 2)`, **not** `c(log(0.75), log(2))`).

**`ref_dose`**: any value within the dose range is accepted, not just zero.
Bubbles are automatically re-centred to align with the shifted curve.

**Sensitivity overlay**: `plot()` returns a `ggplot` object.  Add further
curves with `+ lines()`:
```r
plot(res_primary, col = "black") +
  lines(res_s1, col = "tomato",    lty = 2) +
  lines(res_s2, col = "steelblue", lty = 3)
```

**Zero-cell correction**: when any arm has zero events, `0.5` is added
automatically to events and n (Haldane-Anscombe correction).  Change with
`zero_add = <value>` only if needed.

---

## Function reference

| Function | Description |
|---|---|
| `drma()` | Fit dose-response model |
| `print.drma()` | Coefficients and model info |
| `summary.drma()` | Full model summary |
| `plot.drma()` | Dose-response curve (CI, bubble, rug) — returns ggplot |
| `lines.drma()` | Add curve layer to existing plot (use with `+`) |
| `predict.drma()` | Raw predictions |
| `predict_table()` | Predictions at specific doses |
| `target_dose()` | ED50 / ED95 / ED100 |
| `league_table()` | Pairwise dose comparison matrix |
| `vpc_plot()` | Variance Partitioning Coefficient plot |

---

## Citation

When using `drma`, please cite the underlying statistical method and package:

**One-stage dose-response method:**

> Crippa A, Discacciati A, Bottai M, Spiegelman D, Orsini N.
> One-stage dose-response meta-analysis for aggregated data.
> *Stat Methods Med Res*. 2019;28(5):1579–1596.
> doi:[10.1177/0962280218773122](https://doi.org/10.1177/0962280218773122)

**`dosresmeta` package:**

> Crippa A, Orsini N.
> Dose-response meta-analysis of differences in means.
> *BMC Med Res Methodol*. 2016;16:91.
> doi:[10.1186/s12874-016-0189-0](https://doi.org/10.1186/s12874-016-0189-0)

---

## Dependencies

[`dosresmeta`](https://cran.r-project.org/package=dosresmeta) ·
[`rms`](https://cran.r-project.org/package=rms) ·
[`ggplot2`](https://cran.r-project.org/package=ggplot2) ·
[`dplyr`](https://cran.r-project.org/package=dplyr)
