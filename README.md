# drma

**Dose-Response Meta-Analysis** — a `meta`/`netmeta`-style R package wrapping
[`dosresmeta`](https://github.com/alecri/dosresmeta).

Pass arm-level data in **long format** and `drma` automatically computes OR,
RR, MD, or SMD before fitting a dose-response curve with restricted cubic
splines (or linear / log-linear / quadratic alternatives).

## Installation

```r
remotes::install_github("ykfrkw/drma")
```

## Quick start

```r
library(drma)

res <- drma(
  data    = mydata,      # long format: one row per arm
  studlab = "studyID",
  dose    = "dose",
  sm      = "OR",
  event   = "n_events",
  n       = "n_total",
  knots   = c(1, 2, 3)
)
```

## Argument overview

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
| `knots` / `knots_type` | RCS knot specification (see below) |

---

## Curve types (`curve =`)

| Value | Formula | Notes |
|---|---|---|
| `"rcs"` (default) | restricted cubic spline | knots required |
| `"linear"` | `effect ~ dose` | |
| `"log"` | `effect ~ log(dose + log_shift)` | `log_shift = 1` by default |
| `"quadratic"` | `effect ~ dose + dose²` | |

```r
# Restricted cubic spline (default)
res_rcs  <- drma(..., curve = "rcs",       knots = c(1, 2, 3))

# Log-linear (safe for dose = 0 arms with log_shift = 1)
res_log  <- drma(..., curve = "log",       log_shift = 1)

# Linear
res_lin  <- drma(..., curve = "linear")

# Quadratic
res_quad <- drma(..., curve = "quadratic")
```

---

## Knot specification for RCS

Three equivalent styles — choose whichever is most natural:

```r
# 1. Actual dose values (any value >= 1 triggers this)
drma(..., knots = c(1, 2, 3))

# 2. Quantile probabilities (all values in (0, 1))
drma(..., knots = c(0.1, 0.5, 0.9))   # 10th / 50th / 90th percentile
drma(..., knots = c(0.25, 0.5, 0.75)) # 25th / 50th / 75th percentile

# 3. Auto-placement: single integer = number of knots
drma(..., knots = 3L)   # 3 knots at 10th / 50th / 90th percentile

# Force interpretation when ambiguous (e.g. 0.5 as a dose value, not quantile)
drma(..., knots = c(0.5, 1, 2), knots_type = "values")
drma(..., knots = c(0.1, 0.5, 0.9), knots_type = "quantile")
```

---

## Standard workflow

```r
# 1. Fit
res <- drma(
  data    = df,
  studlab = "studyID",
  dose    = "dose",
  sm      = "OR",
  event   = "N_responders_arm",
  n       = "N_arm",
  knots   = c(1, 2, 3)
)

# 2. Inspect
print(res)    # model coefficients + study count
summary(res)  # full dosresmeta summary

# 3. Plot
plot(res,
     ylab       = "Response (OR)",
     ylim       = c(0.5, 3),
     ref_dose   = 0,
     bubble     = TRUE,   # study-level data (area ∝ n)
     rug        = TRUE)   # tick marks for evaluated doses

# 4. Target doses
target_dose(res, p = c(0.5, 0.95, 1))   # ED50, ED95, ED100

# 5. Predictions at specific doses (with absolute risk conversion)
predict_table(res,
              doses         = c(0, 1, 2, 3),
              baseline_prop = 0.30)    # 30% baseline response rate

# 6. Pairwise league table
league_table(res, doses = c(0, 1, 2, 3))
```

---

## Continuous outcomes (MD / SMD)

```r
res_md <- drma(
  data    = df,
  studlab = "studyID",
  dose    = "dose",
  sm      = "MD",
  mean    = "mean_arm",
  sd      = "sd_arm",
  n       = "n_arm",
  knots   = c(0.25, 0.5, 0.75)   # percentile-based
)
plot(res_md, ylab = "Mean Difference")
```

---

## Pre-computed log-effects

```r
res_pre <- drma(
  data    = df,
  studlab = "id",
  dose    = "dose",
  sm      = "precomputed",
  yi      = "logor",
  sei     = "se",
  event   = "cases",      # needed for within-study covariance
  n       = "n",
  knots   = c(1.5, 3, 6)
)
```

---

## Overlaying multiple curves

```r
# Primary analysis
res_p  <- drma(data = df, ..., knots = c(1, 2, 3))

# Sensitivity: different knot locations
res_s1 <- drma(data = df, ..., knots = c(0.5, 1.5, 2.5))
res_s2 <- drma(data = df, ..., knots = c(0.25, 0.5, 0.75))

# Sensitivity: log-linear curve
res_s3 <- drma(data = df, ..., curve = "log")

plot(res_p,  col = "black", lwd = 2, ylim = c(0.5, 3),
     ylab = "Response (OR)")
lines(res_s1, col = "tomato",      lty = 2)
lines(res_s2, col = "steelblue",   lty = 3)
lines(res_s3, col = "forestgreen", lty = 4)
legend("topright",
       legend = c("Primary (1,2,3 mg)", "S1 (0.5,1.5,2.5 mg)",
                  "S2 (25/50/75%)", "S3 (log-linear)"),
       col    = c("black","tomato","steelblue","forestgreen"),
       lty    = 1:4, bty = "n")
```

---

## League table

Pairwise comparison at specified doses — analogous to an NMA league table:

```r
# Single model
league_table(res, doses = c(0, 1, 2, 3))
#         0     1               2               3
# 0  "0"   "1.00 (1.00, 1.00)" ...             ...
# 1  "..."  "1"                "..."           "..."
# 2  ...
# 3  ...

# Compare primary vs sensitivity analysis
lt <- league_table(res_p, res_s1,
                   doses  = c(0, 1, 2, 3),
                   labels = c("Primary", "Sensitivity S1"))
lt$Primary
lt$`Sensitivity S1`
```

---

## Function reference

| Function | Description |
|---|---|
| `drma()` | Fit dose-response model |
| `print.drma()` | Coefficients and model info |
| `summary.drma()` | Full model summary |
| `plot.drma()` | Dose-response curve (CI, bubble, rug) |
| `lines.drma()` | Add curve to existing plot |
| `predict.drma()` | Raw predictions |
| `predict_table()` | Predictions at specific doses |
| `target_dose()` | ED50 / ED95 / ED100 |
| `league_table()` | Pairwise dose comparison matrix |
| `vpc_plot()` | Variance Partitioning Coefficient plot |

## Dependencies

[`dosresmeta`](https://cran.r-project.org/package=dosresmeta) ·
[`rms`](https://cran.r-project.org/package=rms) ·
[`dplyr`](https://cran.r-project.org/package=dplyr)
