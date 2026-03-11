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

Column names can be specified **with or without quotes** (meta/netmeta style):

```r
library(drma)

# Quoted strings (classic)
res <- drma(
  data    = mydata,
  studlab = "studyID",
  dose    = "dose",
  sm      = "OR",
  event   = "n_events",
  n       = "n_total",
  knots   = c(1, 2, 3)
)

# Bare names (no quotes — meta/netmeta style)
res <- drma(
  data    = mydata,
  studlab = studyID,
  dose    = dose,
  sm      = "OR",
  event   = n_events,
  n       = n_total,
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
| `knots` | RCS knot specification (see below) |

---

## Curve types (`curve =`)

| Value | Formula | Notes |
|---|---|---|
| `"rcs"` (default) | restricted cubic spline | knots required |
| `"linear"` | `effect ~ dose` | |
| `"log"` | `effect ~ log(dose + log_shift)` | `log_shift = 1` by default |
| `"quadratic"` | `effect ~ dose + dose²` | |

```r
res_rcs  <- drma(..., curve = "rcs",       knots = c(1, 2, 3))
res_log  <- drma(..., curve = "log",       log_shift = 1)
res_lin  <- drma(..., curve = "linear")
res_quad <- drma(..., curve = "quadratic")
```

---

## Knot specification for RCS

| `knots =` | Interpretation |
|---|---|
| `"0.1-0.5-0.9"` | 10th / 50th / 90th percentile of observed doses |
| `"0.25-0.50-0.75"` | 25th / 50th / 75th percentile (default) |
| `c(1, 2, 3)` | Actual dose values |
| `3L` | 3 knots auto-placed at evenly-spaced quantiles |

```r
# Percentile string (recommended for sensitivity analyses)
drma(..., knots = "0.1-0.5-0.9")
drma(..., knots = "0.25-0.50-0.75")

# Actual dose values
drma(..., knots = c(1, 2, 3))
drma(..., knots = c(1.5, 3, 6))

# Auto-placement (3 knots at 10/50/90th percentile)
drma(..., knots = 3L)
```

---

## Standard workflow

```r
# 1. Fit
res <- drma(
  data    = df,
  studlab = studyID,           # bare name OK
  dose    = dose,
  sm      = "OR",
  event   = N_responders_arm,
  n       = N_arm,
  knots   = c(1, 2, 3)
)

# 2. Inspect
print(res)    # coefficients + study count
summary(res)  # full dosresmeta summary

# 3. Plot
plot(res,
     ylab     = "Response (OR)",
     ylim     = c(0.5, 3),
     ref_dose = 0,
     bubble   = TRUE,   # study-level data (area ∝ n)
     rug      = TRUE)   # tick marks for evaluated doses

# 4. Target doses
target_dose(res, p = c(0.5, 0.95, 1))   # ED50, ED95, ED100

# 5. Predictions at specific doses
predict_table(res,
              doses         = c(0, 1, 2, 3),
              baseline_prop = 0.30)    # convert OR to absolute risk (30% baseline)

# 6. Pairwise league table
league_table(res, doses = c(0, 1, 2, 3))
```

---

## Continuous outcomes (MD / SMD)

```r
res_md <- drma(
  data    = df,
  studlab = studyID,
  dose    = dose,
  sm      = "MD",
  mean    = mean_arm,
  sd      = sd_arm,
  n       = n_arm,
  knots   = "0.25-0.50-0.75"
)
plot(res_md, ylab = "Mean Difference")
```

---

## Pre-computed log-effects

```r
res_pre <- drma(
  data    = df,
  studlab = id,
  dose    = dose,
  sm      = "precomputed",
  yi      = logor,
  sei     = se,
  event   = cases,   # needed for within-study covariance
  n       = n,
  knots   = c(1.5, 3, 6)
)
```

---

## Overlaying multiple curves (sensitivity analyses)

```r
res_p  <- drma(data = df, ..., knots = c(1, 2, 3))       # primary
res_s1 <- drma(data = df, ..., knots = "0.1-0.5-0.9")    # S1: 10/50/90%
res_s2 <- drma(data = df, ..., knots = "0.25-0.50-0.75") # S2: 25/50/75%
res_s3 <- drma(data = df, ..., curve = "log")             # S3: log-linear

plot(res_p,  col = "black", lwd = 2,
     ylim = c(0.5, 3), ylab = "Response (OR)")
lines(res_s1, col = "tomato",      lty = 2)
lines(res_s2, col = "steelblue",   lty = 3)
lines(res_s3, col = "forestgreen", lty = 4)
legend("topright",
       legend = c("Primary c(1,2,3)",
                  "S1: 10/50/90%",
                  "S2: 25/50/75%",
                  "S3: log-linear"),
       col = c("black","tomato","steelblue","forestgreen"),
       lty = 1:4, bty = "n")
```

---

## League table

```r
# Single model — pairwise comparison at each dose pair
league_table(res, doses = c(0, 1, 2, 3))

# Compare primary vs sensitivity
lt <- league_table(res_p, res_s1,
                   doses  = c(0, 1, 2, 3),
                   labels = c("Primary", "S1 10/50/90%"))
lt$Primary
lt$`S1 10/50/90%`
```

---

## Notes

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
