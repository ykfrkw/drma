# drma

**Dose-Response Meta-Analysis** — a `meta`/`netmeta`-style R package wrapping
[`dosresmeta`](https://github.com/alecri/dosresmeta).

Pass arm-level data in **long format** and `drma` automatically computes OR,
RR, MD, or SMD before fitting a restricted cubic spline dose-response curve.

## Installation

```r
# From GitHub (requires remotes)
remotes::install_github("YourGitHubUsername/drma")
```

## Quick start

```r
library(drma)

# Binary outcome — OR
res <- drma(
  data    = mydata,          # long format: one row per arm
  study   = "studyID",
  dose    = "dose",
  measure = "OR",
  events  = "n_events",
  n       = "n_total",
  knots   = c(1, 2, 3)       # actual dose values
)

print(res)
plot(res, ylab = "Response (OR)", ylim = c(0.5, 3))
target_dose(res, p = c(0.5, 0.95, 1))
predict_table(res, doses = c(0, 1, 2, 3), baseline_prop = 0.30)
```

## Knot specification

| `knots` value | Interpretation |
|---|---|
| `c(1, 2, 3)` | Actual dose values (any value ≥ 1) |
| `c(0.25, 0.5, 0.75)` | Quantile probabilities of observed doses |
| `3L` | 3 knots at 10th / 50th / 90th percentiles |

## Plotting options

```r
# Confidence band + bubble plot + rug marks
plot(res,
     bubble       = TRUE,   # study-level observations (area ∝ n)
     rug          = TRUE,   # tick marks for evaluated doses
     ref_dose     = 0)

# Overlay multiple curves
plot(res_primary, col = "black", ylim = c(0.5, 3))
lines(res_sensitivity, col = "red", lty = 2)
```

## League table

Pairwise comparison of any two doses, analogous to a network meta-analysis
league table:

```r
league_table(res_primary,
             doses = c(0, 1, 2, 3))

# Compare primary vs sensitivity
league_table(res_primary, res_sensitivity,
             doses  = c(0, 1, 2, 3),
             labels = c("Primary", "Sensitivity S1.1"))
```

## Pre-computed effects

If log-effects are already available (e.g. from an existing dataset):

```r
res <- drma(
  data    = df,
  study   = "id",
  dose    = "dose",
  measure = "precomputed",
  yi      = "logor",
  sei     = "se",
  events  = "cases",
  n       = "n",
  knots   = c(1.5, 3, 6)
)
```

## Key functions

| Function | Description |
|---|---|
| `drma()` | Fit dose-response model |
| `plot.drma()` | Dose-response curve with CI, bubbles, rug |
| `lines.drma()` | Add curve to existing plot |
| `predict_table()` | Predictions at specified doses |
| `target_dose()` | ED50, ED95, ED100 … |
| `league_table()` | Pairwise dose comparison matrix |
| `vpc_plot()` | Variance Partitioning Coefficient plot |

## Dependencies

- [`dosresmeta`](https://cran.r-project.org/package=dosresmeta)
- [`rms`](https://cran.r-project.org/package=rms)
- [`dplyr`](https://cran.r-project.org/package=dplyr)
