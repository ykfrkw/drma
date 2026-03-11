# drma

**Dose-Response Meta-Analysis** — a `meta`/`netmeta`-style wrapper around
[`dosresmeta`](https://cran.r-project.org/package=dosresmeta) for
**one-stage dose-response meta-analysis**.

Pass arm-level data in long format and `drma` automatically computes OR,
RR, MD, or SMD before fitting a dose-response curve with restricted cubic
splines (or linear / log-linear / quadratic alternatives).

> **Statistical engine**: all modelling is performed by `dosresmeta`
> (Crippa & Orsini 2016).  `drma` provides the data-preparation,
> knot-resolution, plot, and utility layer on top of it.

## Installation

```r
remotes::install_github("ykfrkw/drma")
```

## Quick start

Data from Furukawa et al. 2022 (brexpiprazole dose-response MA) is bundled
with the package:

```r
library(drma)
data(brexpiprazole)

res <- drma(
  data    = brexpiprazole,
  studlab = study_id,
  dose    = dose,
  sm      = "OR",
  event   = n_responders,
  n       = n_arm,
  knots   = c(1, 2, 3)
)

plot(res, ylab = "Response (OR)", ylim = c(0.75, 2),
     ref_dose = 0, bubble = TRUE, rug = TRUE)
```

Column names can be **bare (unquoted)** or quoted strings — both work.

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
res_rcs  <- drma(data = df, ..., curve = "rcs",       knots = c(1, 2, 3))
res_log  <- drma(data = df, ..., curve = "log",       log_shift = 1)
res_lin  <- drma(data = df, ..., curve = "linear")
res_quad <- drma(data = df, ..., curve = "quadratic")
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
drma(data = df, ..., knots = "0.1-0.5-0.9")     # 10/50/90th percentile
drma(data = df, ..., knots = "0.25-0.50-0.75")  # 25/50/75th percentile
drma(data = df, ..., knots = c(1, 2, 3))         # actual dose values
drma(data = df, ..., knots = 3L)                 # 3 knots, auto-placed
```

---

## Standard workflow

```r
# 1. Fit
res <- drma(
  data    = brexpiprazole,
  studlab = study_id,
  dose    = dose,
  sm      = "OR",
  event   = n_responders,
  n       = n_arm,
  knots   = c(1, 2, 3)
)

# 2. Inspect
print(res)    # coefficients + study count
summary(res)  # full dosresmeta summary

# 3. Plot
plot(res,
     ylab     = "Response (OR)",
     ylim     = c(0.75, 2),
     ref_dose = 0,
     bubble   = TRUE,   # study-level data (area ∝ n)
     rug      = TRUE)   # tick marks for evaluated doses

# 4. Target doses
target_dose(res, p = c(0.5, 0.95, 1))   # ED50, ED95, ED100

# 5. Predictions at specific doses
predict_table(res,
              doses         = c(0, 1, 2, 3),
              baseline_prop = 0.183)   # convert OR to absolute risk

# 6. Pairwise league table
league_table(res, doses = c(0, 1, 2, 3))
```

---

## Precomputed log-effects (`sm = "precomputed"`)

When log-ORs (or log-RRs, MDs) are already computed — e.g. read directly
from a published table — pass them via `yi` and `sei`:

```r
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
```

---

## Overlaying multiple curves (sensitivity analyses)

```r
res_p  <- drma(data = brexpiprazole, ..., knots = c(1, 2, 3))
res_s1 <- drma(data = brexpiprazole, ..., knots = "0.1-0.5-0.9")
res_s2 <- drma(data = brexpiprazole, ..., knots = "0.25-0.50-0.75")
res_s3 <- drma(data = brexpiprazole, ..., curve = "log")

plot(res_p,  col = "black", lwd = 2,
     ylim = c(0.75, 2), ylab = "Response (OR)", ref_dose = 0)
lines(res_s1, col = "tomato",      lty = 2)
lines(res_s2, col = "steelblue",   lty = 3)
lines(res_s3, col = "forestgreen", lty = 4)
legend("topright",
       legend = c("Primary c(1,2,3)",
                  "S1: 10/50/90%",
                  "S2: 25/50/75%",
                  "S3: log-linear"),
       col = c("black", "tomato", "steelblue", "forestgreen"),
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
[`dplyr`](https://cran.r-project.org/package=dplyr)
