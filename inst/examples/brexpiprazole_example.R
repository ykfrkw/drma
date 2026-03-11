# =============================================================================
# drma package — Real-data example
#
# Data source:
#   Furukawa Y, Oguro S, Obata S, Hamza T, Ostinelli EG, Kasai K.
#   Optimal dose of brexpiprazole for augmentation therapy of
#   antidepressant-refractory depression: A systematic review and
#   dose-effect meta-analysis.
#   Psychiatry Clin Neurosci. 2022;76(9):416-422.
#   doi:10.1111/pcn.13438
#
# Statistical method:
#   Crippa A, Discacciati A, Bottai M, Spiegelman D, Orsini N.
#   One-stage dose-response meta-analysis for aggregated data.
#   Stat Methods Med Res. 2019;28(5):1579-1596.
#   doi:10.1177/0962280218773122
# =============================================================================

library(drma)

data(brexpiprazole)   # built-in dataset

# --------------------------------------------------------------------------- #
# 1.  Efficacy — Response (OR)
#     sm = "OR": drma auto-computes log-OR from arm-level event counts
# --------------------------------------------------------------------------- #
res_e <- drma(
  data    = brexpiprazole,
  studlab = study_id,
  dose    = dose,
  sm      = "OR",
  event   = n_responders,
  n       = n_arm,
  knots   = c(1, 2, 3)        # primary knot placement (1 / 2 / 3 mg)
)

print(res_e)
summary(res_e)


# --------------------------------------------------------------------------- #
# 2.  Tolerability — Dropout due to AEs (precomputed log-OR)
#     sm = "precomputed": pass already-computed yi and sei
# --------------------------------------------------------------------------- #
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


# --------------------------------------------------------------------------- #
# 3.  Acceptability — Dropout for any reason (precomputed log-OR)
# --------------------------------------------------------------------------- #
res_a <- drma(
  data    = brexpiprazole,
  studlab = study_id,
  dose    = dose,
  sm      = "precomputed",
  yi      = acceptability_logor,
  sei     = acceptability_se,
  event   = n_dropout_any,
  n       = n_arm,
  knots   = c(1, 2, 3)
)


# --------------------------------------------------------------------------- #
# 4.  Plots — primary curves
# --------------------------------------------------------------------------- #
par(mfrow = c(1, 3))

plot(res_e,
     ylab     = "Response (OR)",
     ylim     = c(0.75, 2),
     xlab     = "Brexpiprazole (mg)",
     ref_dose = 0,
     bubble   = TRUE,
     rug      = TRUE)

plot(res_t,
     ylab     = "Dropout for AEs (OR)",
     ylim     = c(0.2, 5),
     xlab     = "Brexpiprazole (mg)",
     ref_dose = 0,
     bubble   = TRUE,
     rug      = TRUE)

plot(res_a,
     ylab     = "Dropout for any reason (OR)",
     ylim     = c(0.2, 5),
     xlab     = "Brexpiprazole (mg)",
     ref_dose = 0,
     bubble   = TRUE,
     rug      = TRUE)

par(mfrow = c(1, 1))


# --------------------------------------------------------------------------- #
# 5.  Sensitivity analyses — knot placement (efficacy)
# --------------------------------------------------------------------------- #

# S2.1: c(0.5, 1.5, 2.5)
res_s1 <- drma(
  data = brexpiprazole, studlab = study_id, dose = dose, sm = "OR",
  event = n_responders, n = n_arm,
  knots = c(0.5, 1.5, 2.5)
)

# S2.2: 25th / 50th / 75th percentile
res_s2 <- drma(
  data = brexpiprazole, studlab = study_id, dose = dose, sm = "OR",
  event = n_responders, n = n_arm,
  knots = "0.25-0.50-0.75"
)

# S2.3: 10th / 50th / 90th percentile
res_s3 <- drma(
  data = brexpiprazole, studlab = study_id, dose = dose, sm = "OR",
  event = n_responders, n = n_arm,
  knots = "0.1-0.5-0.9"
)

# S2.4: log-linear (alternative functional form)
res_s4 <- drma(
  data = brexpiprazole, studlab = study_id, dose = dose, sm = "OR",
  event = n_responders, n = n_arm,
  curve = "log", log_shift = 1
)

# Overlay all sensitivity curves on one plot
plot(res_e,
     col  = "black", lwd = 2,
     ylim = c(0.75, 2),
     ylab = "Response (OR)",
     xlab = "Brexpiprazole (mg)",
     ref_dose = 0)

lines(res_s1, col = "tomato",      lty = 2, lwd = 1.5)
lines(res_s2, col = "steelblue",   lty = 3, lwd = 1.5)
lines(res_s3, col = "forestgreen", lty = 4, lwd = 1.5)
lines(res_s4, col = "purple",      lty = 5, lwd = 1.5)

legend("topleft",
       legend = c("Primary c(1,2,3)",
                  "S2.1: c(0.5,1.5,2.5)",
                  "S2.2: 25/50/75th pct",
                  "S2.3: 10/50/90th pct",
                  "S2.4: log-linear"),
       col    = c("black", "tomato", "steelblue", "forestgreen", "purple"),
       lty    = 1:5,
       lwd    = 1.5,
       bty    = "n")


# --------------------------------------------------------------------------- #
# 6.  Predictions at clinically relevant doses
#     baseline_prop: baseline response rate ~18.3% (from the paper)
# --------------------------------------------------------------------------- #
predict_table(
  res_e,
  doses         = c(0, 0.5, 1, 1.5, 2, 3),
  ref_dose      = 0,
  baseline_prop = 0.183
)


# --------------------------------------------------------------------------- #
# 7.  Target doses (ED50 / ED95 / ED100)
# --------------------------------------------------------------------------- #
target_dose(res_e, p = c(0.5, 0.95, 1))


# --------------------------------------------------------------------------- #
# 8.  League table — pairwise dose comparisons
# --------------------------------------------------------------------------- #
league_table(res_e, doses = c(0, 1, 2, 3))

# Compare primary vs sensitivity knot specifications
lt <- league_table(
  res_e, res_s2, res_s3,
  doses  = c(0, 1, 2, 3),
  labels = c("Primary c(1,2,3)", "S2.2 25/50/75%", "S2.3 10/50/90%")
)
lt$`Primary c(1,2,3)`
lt$`S2.2 25/50/75%`
lt$`S2.3 10/50/90%`


# --------------------------------------------------------------------------- #
# 9.  VPC plot — variance partitioning
# --------------------------------------------------------------------------- #
vpc_plot(res_e, main = "Efficacy: VPC")
