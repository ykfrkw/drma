#' drma: Dose-Response Meta-Analysis with a meta/netmeta-Style Interface
#'
#' `drma` provides a user-friendly wrapper around the
#' [`dosresmeta`](https://cran.r-project.org/package=dosresmeta) package for
#' one-stage dose-response meta-analysis.  It is modelled on the interface of
#' [`meta`](https://cran.r-project.org/package=meta) and
#' [`netmeta`](https://cran.r-project.org/package=netmeta): pass arm-level
#' data in long format and `drma` automatically computes OR, RR, MD, or SMD
#' before fitting a dose-response curve with restricted cubic splines (or
#' linear / log-linear / quadratic alternatives).
#'
#' @section Core function:
#' * [drma()] — fit a dose-response model
#'
#' @section Post-fit utilities:
#' * [print.drma()] / [summary.drma()] — inspect the model
#' * [plot.drma()] / [lines.drma()] — draw and overlay curves
#' * [predict.drma()] / [predict_table()] — predictions at specific doses
#' * [target_dose()] — ED50 / ED95 / ED100
#' * [league_table()] — pairwise dose comparison matrix
#' * [vpc_plot()] — Variance Partitioning Coefficient plot
#'
#' @section Citation:
#' When using `drma`, please cite the underlying statistical method:
#'
#' Crippa A, Discacciati A, Bottai M, Spiegelman D, Orsini N.
#' One-stage dose-response meta-analysis for aggregated data.
#' *Stat Methods Med Res*. 2019;28(5):1579–1596.
#' \doi{10.1177/0962280218773122}
#'
#' and the `dosresmeta` package:
#'
#' Crippa A, Orsini N.
#' Dose-response meta-analysis of differences in means.
#' *BMC Med Res Methodol*. 2016;16:91.
#' \doi{10.1186/s12874-016-0189-0}
#'
#' @references
#' Crippa A, Discacciati A, Bottai M, Spiegelman D, Orsini N.
#' One-stage dose-response meta-analysis for aggregated data.
#' *Stat Methods Med Res*. 2019;28(5):1579–1596.
#' \doi{10.1177/0962280218773122}
#'
#' Crippa A, Orsini N.
#' Dose-response meta-analysis of differences in means.
#' *BMC Med Res Methodol*. 2016;16:91.
#' \doi{10.1186/s12874-016-0189-0}
#'
#' @keywords internal
"_PACKAGE"
