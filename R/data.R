#' Brexpiprazole dose-response data (primary analysis)
#'
#' Arm-level data from the primary analysis of a dose-response meta-analysis
#' of brexpiprazole augmentation for treatment-resistant depression.  The
#' dataset contains 6 randomised controlled trials with a total of 14 arms
#' (placebo and 0.15 / 1 / 2 / 3 mg brexpiprazole).
#'
#' @format A `data.frame` with 14 rows and 12 columns:
#' \describe{
#'   \item{study_id}{Character. ClinicalTrials.gov identifier.}
#'   \item{dose}{Numeric. Brexpiprazole dose (mg/day); 0 = placebo.}
#'   \item{n_arm}{Integer. Patients randomised to the arm.}
#'   \item{n_responders}{Numeric. Responders (efficacy outcome).
#'     `NA` when not reported.}
#'   \item{n_dropout_ae}{Numeric. Dropouts due to adverse events
#'     (tolerability outcome).}
#'   \item{n_dropout_any}{Numeric. Dropouts for any reason
#'     (acceptability outcome).}
#'   \item{efficacy_logor}{Numeric. Pre-computed log-OR for efficacy
#'     vs placebo.  `0` for the reference (placebo) arm; `NA` when not
#'     estimable.}
#'   \item{efficacy_se}{Numeric. Standard error of `efficacy_logor`.}
#'   \item{tolerability_logor}{Numeric. Pre-computed log-OR for
#'     tolerability.}
#'   \item{tolerability_se}{Numeric. Standard error of
#'     `tolerability_logor`.}
#'   \item{acceptability_logor}{Numeric. Pre-computed log-OR for
#'     acceptability.}
#'   \item{acceptability_se}{Numeric. Standard error of
#'     `acceptability_logor`.}
#' }
#'
#' @source
#' Furukawa Y, Oguro S, Obata S, Hamza T, Ostinelli EG, Kasai K.
#' Optimal dose of brexpiprazole for augmentation therapy of
#' antidepressant-refractory depression: A systematic review and
#' dose-effect meta-analysis.
#' *Psychiatry Clin Neurosci*. 2022;76(9):416–422.
#' \doi{10.1111/pcn.13438}
#'
#' @examples
#' data(brexpiprazole)
#' head(brexpiprazole)
"brexpiprazole"
