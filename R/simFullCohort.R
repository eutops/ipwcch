#' Simulate a Full Cohort with Competing Risks
#'
#' This function generates a synthetic dataset simulating a full cohort with competing risks, incorporating
#' a continuous covariate and allowing for either administrative or informative censoring.
#'
#' @param n Integer. Number of individuals to simulate. Default is 1000.
#' @param lambda1 Numeric. Baseline hazard rate for cause 1. Default is 0.0008.
#' @param lambda2 Numeric. Baseline hazard rate for cause 2. Default is 0.0003.
#' @param meanX Numeric. Mean of the normally distributed covariate \code{X}. Default is 0.
#' @param sdX Numeric. Standard deviation of the covariate \code{X}. Default is 1.
#' @param censoring Character. Type of censoring: either \code{"administrative"} or \code{"informative"}.
#' @param Tmax Numeric. Maximum follow-up time for administrative censoring. Must be provided if \code{censoring = "administrative"}.
#' @param censor.rate Numeric. Rate parameter for exponential censoring times under informative censoring. Must be provided if \code{censoring = "informative"}.
#'
#' @details
#' This function simulates two competing event times:
#' \itemize{
#'   \item \code{T1} from an exponential distribution with a rate that increases with covariate \code{X}.
#'   \item \code{T2} from an exponential distribution with a rate that decreases with covariate \code{X}.
#' }
#' The observed event time \code{T} is the minimum of the two, and the event cause is determined accordingly.
#'
#' Depending on the \code{censoring} argument, censoring times are applied either:
#' \itemize{
#'   \item \strong{Administratively:} All individuals are followed up to a common time point \code{Tmax}.
#'   \item \strong{Informatively:} Each individual is assigned a censoring time drawn from an exponential distribution with rate \code{censor.rate}.
#' }
#' The observed time is the minimum of the event and censoring time, and status indicates the cause of the event (1 or 2) or 0 for censoring.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{id}{Individual ID}
#'   \item{time}{Observed time (either event or censoring time)}
#'   \item{status}{Event indicator: 0 = censored, 1 = cause 1, 2 = cause 2}
#'   \item{X}{Covariate value}
#' }
#'
#' @examples
#' # Administrative censoring example
#' df_admin <- simFullCohort(n = 1000, Tmax = 10, censoring = "administrative")
#'
#' # Informative censoring example
#' df_info <- simFullCohort(n = 1000, censor.rate = 0.001, censoring = "informative")
#'
#' @export simFullCohort

simFullCohort <- function(n = 1000,
                          lambda1 = 0.0008,
                          lambda2 = 0.0003,
                          meanX = 0,
                          sdX = 1,
                          censoring = c("administrative", "informative"),
                          Tmax = NULL,
                          censor.rate = NULL) {

  censoring <- match.arg(censoring)

  # Validate input depending on censoring type
  if (censoring == "administrative" && is.null(Tmax)) {
    stop("For administrative censoring, please specify Tmax.")
  }
  if (censoring == "informative" && is.null(censor.rate)) {
    stop("For informative censoring, please specify censor.rate.")
  }

  # Generate continuous covariate
  X <- stats::rnorm(n, mean = meanX, sd = sdX)

  # Cause-specific event times
  T1 <- stats::rexp(n, rate = lambda1 * exp(0.7 * X))
  T2 <- stats::rexp(n, rate = lambda2 * exp(-0.3 * X))
  T <- pmin(T1, T2)
  cause <- ifelse(T1 < T2, 1, 2)

  # Censoring logic
  if (censoring == "informative") {
    C <- stats::rexp(n, rate = censor.rate)
    time <- pmin(T, C)
    status <- ifelse(T <= C, cause, 0)  # 0 = censored
  } else {  # Administrative
    time <- pmin(T, Tmax)
    status <- ifelse(T <= Tmax, cause, 0)  # 0 = censored at Tmax
  }

  data <- data.frame(
    id = 1:n,
    time = time,
    status = status,
    X = X
  )

  return(data)
}
