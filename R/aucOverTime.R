#' Calculate AUC Over Time with Confidence Intervals
#'
#' This function calculates the Area Under the Curve (AUC) for a given biomarker over specified time points,
#' along with the corresponding confidence intervals using nonparametric bootstrapping in weighted case-cohort designs with adminstrative censoring.
#'
#' This function depends on `weightedRoc()` and `aucWithCI()` to compute the AUC with CI for survival data with censoring and weights.
#'
#' @param dat A data frame with the same structure required by `weightedRoc()` and `aucWithCI()`, including:
#'   - `entry`: Entry time into risk set
#'   - `exit`: Exit time (event or censoring time)
#'   - `to`: Event type or status code
#'   - `weight`: Inverse probability or sampling weight
#'   - `marker`: A continuous risk score or biomarker column
#' @param times A numeric vector of time points at which AUC should be calculated.
#' @param marker Character string naming the column in `dat` containing the marker/risk score.
#' @param event_code Character or numeric value specifying the event of interest (default = "1").
#' @param eps Numeric buffer applied in the combining step of `weightedRoc()` to handle
#'   near-censoring at `tmax` (default = 0.0). See `weightedRoc()` for details.
#' @param n_boot Number of bootstrap resamples used to estimate the confidence interval (default = 1000).
#' @param seed Integer value used to set the random seed for reproducibility (default = 42).
#'
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item \code{t_eval}: The evaluated time points at which AUC was calculated.
#'   \item \code{AUC}: The calculated Area Under the Curve for each time point.
#'   \item \code{CI_lower}: The lower bound of the confidence interval for the AUC.
#'   \item \code{CI_upper}: The upper bound of the confidence interval for the AUC.
#' }
#'
#' @details This function iterates over the specified time points and calculates the AUC for each point,
#' along with the 95% confidence interval using bootstrapping.
#' If there are no individuals at risk at a particular time point, that time point is skipped, and a warning is issued.
#'
#' For weighted case-cohort designts, it is assumed `dat` includes appropriate inverse probability sample weights.
#' Weights should have been assigned, for this may use this package function `assignWeightsCch()`.
#'
#' @examples
#' # Create example data frame with required columns
#' dat <- data.frame(
#'   entry = rep(0, 30),                     # entry time (could be 0 for all)
#'   exit = rep(seq(1, 10, length.out = 30)),  # exit times evenly spaced
#'   to = sample(c(0, 1), 30, replace = TRUE), # event indicator 0/1
#'   weight = runif(30, 0.5, 1.5),           # random weights
#'   biomarker_column_name = rnorm(30)        # marker values
#' )
#' times <- seq(1, 5, by = 1)
#' marker <- "biomarker_column_name"
#' auc_results <- aucOverTime(dat = dat, times = times, marker = marker)
#' print(auc_results)
#'
#' @seealso \code{\link{weightedRoc}}, \code{\link{aucWithCI}}
#' @export aucOverTime
aucOverTime <- function(dat, times, marker, event_code = "1", eps = 0.0, n_boot = 1000, seed = 42) {

  # Ensure times are sorted
  times <- sort(times)

  # Preallocate results
  results <- data.frame(
    t_eval = numeric(length(times)),
    AUC = numeric(length(times)),
    CI_lower = numeric(length(times)),
    CI_upper = numeric(length(times)),
    stringsAsFactors = FALSE
  )

  tmax <- max(dat$exit)

  # Loop through each time point
  for (i in seq_along(times)) {

    t <- times[i]

    # Use buffer logic for near-tmax evaluations
    if (t >= tmax - eps) {
      # All individuals will be considered at risk
      has_data <- TRUE
    } else {
      # Check if there are any events or controls at risk
      is_case <- dat$exit <= t & dat$to == event_code
      is_control <- dat$entry <= t & dat$exit > t
      has_data <- any(is_case) && any(is_control)
    }

    if (!has_data) {
      warning(sprintf("Skipping t = %s: No events or controls at risk\n", t))
      results[i, ] <- c(t, NA, NA, NA)
      next
    }

    auc_ci <- aucWithCI(dat, t, marker, event_code, eps, n_boot, seed + i)
    results[i, ] <- c(t, stats::median(auc_ci$AUC_bootstrap), auc_ci$CI_lower, auc_ci$CI_upper)
  }

  return(results)
}
