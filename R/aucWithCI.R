#' Compute AUC with Confidence Interval via Bootstrapping
#'
#' Estimates the area under the weighted ROC curve (AUC) at a specified time point `t_eval`
#' along with a 95% confidence interval using nonparametric bootstrapping.
#'
#' This function depends on `weightedRoc()` to compute the AUC for survival data with censoring and weights,
#' and supports designs such as case-cohort via optional weighting.
#'
#' @param dat A data frame with the same structure required by `weightedRoc()`, including:
#'   - `entry`: Entry time into risk set
#'   - `exit`: Exit time (event or censoring time)
#'   - `to`: Event type or status code
#'   - `weight`: Inverse probability or sampling weight
#'   - `marker`: A continuous risk score or biomarker column
#' @param t_eval Numeric value for the time point at which AUC is evaluated.
#' @param marker Character string naming the column in `dat` containing the marker/risk score.
#' @param event_code Character or numeric value specifying the event of interest (default = "1").
#' @param eps Numeric buffer applied in the combining step of `weightedRoc()` to handle
#'   near-censoring at `tmax` (default = 0.0). See `weightedRoc()` for details.
#' @param n_boot Number of bootstrap resamples used to estimate the confidence interval (default = 1000).
#' @param seed Integer value used to set the random seed for reproducibility (default = 42).
#'
#' @return A list with following components:
#' \itemize{
#'   \item AUC: The observed AUC at time `t_eval` on the original data.
#'   \item CI_lower: Lower bound of the 95% bootstrap percentile confidence interval.
#'   \item CI_upper: Upper bound of the 95% bootstrap percentile confidence interval.
#'   \item AUC_bootstrap: Vector of AUC values from bootstrap replicates (after removing NAs).
#'   \item t_eval: The evaluation time at which the AUC was computed.
#'   \item n_boot: Number of successful bootstrap samples (after removing any with errors).
#' }
#'
#' @details
#' Bootstrap resampling is done by drawing rows with replacement from `dat`.
#' For weighted case-cohort designts, this assumes `dat` includes appropriate inverse probability sample weights.
#' Weights should have been assigned, for this may use this package function `assignWeightsCch()`.
#'
#' If a bootstrap sample results in no events or non-events (causing AUC computation to fail),
#' the result is skipped (NA removed).
#'
#' @examples
#' \dontrun{
#' #' aucWithCI(dat = mydata, t_eval = 365, marker = "risk_score", eps = 5)
#' }
#' @seealso \code{\link{weightedRoc}}
#' @export aucWithCI
aucWithCI <- function(dat, t_eval, marker, event_code = "1", eps = 0.0, n_boot = 1000, seed = 42) {
  set.seed(seed)

  # Step 1: AUC on original data
  roc_obj <- weightedRoc(dat, t_eval, marker, event_code)
  auc_obs <- roc_obj$auc

  # Step 2: Bootstrap
  auc_boot <- numeric(n_boot)

  for (i in seq_len(n_boot)) {
    # Resample: Sample with replacement within subcohort or all "at risk"
    boot_sample <- dat[sample(1:nrow(dat), replace = TRUE), ]

    # Try-catch to handle rare edge cases (e.g., 0 events in a resample)
    auc_boot[i] <- tryCatch({
      weightedRoc(boot_sample, t_eval, marker, event_code, eps)$auc
    }, error = function(e) NA)
  }

  # Remove NA (e.g., from failed resamples)
  auc_boot <- auc_boot[!is.na(auc_boot)]

  # Step 3: CI from bootstrap percentiles
  ci <- stats::quantile(auc_boot, probs = c(0.025, 0.975))

  return(list(
    AUC = auc_obs,
    CI_lower = ci[1],
    CI_upper = ci[2],
    AUC_bootstrap = auc_boot,
    t_eval = t_eval,
    n_boot = length(auc_boot)
  ))
}
