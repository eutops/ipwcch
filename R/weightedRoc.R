#' Compute a Weighted ROC Curve for Time-to-Event Data
#'
#' Calculates a weighted ROC curve for a specified evaluation time (`t_eval`)
#' using time-to-event data with possible right-censoring and precalculated weights.
#' Designed for use in case-cohort where dynamic inverse sampling probability weights are applied.
#' Double-weighting (IPWC) could be done by multiplying IPW * IPCW weights, enter these for dat$weight.
#' When only IPW weights are calculated (administrative censoring), the size of the weights will have no effect
#' on the ROC shape nor the AUC, as such weights are constant for the case and control group at a given time point.
#'
#' @param dat A data frame containing survival data with the following required columns:
#'   - `entry`: Entry time (e.g., time at risk begins)
#'   - `exit`: Exit time (e.g., time of event or censoring)
#'   - `to`: Event type or status code
#'   - `weight`: Sampling or inverse-probability weight
#'   - `marker`: A continuous biomarker or risk score used for ROC computation
#' @param t_eval Numeric value indicating the evaluation time point at which the ROC curve is computed.
#' @param marker Character string specifying the name of the column in `dat` containing the risk score or marker.
#' @param event_code Character or numeric value indicating the event of interest (default = "1").
#' @param eps Numeric value indicating the epsilon time buffer applied near the maximum follow-up time (`tmax`)
#'   to avoid excluding individuals due to minor censoring timing variations (default = 0.0).
#'
#' @return A list containing:
#'   - `tpr`: True positive rates at each threshold
#'   - `fpr`: False positive rates at each threshold
#'   - `thresholds`: Sorted marker values corresponding to ROC points
#'   - `auc`: Area under the ROC curve
#'   - `call`: Function call
#'   - `t_eval`: Evaluation time used
#'
#' @details
#' If `t_eval` is within `eps` of the maximum follow-up time (`tmax`), all observations are included
#' to ensure stability of the risk set and avoid excluding individuals who were censored slightly before `tmax`.
#'
#' Ties in the marker variable are randomly broken.
#'
#' @examples
#' # Example usage with a survival dataset:
#' # weightedRoc(dat = mydata, t_eval = 365, marker = "risk_score", eps = 5)
#'
#' @export
weightedRoc <- function(dat, t_eval, marker, event_code = "1", eps = 0.0) {

  # CASE: Event of interest happened before or at t_eval
  is_case <- dat$entry <= t_eval & dat$exit <= t_eval & dat$to == event_code

  # CONTROL: Did not have event of interest before t_eval
  # Includes censoring or competing events
  #is_control <- !(dat$exit <= t_eval & dat$to == event_code)
  is_control <- dat$entry <= t_eval & dat$exit > t_eval

  # Combine
  tmax <- max(dat$exit)

  # Combine
  if (t_eval >= tmax - eps) {
    message(sprintf("t_eval = %s is within eps of tmax = %s; including all observations to avoid excluding near-at-risk individuals.", t_eval, tmax))
    # Include all observations to avoid empty at-risk set, but remove duplicated rows for cases inside subcohort
    # Buffer applied to avoid excluding individuals due to small variations in follow-up times near maximum follow-up
    at_risk <- dat[!(dat$status == 1 & dat$to == 0), ]
  } else {
    at_risk <- dat[is_case | is_control, ]
  }

  # Check for duplicated IDs
  if (any(duplicated(at_risk$id))) {
    warning("There are duplicated 'id' values in the 'at_risk' data.")
  }

  # Define event: only the event of interest counts as 1
  event <- ifelse(at_risk$to %in% event_code, 1, 0)

  # Marker and weight
  marker_vals <- at_risk[[marker]]
  weights <- at_risk$weight

  # Order by marker descending (high marker = high risk)
  ord <- order(-marker_vals, stats::runif(length(marker_vals)))  # Tie-breaking
  event <- event[ord]
  weights <- weights[ord]
  marker_vals <- marker_vals[ord]

  # Cumulative weighted sums
  tp_weights <- cumsum(weights * event)
  fp_weights <- cumsum(weights * (1 - event))

  # Normalize by total weights
  total_event_weight <- sum(weights * event)
  total_nonevent_weight <- sum(weights * (1 - event))

  # Avoid div-by-zero
  if (total_event_weight == 0 || total_nonevent_weight == 0) {
    stop("No events or non-events at time t_eval. Check your dat or t_eval.")
  }

  tpr <- tp_weights / total_event_weight
  fpr <- fp_weights / total_nonevent_weight

  auc = sum(diff(fpr) * (utils::head(tpr, -1) + utils::tail(tpr, -1)) / 2) # Interpolation-based AUC

  roc_obj <- list(
    tpr = tpr,
    fpr = fpr,
    thresholds = marker_vals,
    auc = auc,
    call = match.call(),
    t_eval = t_eval
  )

  return(roc_obj)
}
