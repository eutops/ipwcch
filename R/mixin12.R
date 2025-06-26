#' Combine CoxPH Hazard Estimates with Exit Times
#'
#' This function merges a cumulative hazard function (from a `coxph` model) with
#' a vector of event or exit times (from a multi-state model), filling in missing
#' hazard values via last observation carried forward (LOCF).
#'
#' This is typically used to prepare cumulative hazard data for use in
#' transition probability calculations (e.g., with `msfit`).
#'
#' @param coxphData A `data.frame` containing:
#'   - `time`: time points at which cumulative hazards are estimated
#'   - `hazard`: corresponding cumulative hazard values
#' @param exitTime A numeric vector of exit or transition times for all observed states.
#'
#' @return A `data.frame` with:
#'   - `time`: unique sorted time points from both `coxphData$time` and `exitTime`
#'   - `hazard`: cumulative hazard values at each time point, with missing ones filled by LOCF
#'
#' @examples
#' # Example with dummy data
#' cox_data <- data.frame(time = c(1, 3, 5), hazard = c(0.1, 0.4, 0.9))
#' exit_times <- c(2, 4, 5, 6)
#' mixin12(cox_data, exit_times)
#'
#' @export mixin12
mixin12 <- function(coxphData, exitTime) {

  if (!all(c("time", "hazard") %in% names(coxphData))) {
    stop("coxphData must contain 'time' and 'hazard' columns.")
  }

  # Combine and sort all time points
  allTimes <- sort(unique(c(coxphData$time, exitTime)))

  # Initialize result with placeholder hazards
  result <- data.frame(time = allTimes, hazard = rep(NA_real_, length(allTimes)))

  # Match existing hazard values to time points
  hazard_match <- match(coxphData$time, result$time)
  result$hazard[hazard_match] <- coxphData$hazard

  # Fill in missing hazards via last observation carried forward (LOCF)
  for (i in seq_along(result$hazard)) {
    if (is.na(result$hazard[i])) {
      if (i == 1) {
        result$hazard[i] <- 0  # Assume cumulative hazard = 0 at time 0
      } else {
        result$hazard[i] <- result$hazard[i - 1]
      }
    }
  }

  added_times <- length(allTimes) - nrow(coxphData)
  message("Result: ",
          nrow(result),
          " rows. Added from exitTime: ",
          added_times)

  return(result)
}
