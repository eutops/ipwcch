#' Assign Dynamic Weights for Case-Cohort Analysis
#'
#' @md
#'
#' @description
#' Assigns dynamic time-dependent weights to a case-cohort dataset for survival analysis,
#' preparing it for use with Cox proportional hazards models. The function transforms the input
#' dataset by assigning weights based on case/subcohort status and event time, duplicating rows
#' where necessary to represent changes in weight just before an event.
#'
#' @param cohort_dat A `data.frame` with one row per individual. Required columns:
#'   - `id`: Unique individual identifier.
#'   - `time`: Time of event or censoring.
#'   - `status`: 1 if the event occurred, 0 if censored.
#'   - `subcohort`: 1 if in the randomly selected subcohort, 0 otherwise.
#' @param Psub A numeric value representing the sampling probability of being included in the subcohort (e.g., 0.1 for 10%).
#' @param eps A small positive numeric value used to offset event times slightly for weight transitions (default is `0.00005`).
#'
#' @details
#' The function creates an entry-exit format dataset with dynamic weights:
#'   - Individuals in the subcohort are given a weight of `1/Psub` before any event.
#'   - Cases in the subcohort have their weight changed to 1 just before the event.
#'   - Cases not in the subcohort are assigned a weight of 1 just before the event.
#'   - Controls not in the subcohort are not included (weight remains 0).
#'
#' The output dataset can be used directly with Cox proportional hazards models
#' that accept start-stop format (e.g., `coxph(Surv(entry, exit, event) ~ ...)`)
#'
#' @return A `data.frame` with additional columns:
#'   - `entry`, `exit`: Time intervals for dynamic risk modeling.
#'   - `from`, `to`: Transition states (used for multi-state models).
#'   - `weight`: Assigned analytical weights for each interval.
#'
#' @source
#' O'Brien et al. (2022). *The Case for Case-Cohort: An Applied Epidemiologist's Guide to Reframing Case-Cohort Studies to Improve Usability and Flexibility*.
#' \doi{10.1097/EDE.0000000000001469}
#'
#' @examples
#' \dontrun{
#'Psub <- 0.1  # Example: 10% subcohort sampling
#' cohort_weighted <- assignWeightsCch(my_data, Psub, eps = 1e-5)
#' }
#'
#' @export

assignWeightsCch <- function(cohort_dat, Psub, eps = 0.00005) {

  # Create new variables to deal with dynamic weighting
  cohort_dat$entry <- 0
  cohort_dat$exit <- cohort_dat$time
  cohort_dat$from <- 0
  cohort_dat$to <- cohort_dat$status

  # Until incidence, cases inside of the sub-cohort are controls
  cohort_dat[cohort_dat$subcohort == 1 & cohort_dat$to == 1, ]$to <- 0

  # Individuals not in the sub-cohort get a weight of zero
  cohort_dat$weight <- 0

  # Cases and controls inside the sub-cohort get 1/sampling probability
  cohort_dat[cohort_dat$subcohort == 1, ]$weight <- 1 / Psub

  # Cases in the subcohort change their weight from 1/Psub to 1 shortly before incidence; rows are doubled
  cohort_dat[cohort_dat$subcohort == 1 & cohort_dat$status == 1, ]$exit <- cohort_dat[cohort_dat$subcohort == 1 & cohort_dat$status == 1, ]$exit - eps
  cases.in.subcohort <- cohort_dat[cohort_dat$subcohort == 1 & cohort_dat$status == 1, ]
  cases.in.subcohort$weight <- 1
  cases.in.subcohort$entry <- cases.in.subcohort$exit
  cases.in.subcohort$exit <- cases.in.subcohort$entry + eps
  cases.in.subcohort$to <- 1

  # Cases not in the subcohort are assigned a weight of 1 shortly before incidence (enter CC); rows are doubled
  cases.no.subcohort <- cohort_dat[cohort_dat$subcohort == 0 & cohort_dat$to == 1, ]
  cases.no.subcohort$weight <- 1
  cases.no.subcohort$entry <- cases.no.subcohort$exit - eps

  # Combine data
  cohort_weighted <- rbind(cohort_dat, cases.in.subcohort, cases.no.subcohort)
  cohort_weighted <- cohort_weighted[cohort_weighted$weight != 0, ] # Cases outside of subcohort only enter shortly before incidence, doubled rows removed again
  cohort_weighted <- cohort_weighted[order(cohort_weighted$id), ]

  return(cohort_weighted)

}
