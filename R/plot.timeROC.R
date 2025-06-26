#' Plot ROC Curve from timeROC Output Using ggplot2
#'
#' This helper function creates a ROC plot from the output of the `timeROC` function,
#' using either Definition 1 or Definition 2 for handling competing risks.
#' Note that it only works for one time_index = 2, with times = c() a vector of length 2, and the first element = 0.
#'
#' @param timeROC_obj Output object from the `timeROC()` function.
#' @param ci_df A `data.frame` from `confint(timeROC_obj)` representing AUC confidence intervals.
#' @param definition Integer; either 1 or 2.
#'   \itemize{
#'     \item `1`: Competing events do not contribute to AUC (Definition 1).
#'     \item `2`: Competing risks treated as controls (Definition 2).
#'   }
#' @param time_index Integer index of the `times` value in the `timeROC` object to use for plotting (default = 2).
#' @param ... Additional arguments passed to `ggplot2::theme()` for plot customization.
#'
#' @return A `ggplot` object of the ROC curve with annotated AUC and 95% CI.
#'
#' @examples
#' \dontrun{
#' test <- timeROC(T = cohort$time, delta = cohort$status, marker = cohort$X,
#'                 cause = 1, times = c(0, 1826.25), weighting = "marginal")
#' result.ci <- confint(test, level = 0.95, n.sim = 1000) %>% as.data.frame()
#' plot_timeROC(test, result.ci, definition = 2)
#' }
#'
#' @import ggplot2
#' @export
plot_timeROC <- function(timeROC_obj, ci_df, definition = 2, time_index = 2, ...) {

  # Sanity checks
  if (!definition %in% c(1, 2)) stop("Only definitions 1 and 2 are supported.")

  # Extract FPR and TPR based on chosen definition
  fpr <- if (definition == 1) {
    timeROC_obj$FP_1[, time_index]
  } else {
    timeROC_obj$FP_2[, time_index]
  }

  tpr <- timeROC_obj$TP[, time_index]

  # Extract AUC and CI
  auc <- if (definition == 1) timeROC_obj$AUC_1[time_index] else timeROC_obj$AUC_2[time_index]
  ci_start_index <- if (definition == 1) 1 else 6  # Index shift for CI columns
  ci_lower <- unname(ci_df[ci_start_index] / 100)
  ci_upper <- unname(ci_df[ci_start_index + 1] / 100)

  # Create data frame for ggplot
  roc_df <- data.frame(
    FPR = fpr,
    TPR = tpr
  )

  # Construct the ggplot
  p <- ggplot(roc_df, aes(x = .data$FPR, y = .data$TPR)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = paste0("t = ", timeROC_obj$times[time_index],
                     " (Definition ", definition, ")"),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    annotate(
      "text",
      x = 1, y = 0.1,
      label = sprintf("AUC = %.3f (%.3f - %.3f)", auc, ci_lower, ci_upper),
      hjust = 1, size = 4.5
    ) +
    theme_classic() +
    theme(...)

  return(p)
}
