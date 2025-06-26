
#' Plot Weighted ROC Curve using ggplot2
#'
#' Generates a ggplot of the weighted ROC curve at a specified evaluation time `t_eval`,
#' using the output from `weightedRoc()` and `aucWithCI()`.
#'
#' @param x A list object returned by `weightedRoc()` containing TPR, FPR, AUC, etc.
#' @param y A list object returned by `aucWithCI()` containing AUC and its confidence interval.
#' @param t_eval Numeric value for the time point at which AUC is evaluated.
#' @param ... Additional arguments passed to `ggplot2::theme()` for plot customization.
#'
#' @return A `ggplot` object showing the ROC curve with AUC and confidence interval annotated.
#' @examples
#' \dontrun{
#' roc_obj <- weightedRoc(data, t_eval = 5, marker = "biomarker")
#' auc_ci  <- aucWithCI(data, t_eval = 5, marker = "biomarker")
#' plot_ROC(roc_obj, auc_ci, 5)
#' }
#'
#' @import ggplot2
#' @export
#'
plot_ROC <- function(x, y, t_eval, ...) {

  # Create data frame for ROC curve
  roc_df <- data.frame(
    FPR = x$fpr,
    TPR = x$tpr
  )

  # Construct the ggplot object
  p <- ggplot(roc_df, aes(x = .data$FPR, y = .data$TPR)) +
    geom_line() +  # ROC curve
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(
      title = paste0("t = ", x$t_eval),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    annotate(
      "text",
      x = 1, y = 0.1,
      label = sprintf("AUC = %.3f (%.3f - %.3f)", stats::median(y$AUC_bootstrap), y$CI_lower, y$CI_upper),
      hjust = 1, size = 4.5
    ) +
    theme_classic() +
    theme(...)

  return(p)
}
