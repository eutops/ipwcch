
#' Plot a comparison of weighted ROC Curve using ggplot2
#'
#' Generates a ggplot of the two or more weighted ROC curves at a specified evaluation time `t_eval`,
#' using the output from `weightedRoc()` and `aucWithCI()`.
#'
#' @param x A nested list, with for each comparison a list object returned by `weightedRoc()` containing TPR, FPR, AUC, etc.
#' @param y A nested list, with for each comparison a list object returned by `aucWithCI()` containing AUC and its confidence interval.
#' @param label Group labels to be used in the plot
#' @param names Group names
#' @param colors Line colors for the different groups
#' @param ... Additional arguments passed to `ggplot2::theme()` for plot customization.
#'
#' @return A `ggplot` object showing the ROC curve with AUC and confidence interval annotated.
#' @examples
#' \dontrun{
#' #' roc_obj1 <- weightedRoc(data, t_eval = 5, marker = "biomarker")
#' auc_ci1  <- aucWithCI(data, t_eval = 5, marker = "biomarker")
#' roc_obj2 <- weightedRoc(data, t_eval = 5, marker = "biomarker")
#' auc_ci2  <- aucWithCI(data, t_eval = 5, marker = "biomarker")
#' roc_obj_list <- list(roc_obj1,roc_obj2)
#' auc_ci_list <- list(auc_ci1,auc_ci2)
#' plot_compareROC(roc_obj_list, auc_ci_list, t_eval = 5,
#'                names = c(cohort1, cohort2), colors = c("purple", "pink"))
#' }
#' @import ggplot2
#' @export
#'
plot_compareROC <- function(x, y, label, names, colors, ...) {

  roc_df_combined <- c()

  # Create a combined dataframe for each roc curve
  for (i in 1:length(x)){
    roc_df_combined[[i]] <- data.frame(
      FPR = x[[i]]$fpr,
      TPR = x[[i]]$tpr,
      name = names[i]
    )
  }

  roc_df_combined <- bind_rows(roc_df_combined)
  roc_df_combined$name <- factor(roc_df_combined$name, levels = names)

  # Construct the ggplot object
  p <- ggplot(roc_df_combined, aes(x = .data$FPR, y = .data$TPR, color = .data$name)) +
    geom_path() +  # ROC curve
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    scale_color_manual(values = colors) +
    annotate("label", x = 0.2, y = 0.95,
             label = label,
             color = "black",
             fill = "white",  # Background color of the box
             size = 4,
             label.size = 0.25) +
    theme_classic() +
    theme(...)

  # Add the AUC values
  y_position <- 0.4

  for (i in 1:length(y)){
    p <- p +
      annotate(
        "text",
        x = 1, y = y_position,
        label = sprintf("AUC = %.2f (%.2f - %.2f)", stats::median(y[[i]]$AUC_bootstrap), y[[i]]$CI_lower, y[[i]]$CI_upper),
        color = colors[i],
        hjust = 1, size = 3
      )
    y_position <- y_position - 0.1
  }

  return(p)
}
