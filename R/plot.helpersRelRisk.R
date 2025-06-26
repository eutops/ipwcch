#' Plot Cause-Specific Cumulative Hazard from Stratified Cox Model
#'
#' This function generates a plot of the cause-specific cumulative hazard for stratified groups
#' based on a fitted Cox model. It includes pointwise 95% confidence intervals and optional
#' p-values, if available in the Cox model output.
#'
#' @param cohort A data frame containing the original cohort data, including a stratification variable named `strat`. The `strat` variable should be a factor.
#' @param cox_fit A Cox model object or an object with cumulative hazard estimates, standard errors, strata info, and optional p-values. It must contain elements `time`, `cumhaz`, `std.err`, `n.risk`, `n.event`, and `strata`. If `pvalue` is included, it will be annotated on the plot.
#' @param cols A character vector of colors used for plotting different strata. Defaults to `c("blue", "pink", "red")`. The number of colors should match the number of strata in `cox_fit`.
#'
#' @return A `ggplot2` object showing the stratified cumulative hazards with confidence bands and optional p-values.
#'
#' @import ggplot2
#' @export plot_RelRiskCox
plot_RelRiskCox <- function(cohort,cox_fit, cols = default_cols) {

  default_cols <- c("blue","pink","red")

  # Prepare dataframe
  df <- data.frame(
    time = cox_fit$time,
    estimate = cox_fit$cumhaz,
    conf.low = cox_fit$cumhaz - 1.96 * cox_fit$std.err,
    conf.high = cox_fit$cumhaz + 1.96 * cox_fit$std.err,
    n.risk = cox_fit$n.risk,
    n.event = cox_fit$n.event
  )

  # Extract and clean strata names
  strata_names <- gsub("strat=", "", rep(names(cox_fit$strata), cox_fit$strata))
  df$Stratum <- factor(strata_names, levels = levels(cohort$strat))

  # Format p-values
  pvalues <- c()
  if (!is.null(cox_fit$pvalue)) {
    for (i in 1:length(cox_fit$pvalue)) {
      if (cox_fit$pvalue[i] < 0.01) {
        pval <- paste('p = ',
                      format(
                        cox_fit$pvalue[i],
                        scientific = TRUE,
                        digits = 2
                      ),
                      sep = '')
      } else {
        pval <- paste('p = ', round(cox_fit$pvalue[i], digits = 2), sep = '')
      }

      pvalues <- c(pvalues,pval)
    }

  } else {
    pval <- NA  # Assign NA if p-value is not available
  }

  # Plot base
  p <- ggplot(data = df, aes(x = .data$time, y = .data$estimate, colour = .data$Stratum)) +
    geom_line(size = 0.7) +
    geom_ribbon(
      aes(ymin = .data$conf.low, ymax = .data$conf.high, fill = .data$Stratum),
      alpha = 0.1,
      linetype = "dotted"
    ) +
    ylab("Cause-specific cumulative hazard\n(semi-parametric)") +
    scale_colour_manual(values = cols[1:length(unique(cox_fit$strata))]) +
    scale_fill_manual(values = cols[1:length(unique(cox_fit$strata))]) +
    theme_classic() +
    theme(legend.title = element_blank(), aspect.ratio = 1)

  # Annotate p-values if provided

  if (!is.null(pvalues) && !is.null(strata_names)) {
    for (i in seq_along(pvalues)) {
      p <- p + annotate(
        "text",
        x = max(df$time) * 0.95,
        y = (0.2 - (i * 0.1)) * max(df$estimate),
        label = pvalues[i],
        colour = cols[i + 1],
        hjust = 1,
        size = 4
      )

    }
  }

  return(p)
}
