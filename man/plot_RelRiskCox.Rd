% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.helpersRelRisk.R
\name{plot_RelRiskCox}
\alias{plot_RelRiskCox}
\title{Plot Cause-Specific Cumulative Hazard from Stratified Cox Model}
\usage{
plot_RelRiskCox(cohort, cox_fit, cols = default_cols)
}
\arguments{
\item{cohort}{A data frame containing the original cohort data, including a stratification variable named `strat`. The `strat` variable should be a factor.}

\item{cox_fit}{A Cox model object or an object with cumulative hazard estimates, standard errors, strata info, and optional p-values. It must contain elements `time`, `cumhaz`, `std.err`, `n.risk`, `n.event`, and `strata`. If `pvalue` is included, it will be annotated on the plot.}

\item{cols}{A character vector of colors used for plotting different strata. Defaults to `c("blue", "pink", "red")`. The number of colors should match the number of strata in `cox_fit`.}
}
\value{
A `ggplot2` object showing the stratified cumulative hazards with confidence bands and optional p-values.
}
\description{
This function generates a plot of the cause-specific cumulative hazard for stratified groups
based on a fitted Cox model. It includes pointwise 95% confidence intervals and optional
p-values, if available in the Cox model output.
}
