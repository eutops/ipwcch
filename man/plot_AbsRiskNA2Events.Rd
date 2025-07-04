% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.helpersAbsRisk.R
\name{plot_AbsRiskNA2Events}
\alias{plot_AbsRiskNA2Events}
\title{Plot Cause-Specific Cumulative Hazards (Non-Parametric)}
\usage{
plot_AbsRiskNA2Events(comhaz.out, cols = default_cols)
}
\arguments{
\item{comhaz.out}{List output of cumulative hazard estimates from a non-parametric method (e.g., `prodlim` or `survfit`).}

\item{cols}{Vector of colors for each transition. Default is `c("black", "darkgrey")`.}
}
\value{
A `ggplot` object.
}
\description{
Generates a non-parametric plot of cause-specific cumulative hazards with pointwise 95% confidence intervals.
}
