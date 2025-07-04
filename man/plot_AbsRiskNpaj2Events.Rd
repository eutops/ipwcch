% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.helpersAbsRisk.R
\name{plot_AbsRiskNpaj2Events}
\alias{plot_AbsRiskNpaj2Events}
\title{Plot Non-Parametric Transition Probabilities (Aalen-Johansen)}
\usage{
plot_AbsRiskNpaj2Events(transprob, cols = default_cols)
}
\arguments{
\item{transprob}{Output from a multistate non-parametric estimator such as `etm::etm`.}

\item{cols}{Vector of colors for each transition. Default is `c("black", "darkgrey")`.}
}
\value{
A `ggplot` object.
}
\description{
Generates a plot of non-parametric Aalen-Johansen estimates of transition probabilities over time.
}
