% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.helpersAbsRisk.R
\name{makeCox2eventsDf}
\alias{makeCox2eventsDf}
\title{Prepare Cox Model Output for Plotting}
\usage{
makeCox2eventsDf(cox.out1, cox.out2)
}
\arguments{
\item{cox.out1}{Output from a Cox model (e.g., `survfit`) for event 1.}

\item{cox.out2}{Output from a Cox model for event 2.}
}
\value{
A `data.frame` containing combined cumulative hazard data for both transitions.
}
\description{
Creates a tidy data frame from two Cox model outputs, including cumulative hazard estimates, standard errors, and confidence intervals.
}
