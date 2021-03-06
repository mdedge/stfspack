% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_functions.R
\name{many.outcome.sim}
\alias{many.outcome.sim}
\title{Simulate a multiple-testing scenario}
\usage{
many.outcome.sim(n, nsim, n.dv = 3, correl = 0.5)
}
\arguments{
\item{n}{The number of people measured in each group per simulation.}

\item{nsim}{The number of simulated samples to run.}

\item{n.dv}{The number of outcome measurements (i.e. dependent variables or DVs) taken per person.}

\item{correl}{The correlation of the (normally distributed) outcome measurements---all pairs of measurements have the same correlation.}
}
\value{
A matrix of p values with n.dv columns and nsim rows. Each row corresponds to a different sample, and the p values from the respective outcome measures are in distinct columns.
}
\description{
Simulates type I error inflation from multiple testing for a test comparing two group means. The groups are compared on n.dv outcome measurements, each mutually correlated at a user-specified level. A matrix of p values is returned with n.dv columns and nsim rows.
}
\examples{
sims <- many.outcome.sim(100, nsim = 1000, n.dv = 5, correl = .4)
mean(apply(sims, 1, min) < .05) #familywise error rate with no correction
mean(apply(sims, 1, min) < .05/5) #familywise error rate with Bonferonni correction for 5 comparisons
}
\keyword{data,}
\keyword{hacking}
\keyword{multiple}
\keyword{p}
\keyword{simulated}
\keyword{testing,}
