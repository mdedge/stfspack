% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_functions.R
\name{mat.samps}
\alias{mat.samps}
\title{Generate a matrix of samples from a chosen distribution.}
\usage{
mat.samps(n, nsim = 10000, rx = rnorm, ...)
}
\arguments{
\item{n}{The number of independent observations to include in each sample (i.e. each row of the output matrix).}

\item{nsim}{The number of samples to generate (i.e. the number of rows in the output matrix).}

\item{rx}{A function generating random data (rnorm by default).}

\item{...}{Additional arguments to rx.}
}
\value{
A matrix of independent random numbers with nsim rows and n columns.
}
\description{
Draws random data and formats them into a matrix, where each row contains a sample. By default, the random data are drawn form a normal distribution, but the user can supply an alternative.
}
\examples{
mat.samps(3, 5) #data from standard normal
mat.samps(3, 5, mean = 10, sd = 2) #normal with expectation 10 and sd 2.
mat.samps(3, 5, rx = rexp, .1) #exponential with rate .1
mat.samps(3, 5, rx = rnorm.contam, contam.p = .2) #standard normal contaminated with 20\% outliers
}
\keyword{matrix}
\keyword{samples,}
\keyword{simulation,}
