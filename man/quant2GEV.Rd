% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant2GEV.R
\name{quant2GEV}
\alias{quant2GEV}
\title{Find GEV Parameters from Given Quantiles}
\usage{
quant2GEV(q, p, lower.tail = TRUE, plot = FALSE, eps = 1e-07, trace = 0)
}
\arguments{
\item{q}{A numeric vector with length \eqn{3} containing distinct
values for the GEV quantiles.}

\item{p}{A numeric vector with length \eqn{3} containing distinct
values for the probabilities.}

\item{lower.tail}{Logical. If \code{TRUE} the values in \code{p}
are for the probability of exceedance, else they are for the
probability of non-exceedance.}

\item{plot}{Logical. If \code{TRUE} a simple plot will illustrate
the zero-finding used to find the GEV shape \eqn{\xi}.}

\item{eps}{A small numeric number used to decide when the shape
\eqn{xi} is close enough to zero.}

\item{trace}{Integer level of verbosity.}
}
\value{
A named numeric vector containing the GEV parameters.
}
\description{
Find the vector of the three GEV parameters corresponding to a
vector of three quantiles for three given distinct probabilities.
}
\details{
Given three distinct probabilities, there is a one-to-one
correspondance between the vector of the three corresponding
quantiles and the vector of the three GEV parameters so the vector
of quantiles can be used to re-parameterise the GEV distribution.
The quantile parameterisation can be preferred to define
informative priors based on expert knowledge, see Coles and Tawn.
}
\examples{

co <- quant2GEV(p = c(0.1, 0.01, 0.001),
                q = c(60, 80, 120), lower.tail = FALSE)

## check the result
nieve::qGEV(p = c(0.1, 0.01, 0.001),
            loc = co["loc"], scale = co["scale"], shape = co["shape"],
            lower.tail = FALSE)

}
\references{
Coles S. and Tawn J. (1996). A Bayesian Analysis of Extreme
Rainfall Data \emph{Appl. Statist.} 45 (4), pp. 463-478.
}
