% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poissonGPBayes.R
\name{predict.poisGPBayes}
\alias{predict.poisGPBayes}
\title{Predictive Quantiles or Return Levels for a EV Model
of class \code{"poisGPBayes"}.}
\usage{
\method{predict}{poisGPBayes}(object, ...)
}
\arguments{
\item{object}{An object with class \code{"poisGPBayes"}
representing a Poisson-GP fitted model with Bayesian inference
results.}

\item{...}{Other arguments passed to
\code{\link{predict.poisGPBayes0}}.}
}
\value{
An object with class \code{"predRL"} inheriting from
\code{"data.frame"}.
}
\description{
Prediction of a Bayesian Extreme-Value model of type Poisson-GP.
}
\details{
This method actually call \code{\link{predict.poisGPBayes0}}.
\code{object} is first coerced into the \code{"poisGPBayes0"}
simpler class for which the \code{predict} method is fully
implemented.
}
\seealso{
\code{\link{predict.poisGPBayes0}}.
}
