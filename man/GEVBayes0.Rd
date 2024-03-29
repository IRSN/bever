% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GEVBayes0.R
\name{GEVBayes0}
\alias{GEVBayes0}
\title{Create a Posterior for a GEV Model}
\usage{
GEVBayes0(MCMC, blockDuration = 1.0,
          MAP = NULL,
          yMax = NULL,
          nMax = length(yMax))
}
\arguments{
\item{MCMC}{An object that can be coerced into a matrix containing
the MCMC iterates. It should have the burnin period removed
and be thinned if necessary.}

\item{blockDuration}{The block duration given as a single positive
numeric value. The GEV distribution which parameters are
sampled in \code{MCMC} refers to the maximum on a period with
duration \code{blockDuration}.}

\item{MAP}{An optional vector of Maximum A Posteriori for the
parameter vector. Should be named with names matching the
colnames of \code{MCMC}.}

\item{yMax}{An optional vector of observations.}

\item{nMax}{An optional number of observations. Useful only when
\code{yMax} is not given.}
}
\value{
An object with class \code{"GEVBayes0"} inheriting from
    \code{"Bayes0"}. This object can be used to produce RL plots.
}
\description{
Create a "Poor Man's" Posterior for a GEV model using MCMC
iterates.
}
\note{
The argument \code{yMax} is intended for the classical
    framework where block maxima are used corresponding to a
    constant block duration. This is equivalent to using the
    \code{potData} argument with the value

\code{potData(MAX.data = as.list(yMax), MAX.effDuration =
rep(blockDuration, length(yMax))}.
}
\examples{
require(revdbayes)
## ========================================================================
## Portpirie data. Note that 'yMax' is only used for graphics later
## ========================================================================
prior <- set_prior(prior = "flatflat", model = "gev")
post <- rpost_rcpp(n = 10000, model = "gev", prior = prior,
                       data = portpirie)
## retrieve the MAP within the object
MAP <- post$f_mode
names(MAP) <- c("loc", "scale", "shape")

postGEV0 <- GEVBayes0(MCMC = post$sim_vals, yMax = portpirie, MAP = MAP)

## ========================================================================
## some methods
## ========================================================================
summary(postGEV0)
coef(postGEV0)
vcov(postGEV0)

## ========================================================================
## RL plot
## ========================================================================
RL0 <- RL(postGEV0)
autoplot(postGEV0) + ggtitle("GEV fit to Portpirie data")

## ========================================================================
## predictive distribution
## ========================================================================
pred <- predict(postGEV0)
autoplot(pred) 
autoplot(predict(postGEV0, newDuration = 100)) +
    ggtitle("Prediction for a 'new' period of 100 years")

}
\seealso{
\code{\link{RL}} method to generate a data frame of
    "classical" return levels (as shown on a classical RL plot),
    \code{\link{predict.GEVBayes0}} to generate a data frame of
    predictive return levels (as shown on a predictive RL plot).
}
