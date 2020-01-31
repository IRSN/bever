\name{bever-package}
\alias{bever-package}
\alias{bever}
\docType{package}
\title{
\packageTitle{bever}
}
\description{
\packageDescription{bever}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{bever}
\packageIndices{bever}

}
\author{
\packageAuthor{bever}

Maintainer: \packageMaintainer{bever}
}

%% \references{ }
%% \keyword{ package }
%% \seealso{ }
\examples{
## =========================================================================
## 'portirie' example.  Wrap the results given by 'revdbayes',
## retrieving the MAP within the object
## =========================================================================

prior <- set_prior(prior = "flatflat", model = "gev")
post <- rpost_rcpp(n = 10000, model = "gev", prior = prior,
                            data = portpirie)
MAP <- post$f_mode
names(MAP) <- c("loc", "scale", "shape")

postGEV0 <- GEVBayes0(MCMC = post$sim_vals, yMax = portpirie, MAP = MAP)

## ========================================================================
## RL plot
## ========================================================================

RL0 <- RL(postGEV0)
autoplot(postGEV0) + ggtitle("GEV fit to Portpirie data")

## ========================================================================
## predictive distribution: compute and autoplot for a given 'new'
## duration
## ========================================================================

pred <- predict(postGEV0)
autoplot(pred, xVar = "T")
autoplot(predict(postGEV0, newDuration = 100)) +
    ggtitle("Prediction for a 'new' period of 100 years")

## ========================================================================
## Compare the predictive distribution corresponding to a new duration
## identical to the (constant) block duration.
##
## CAUTION: Models with blockDuration != 1.0 have not been fully
## checked yet.
## ========================================================================

g <- autoplot(postGEV0, xVar = "p")
g <- g + autolayer(pred, aes = TRUE)
g 

}
