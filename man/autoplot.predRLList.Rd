% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predRL.R
\name{autoplot.predRLList}
\alias{autoplot.predRLList}
\title{\code{autoplot} Method for a List of Prediction Results}
\usage{
\method{autoplot}{predRLList}(object, groupName = "", ...)
}
\arguments{
\item{object}{A named list with all its elements of class
\code{"predRL"}.}

\item{groupName}{A string used to label the legend.}

\item{...}{Not used yet.}
}
\value{
An object of class \code{"gg"} inheriting from
\code{"ggplot"}.
}
\description{
Autoplot a list of prediction results, usually to compare them.
}
\details{
The list must be given the S3 class \code{"predRLList"} in order
to call the \code{autoplot} method. The method can be used to draw
several curves on the same plot or to use facets thanks to the
functions \code{facet_grid} or \code{facte_wrap}.
}
\examples{
prior <- set_prior(prior = "flatflat", model = "gev")
post <- rpost_rcpp(n = 10000, model = "gev", prior = prior,
                   data = portpirie)
pfitGEV0 <- GEVBayes0(MCMC = post$sim_vals, yMax = portpirie)
pL <- list("1 year" = predict(pfitGEV0),
           "10 year"= predict(pfitGEV0, newDuration = 10))
class(pL) <- "predRLList"
g <- autoplot(pL, groupName = "period") +
         ggtitle("Predictive RL plot : Bayesian GEV for portpirie") +
             theme_gray()
g
g + facet_grid(NewDuration ~ ., scales = "free_y")

}
