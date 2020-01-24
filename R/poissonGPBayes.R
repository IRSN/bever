## ****************************************************************************
##' Create an object representing a Poisson-GP model with Bayesian
##' inference result using the \strong{revdbayes} package.
##' 
##' This function mainly relies on the \pkg{revdbayes} package which
##' is used to produce the results of the Bayesian inference using
##' MCMC. However the Maximum A Posteriori (MAP) for the parameter
##' vector is computed by maximising the posterior function using a
##' numerical optimisation. This allows a better comparison with the
##' frequentist ML approach as used in \code{\link{poisGPML}}.
##' 
##' @title Poisson-GP Model with Bayesian Inference Results
##' 
##' @param y Numeric vector containing the marks.
##'
##' @param threshold The threshold.
##'
##' @param duration The duration of the observation period.
##'
##' @param a0,b0 Shape and rate for a gamma prior on the Poisson rate
##' \code{lambda}. Note that by choosing \code{a0 = 1} and \code{b0 =
##' 0} we define an improper flat prior on \eqn{(0, \infty)}{(0,
##' Inf)}. While \code{a0} is dimensionless and compares to \eqn{1 +} a
##' number of observations, \code{b0} has the dimension of a duration
##' and compares to an observation duration.
##'
##' @param priorGP A character defining which prior is used for the GP
##' part of the model as in the \bold{revdbayes} package.
##'
##' @param lowerGP,upperGP Bounds on the GP parameters.
##'
##' @param trace Level of verbosity.
##'
##' @return An object with S3 class \code{"poisGPBayes"}.
##'
##' @seealso \code{\link{poisGPML}} for a comparable objec with
##' frequentist inference results.
##'
##' @examples
##' ## ========================================================================
##' ## Use the Garonne data from Renext
##' ## ========================================================================
##' fit <- poisGPBayes(y = Garonne$OTdata$Flow,
##'                    threshold = 2500, duration = 64)
##' 
##' ## ========================================================================
##' ## Some S3 methods: RL for Return Levels, ...
##' ## ========================================================================
##' class(fit)
##' coef(fit)
##' RL(fit)
##' autoplot(fit)
##' 
poisGPBayes <- function(y,
                        threshold,
                        duration,
                        a0 = 1, b0 = 0,
                        priorGP = "flatflat",
                        lowerGP = c("scale" = 0, "shape" = -0.9),
                        upperGP = c("scale" = Inf, "shape" = Inf),
                        trace = 0) {
    
    yOT <- y[y > threshold] - threshold
    nOT <- length(yOT)
    an <- a0 + nOT
    bn <- b0 + duration

    priorObj <- set_prior(prior = priorGP, model = "gp",
                          min_xi = lowerGP["shape"])

    
    postObj <- rpost_rcpp(n = 10000, model = "gp",
                          prior = priorObj,
                          thresh = threshold,
                          data = y)
    
    ## manual transformation into a mcmc object
    MCMC <- mcmc(postObj$sim_vals[-(1:1000),  ])
    
    MCMC <- cbind(lambda = rgamma(nrow(MCMC), shape = an, rate = bn),
                  scale = MCMC[ , 1], shape = MCMC[ , 2])
    
    lower <- c("lambda" = 0, lowerGP)
    upper <- c("lambda" = Inf, upperGP) 
    
    logPriorGPFun <- match.fun(paste0("gp_", priorGP))

    rlFun <- function(par, period, deriv = 0) {
        qGPD(p = 1 / par[1] / period, loc = 0, scale = par[2], shape = par[3])
    }
    
    logPriorFun <- function(par) {
        if (b0 > 0) {
            logPrior <-  dgamma(par[1], shape = a0, rate = b0, log = TRUE)
        } else {
            logPrior <- 0.0
        }
        logPrior + logPriorGPFun(par[2:3],
                                 min_xi = lower[3], max_xi = upper[3])
        logPrior
    }
    
    negLogLikFun <- function(par) {
        -logLikFun(par)
    }
    
    logLikFun <- function(par) {
        logL <- 0.0
        logL <- logL + dpois(nOT, lambda = par[1] * duration, log = TRUE)
        logL <-  logL + sum(dGPD(yOT, loc = 0,
                                 scale = par[2], shape = par[3],
                                 log = TRUE))
        logL
    }

    logPostFun <- function(par) {
        logPriorFun(par) + logLikFun(par) 
    }

    negLogPostFun <- function(par) {
        -logPostFun(par)
    }

    ## fit <- fGPD(x = yOT)
    ## parSol <- c("lambda" = nOT / w, fit$estimate)
    parIni <- c("lambda" = nOT / duration, "scale" = mean(yOT), "shape" = 0.0)
 
    opts <- list("algorithm" = "NLOPT_LN_COBYLA",
                 "ftol_abs" = 1e-3,
                 "maxeval" = 1000, "print_level" = trace)
        
    opt <- nloptr::nloptr(x0 = parIni, eval_f = negLogPostFun,
                          lb = lower, ub = upper, opts = opts)

    if (opt$status %in% c(0, 1, 3, 4)) {
        MAP <- opt$solution
        names(MAP) <- c("lambda", "scale", "shape")
    } else {
        MAP <- NULL
    }
    
    obj <- list(duration = duration,
                y = y,
                yOT = yOT,
                nOT = nOT,
                threshold = threshold,
                a0 = a0,
                b0 = b0,
                rlFun = rlFun,
                logPriorFun = logPriorFun,
                logLikFun = logLikFun,
                negLogLikFun = negLogLikFun,
                logPostFun = logPostFun,
                negLogPostFun = negLogPostFun,
                estimate = MAP,
                MAP = MAP,
                MCMC = MCMC)

    class(obj) <- "poisGPBayes" ##  "poisGP"

    obj 
    
}

##' *******************************************************************************
##' Coercion method. 
##'
##' @title Coercion Method
##' 
##' @param object An object with class \code{"poisGPBayes0"} 
##'
##' @param ... Further arguments to be passed to \code{poisGPBayes0}.
##'
##' @return An object with class \code{"poisGPBayes0"}.
##'
##' @note Inasmuch the classes \code{"poisGPBayes"} and
##' \code{"poisGPBayes0"} contain named list with consistent naming
##' rules this method nearly only changes the class of the object. It
##' could be made to work in this way in the future.
##' 
as.poisGPBayes0.poisGPBayes <- function(object, ...) {

    poisGPBayes0(MCMC = object$MCMC,
                 threshold = object$threshold,
                 obsDuration = object$duration,
                 nOT = object$nOT,
                 yOT = object$yOT,
                 a0 = object$a0,
                 b0 = object$b0,
                 MAP = object$MAP,
                 ...)
    
}

## ****************************************************************************
##' Prediction of a Bayesian Extreme-Value model of type Poisson-GP.
##'
##' \code{object} is coerced into the \code{"poisGPBayes0"} simpler
##' class for which the \code{predict} method is fully implemented.
##'
##' @title Predictive Quantiles or Return Levels for a EV Model
##' of class \code{"poisGPBayes"}.
##'
##' @param object An object with class \code{"poisGPBayes"}
##' representing a Poisson-GP fitted model with Bayesian inference
##' results.
##'
##' @param ... Other arguments passed to
##' \code{\link{predict.poisGPBayes0}}.
##'
##' @return An object with class \code{"predRL"} inheriting from
##' \code{"data.frame"}.
##'
##' @seealso \code{\link{predict.poisGPBayes0}}.
##' 
predict.poisGPBayes <- function(object, ...) {
    predict(as.poisGPBayes0(object, ...))
}

coef.poisGPBayes <- function(object, ...) {
    coef(as.poisGPBayes0(object, ...))
}

vcov.poisGPBayes <- function(object, ...) {
    vcov(as.poisGPBayes0(object, ...))
}

summary.poisGPBayes <- function(object, ...) {
    summary(as.poisGPBayes0(object, ...))
}





