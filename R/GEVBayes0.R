# ****************************************************************************
##' Create a "Poor Man's" Posterior for a GEV model using MCMC
##' iterates.
##'
##' @title Create a Posterior for a GEV Model
##'
##' @usage
##'
##' GEVBayes0(MCMC, blockDuration = 1.0,
##'           MAP = NULL,
##'           yMax = NULL,
##'           nMax = length(yMax)) 
##' 
##' @param MCMC An object that can be coerced into a matrix containing
##'     the MCMC iterates. It should have the burnin period removed
##'     and be thinned if necessary.
##'
##' @param blockDuration The block duration given as a single positive
##'     numeric value. The GEV distribution which parameters are
##'     sampled in \code{MCMC} refers to the maximum on a period with
##'     duration \code{blockDuration}.
##' 
##' @param MAP An optional vector of Maximum A Posteriori for the
##'     parameter vector. Should be named with names matching the
##'     colnames of \code{MCMC}.
##'
##' @param yMax An optional vector of observations.
##'
##' @param nMax An optional number of observations. Useful only when
##'     \code{yMax} is not given.
##'
## @param potData an object of class \code{"potData"} describing the
## data that have been used to produce the MCMC iterates. This object
## can heterogeneous data with historical data. See
## \code{\link[potomax]{potData}}.
##' 
##' @return An object with class \code{"GEVBayes0"} inheriting from
##'     \code{"Bayes0"}. This object can be used to produce RL plots.
##'
##' @seealso \code{\link{RL}} method to generate a data frame of
##'     "classical" return levels (as shown on a classical RL plot),
##'     \code{\link{predict.GEVBayes0}} to generate a data frame of
##'     predictive return levels (as shown on a predictive RL plot).
##'
##' @note The argument \code{yMax} is intended for the classical
##'     framework where block maxima are used corresponding to a
##'     constant block duration. This is equivalent to using the
##'     \code{potData} argument with the value
##'
##' \code{potData(MAX.data = as.list(yMax), MAX.effDuration =
##' rep(blockDuration, length(yMax))}.
##'
##' @export
##' 
##' @examples
##' require(revdbayes)
##' ## ========================================================================
##' ## Portpirie data. Note that 'yMax' is only used for graphics later
##' ## ========================================================================
##' prior <- set_prior(prior = "flatflat", model = "gev")
##' post <- rpost_rcpp(n = 10000, model = "gev", prior = prior,
##'                        data = portpirie)
##' ## retrieve the MAP within the object
##' MAP <- post$f_mode
##' names(MAP) <- c("loc", "scale", "shape")
##'
##' postGEV0 <- GEVBayes0(MCMC = post$sim_vals, yMax = portpirie, MAP = MAP)
##'
##' ## ========================================================================
##' ## some methods
##' ## ========================================================================
##' summary(postGEV0)
##' coef(postGEV0)
##' vcov(postGEV0)
##' 
##' ## ========================================================================
##' ## RL plot
##' ## ========================================================================
##' RL0 <- RL(postGEV0)
##' autoplot(postGEV0) + ggtitle("GEV fit to Portpirie data")
##' 
##' ## ========================================================================
##' ## predictive distribution
##' ## ========================================================================
##' pred <- predict(postGEV0)
##' autoplot(pred) 
##' autoplot(predict(postGEV0, newDuration = 100)) +
##'     ggtitle("Prediction for a 'new' period of 100 years")
##'
GEVBayes0 <- function(MCMC,
                      blockDuration = 1.0,
                      MAP = NULL,
                      yMax = NULL,
                      nMax = length(yMax)) {

    if (!is.finite(blockDuration) || blockDuration < 0.0) {
        stop("'blockDuration' must be a finite postive number")
    }
    
    if (ncol(MCMC) != 3) {
        stop("'mcmc' does not have 3 columns")
    }
    
    cpn <-  checkParNames(parNames = colnames(MCMC), model = "GEV")
    if (length(cpn$indOut) != 0) {
        stop("'MCMC' colnames do not correspond to a GEV model")
    }
    nSim <- nrow(MCMC)
    
    colnames(MCMC) <- parNames <- cpn$parNames[cpn$indIn]

    nMax <- length(yMax)
    if (nMax == 0) {
        pd <- NULL
        yMax <- NULL
    } else {
        pd <- potomax::potData(MAX.data = as.list(yMax),
                               MAX.effDuration = rep(blockDuration, nMax))
    }

    meanPost <- apply(MCMC, 2, mean)
    sdPost <- apply(MCMC, 2, sd)
    medianPost <- apply(MCMC, 2, median)
    covPost <- cov(MCMC)

    if (missing(MAP)) {
        MAP <- rep(NA, length(parNames))
        names(MAP) <- parNames           
    } else {
        if (length(MAP) != length(parNames)) {
            stop("'MAP is expected to be of length ", length(parNames))
        }
    }           
    
    res <- list(MCMC = MCMC,
                blockDuration = blockDuration,
                model = "GEV",
                nMax = nMax,
                yMax = yMax,
                ## potData = pd,
                meanPost = meanPost,
                sdPost = sdPost,
                medianPost = medianPost,
                MAP = MAP,
                covPost = covPost)
    
    class(res) <- c("GEVBayes0", "Bayes0")
    res

}

## ****************************************************************************
##' Compute the predictive quantiles or return levels for a GEV
##' Model. The quantiles are those for the maximum on a "new" period
##' of time \code{newDuration} years.
##'
##' @method predict GEVBayes0
##'
##' @usage
##' \method{predict}{GEVBayes0}(object, newDuration = 1.0, prob,
##'         type = "RL",
##'         approx = FALSE,
##'         trace = 0, ...)
##' 
##' @title Predictive Quantiles or Return Levels for a GEV Model,
##' typically a Model for Block Maxima
##'
##' @param object a \code{GEVBayes0} object.
##' 
##' @param newDuration The duration of the 'new' period for which the
##' maximum is to be predicted. The newduration is expressed by using
##' the block duration in \code{object} as unit. So \code{newDuration
##' = 10} means a duration of \code{10 * object$blockDuration}.
##' 
##' @param prob A vector of exceedance probabilities. The default
##' value contains such as \eqn{0.01} and \code{0.001}.
##' 
##' @param type The type of prediction wanted. Remind that the
##' \code{predict} method of the \strong{revdbayes} package allows
##' several types of prediction.
##'
##' @param approx Logical. For the default \code{FALSE}, each value of
##' the tail quantile function is computed by zero-finding. For
##' \code{TRUE}, the quantiles are computed by using a fine of pairs
##' (argument, value) of the distribution function. This is likely to
##' be faster than \code{approx = FALSE} when \code{length(prob)} is
##' large.
##'
##' @param trace Integer level of verbosity.
##'
##' @param ... Not used yet.
##'
##' @return A data frame with three columns
##' \item{NewDuration}{
##'
##' A vector indicating the duration of the predicted period.
##'
##' }
##' \item{Prob}{
##'
##' The vector of exceedance probabilities.
##'
##' }
##' \item{Quant}{
##'
##' The vector of quantile or return levels.
##'
##' }
##'
##' The dataframe is given the S3 class \code{"predRL"} and it
##' receives several attributes such as the names of the factor
##' columns.
##'
##' @note The duration \code{newDuration} is understood as given in
##' years and not not in block duration. Thus it sould be kept
##' constant when comparing Block Maxima models with different block
##' durations, e.g. one and two years.
##'
##' @seealso \code{\link{GEVBayes0}}
##'
##' @method predict GEVBayes0
##' @export
##' 
predict.GEVBayes0 <- function(object,
                              newDuration = 1.0,
                              prob,
                              type = "RL",
                              approx = FALSE,
                              trace = 0,
                              ...) {
    
    type <-  match.arg(type)  
    
    if (missing(prob)) {
        prob <- c(0.800, 0.600, 0.500, 0.250, 0.200, 0.150, 0.100,
                  0.050, 0.030, 0.025, 0.020, 0.015,
                  0.010, 0.008, 0.005, 0.004, 0.003, 0.002, 0.001)
    } else {
        if (any(prob < 0.0) || any(prob > 1.0)) {
            stop("'prob' must be a numeric vector with values between ",
                 "0.0 and 1.0")
        }
        prob <- sort(prob, decreasing = TRUE)
    }
    nProb <- length(prob)
    T <- 1 / prob
  
    newDuration <- sort(newDuration)
    mStar <- newDuration / object$blockDuration

    ## =======================================================================
    ## Find a interval to localise the largest and the smallest
    ## quantiles
    ## =======================================================================
    
    probTestU <- c(0.900, 0.950, 0.975, 0.980, 0.990, 0.999)
    probTestL <- c(0.300, 0.200, 0.100, 0.050, 0.010, 0.001)
    
    for (iY in length(newDuration):1) {

        quant <- rep(NA, nProb)
        
        ## =====================================================================
        ## The posterior cumulative distribution function 
        ## =====================================================================
        FTilde <- function(q) {
            res <- rep(NA, length(q))
            for (i in seq_along(res)) {
                res[i] <- mean(nieve::pGEV(q = q[i],
                                           loc = object$MCMC[ , 1],
                                           scale = object$MCMC[ , 2],
                                           shape = object$MCMC[ , 3],
                                           lower.tail = TRUE)^mStar[iY])
            }
            res
        }
        if (trace) {
            cat("Finding quantiles 'qL' and 'qU' corresponding to the greatest",
                " and the smallest exceedance probability.\n")
        }
        
        OKL <- OKU <-  FALSE
        
        for (i in seq_along(probTestU)) {

            if (!OKU) {
                yU <- quantile(nieve::qGEV(p = min(prob) / mStar[iY],
                                           loc = object$MCMC[ , 1],
                                           scale = object$MCMC[ , 2],
                                           shape = object$MCMC[ , 3],
                                           lower.tail = FALSE),
                               prob = probTestU[i])
            }
            if (!OKL) { 
                yL <- quantile(nieve::qGEV(p = max(prob) / mStar[iY],
                                           loc = object$MCMC[ , 1],
                                           scale = object$MCMC[ , 2],
                                           shape = object$MCMC[ , 3],
                                           lower.tail = FALSE),
                               prob = probTestL[i])
            }
            
            if (!OKU) {
                qU <- try(uniroot(f = function(q) { FTilde(q) - 1 + min(prob) },
                                  interval = c(yL, yU)),
                          silent = TRUE)
                if (!inherits(qU, "try-error") && qU$estim.prec < 1e-4) {
                    OKU <-  TRUE
                    qU <- qU$root
                    if (trace) {
                        cat("At trial #", i, "'qU' is successfully",
                            "localised.\n") 
                    }
                }
                
                
            }
            
            if (!OKL) {
                qL <- try(uniroot(f = function(q) {FTilde(q) - 1 + max(prob)},
                                  interval = c(yL, yU)),
                          silent = TRUE)
                if (!inherits(qL, "try-error") && qL$estim.prec < 1e-4) {
                    OKL <-  TRUE
                    qL <- qL$root
                    if (trace) {
                        cat("At trial #", i, "'qL' is successfully",
                            "localised.\n") 
                    }
                }
            }
        }
        if (!OKL || !OKU) {
            stop("one of the minimal/maximal quantiles could not be localised")
        }
        quant[1L] <-  qL
        quant[nProb] <- qU
        
        if (nProb > 2) {
            
            ind <- (nProb - 1):2
            if (approx) {
                yGrid <- seq(from = qL, to = qU, length.out = 200)
                pGrid <- FTilde(yGrid)
                quant[ind] <- approx(x = pGrid, y = yGrid, xout = prob[ind])$y
            } else {
                
                for (i in ind) {
                    q <- try(uniroot(f = function(q){FTilde(q) - 1 + prob[i]},
                                     interval = c(qL, qU)),
                             silent = TRUE)
                    if (!inherits(q, "try-error") && q$estim.prec < 1e-2) {
                        OKL <-  TRUE
                        quant[i] <- q$root
                    } else {
                        stop("problem in quantile determinantion. Try",
                             " 'approx = TRUE'")
                    }
                }
            }
        }

        if (iY == length(newDuration)) {
            res <- data.frame(NewDuration = newDuration[iY], Prob = prob,
                              Quant = quant)
        } else {
            res <-
                dplyr::bind_rows(res,
                                 data.frame(NewDuration = newDuration[iY], Prob = prob,
                                            Quant = quant))
        }
    }
    res <- within(res, NewDuration <- factor(NewDuration))

    attr(res, "yMax") <- object$yMax
    attr(res, "newDuration") <- newDuration
    attr(res, "blockDuration") <- object$blockDuration
    attr(res, "model") <- "GEV"
    class(res) <-  c("predRL", "data.frame")
    res
    
}
