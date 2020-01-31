## ****************************************************************************
##' Compute return levels along with credible bounds for a Poisson-GP
##' model with Bayesian inference results.
##'
##' Three return levels are computed, named \code{"mode"},
##' \code{"median"} and \code{"mean"} in relation to the posterior
##' distribution of each return level. The column \code{"mode"} does
##' not contain the value of the marginal mode of the return levels,
##' but rather the modal curve obtained by plugging the MAP parameter
##' into the quantile function. This leads to a smooth curve; the
##' value for a given probability or period is (usually slightly)
##' different from the posterior marginal mode.
##'
##' @method RL poisGPBayes
##' 
##' @title Return Levels and Credible Intervals for a Poisson-GP Model
##'
##' @param object An object with class \code{\link{poisGPBayes}} representing
##' the inference results for a Poisson-GP model.
##'
##' @param period A vector of periods for which the return levels will
##' be computed.
##'
##' @param level The credible level. 
##'
##' @param smooth Logical. If \code{TRUE} the bounds of the credible
##' intervals are smoothed against the period.
##'
##' @param ... Not used yet.
##'
##' @return An object with class \code{"RL.poisGPBayes"} inheriting
##' from \code{"data.frame"}. 
##'
##' @seealso The \code{\link[potomax]{RL}} generic function.
##' 
##' @examples
##' fit <- poisGPBayes(data = Garonne$OTdata$Flow,
##'                    threshold = 2500,
##'                    effDuration = Garonne$OTinfo$effDuration)
##' ## Return Levels
##' RL(fit)
##' autoplot(fit)
RL.poisGPBayes <- function(object,
                           period = NULL,
                           level = 0.70,
                           smooth = TRUE,
                           ...) {
    
    if (length(level) > 1) {
        warning("'level' can only be of length 1 for now. ",
                "Only the first element will considered")
        level <- level[1]
    }
    
    eps <-  1e-4
    fLevel <- formatLevel(level)
    
    nm <- c("lambda", "scale", "shape")
    if ((ncol(object$MCMC) != 3) || (!identical(colnames(object$MCMC), nm))) {
        stop("'MCMC' must be a matrix with 3 columns and  ",
             " colamnes c(\"lambda\", \"scale\", \"rate\")")
    }
    
    if (is.null(period)) {
        period <- c(1.1, 1.5, 2, 5, 10, 20, 50, 75, 100,
                    125, 150, 175, 200, 250, 300, 500, 700, 1000)
    }

    res <- array(NA, dim = c(length(period), 6),
                 dimnames = list(NULL, c("Period", "Mode", "Median", "Mean", "L", "U")))
    
    for (i in seq_along(period)) {

        res[i, "Period"] <- period[i]

        prov <- period[i] * object$MCMC[ , "lambda"]
        x <- rep(NA, length(prov))
        ind <- (prov > 1.0)
        
        if (sum(ind) > nrow(object$MCMC) / 2) {
            x[ind]  <- object$threshold + potomax::qGPD2(1 / prov[ind],
                                                         scale = object$MCMC[ind, "scale"],
                                                         shape = object$MCMC[ind , "shape"],
                                                         lower.tail = FALSE)
            
            res[i, "Mean"] <- mean(x, na.rm = TRUE)
            res[i, "Median"] <- median(x, na.rm = TRUE)
            LU <- credInt(x[ind], level = level)
            res[i, c("L", "U")] <- LU
        }
            
        if (!is.null(object$MAP)) {
            res[i, "Mode"] <- object$threshold +
                potomax::qGPD2(1 / period[i] / object$MAP["lambda"],
                      scale = object$MAP["scale"], shape = object$MAP["shape"],
                      lower.tail = FALSE) 
            
        } 

    }
    
    if (smooth) {
        for (nm in c("L", "U")) {
            ss <- smooth.spline(x = log(res[ , "Period"]), y = res[ , nm], df = 5)
            res[ , nm] <- predict(ss)$y
        }
    }
    
    res <- as.data.frame(res)
    res <- data.frame(Period = res$Period, Level = fLevel,
                      Mode = res$Mode, Median = res$Median, Mean = res$Mean,
                      L = res$L, U = res$U)

    for (nm in c("nOT", "data", "estimate", "threshold")) {
        attr(res, nm) <- object[[nm]]
    }
    
    class(res) <- c("RL.poisGPBayes", "data.frame")
    
    res
}

## ****************************************************************************
##' Compute return levels along with credible bounds for a Poisson-GP
##' model with Bayesian inference results, "Poor Man's version".
##'
##' When in the MCMC object the rate of the Poisson process is not
##' sampled (hence there is no column \code{"lambda"} in the matrix of
##' MCMC iterates) a new suitable column is added by sampling from the
##' marginal posterior, and the \code{RL} method for the
##' \code{"poisGPBayes"} class is used.
##'  
##' @title Return Levels and Credible Intervals for a Poisson-GP Model 
##' "Poor Man's version".
##' 
##' @param object The \code{poisGPBayes0} object to use.
##'
##' @param period See  \code{\link{RL.poisGPBayes}}.
##'
##' @param level See  \code{\link{RL.poisGPBayes}}.
##'
##' @param smooth See \code{\link{RL.poisGPBayes}}.
##'
##' @param ... Other arguments to be passed to \code{\link{RL.poisGPBayes}}.
##'
##' @return An object with class \code{"RL.poisGPBayes"} inheriting
##' from \code{"data.frame"}.
##'
RL.poisGPBayes0 <- function(object,
                            period,
                            level = 0.70,
                            smooth = TRUE,
                            ...) {
    
    if (length(level) > 1) {
        warning("'level' can only be of lenght 1 for now. ",
                "Only the first element will considered")
        level <- level[1]
    }
    
    ## Patch: could this be done in a better way? 
    if (object$lambdaOut) {
        
        if (is.na(object$nOT) || length(object$nOT) == 0) {
            stop("Since 'lambda' is not sampled in the MCMC iterates ",
                 "'object' must contain a valid 'nOT' element")
        }
        if (is.na(object$effDuration) || length(object$effDuration) == 0) {
            stop("Since 'lambda' is not sampled in the MCMC iterates ",
                 "'object' must contain a valid 'effDuration' element")
        }
        
        object$MCMC <-  cbind(lambda = rgamma(nrow(object$MCMC),
                                  shape = 1 + object$nOT,
                                  rate = 0 + object$effDuration),
                              object$MCMC)
    }
    RL.poisGPBayes(object, level = level, ...)
    
}

