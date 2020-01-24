## ****************************************************************************
##' Compute return levels along with credible bounds for a minimal
##' GEV model with Bayesian inference results.
##'
##' @title Return Levels and Credible Intervals for a Poisson-GP Model
##'
##' @param object An object with class \code{\link{poisGPBayes}} representing
##' the inference results for a Poisson-GP model.
##'
##' @param period A vector of periods for which the return levels will
##' be computed.
##'
##' @param level The credible level 
##'
##' @param smooth Logical. If \code{TRUE} the bounds of the credible
##' intervals are smoothed against the period.
##'
##' @param ... Not used yet.
##'
##' @return An object with class \code{"RL.GEVBayes"} inheriting from
##' \code{"data.frame"}. 
##'
##' 
RL.GEVBayes0 <- function(object,
                         period,
                         level = 0.70,
                         smooth = TRUE,
                         ...) {
    
    eps <-  1e-4
    fLevel <- formatLevel(level)
    
    nm <- c("loc", "scale", "shape")
    if ((ncol(object$MCMC) != 3) || (!identical(colnames(object$MCMC), nm))) {
        stop("'MCMC' must be a matrix with 3 columns and  ",
             " colamnes c(\"loc\", \"scale\", \"rate\")")
    }
    
    if (missing(period)) {
        period <- c(1.8, 2, 5, 10, 20, 50, 75, 100,
                    125, 150, 175, 200, 250, 300, 500, 700, 1000)
    }

    res <- array(NA, dim = c(length(period), 6),
                 dimnames = list(NULL, c("Period", "Mode", "Median", "Mean", "L", "U")))
    
    for (i in seq_along(period)) {

        m <- period[i] / object$blockDuration
        
        x <- NSGEV::qGEV(rep(1 / m, nrow(object$MCMC)),
                         loc = object$MCMC[ , "loc"],
                         scale = object$MCMC[ , "scale"],
                         shape = object$MCMC[ , "shape"],
                         lower.tail = FALSE)
        
        res[i, "Period"] <- period[i]
        res[i, "Mean"] <- mean(x)
        if (!is.null(object$MAP) && !any(is.na(object$MAP))) {
            res[i, "Mode"] <-
                NSGEV::qGEV(1 / m,
                            loc = object$MAP["loc"],
                            scale = object$MAP["scale"],
                            shape = object$MAP["shape"],
                            lower.tail = FALSE) 
            
        } 
        res[i, "Median"] <- median(x)
        LU <- credInt(x, level = level)
        res[i, c("L", "U")] <- LU
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

    for (nm in c("nMax", "yMax", "estimate", "blockDuration")) {
        attr(res, nm) <- object[[nm]]
    }
    
    class(res) <- c("RL.GEVBayes", "data.frame")
       
    res
       

}
