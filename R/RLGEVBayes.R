## ****************************************************************************
##' Compute return levels along with credible bounds for a minimal
##' GEV model with Bayesian inference results.
##'
##' @title Return Levels and Credible Intervals for a Poisson-GP Model
##'
##' @param object An object with class \code{\link{GEVBayes0}}
##' representing the inference results for a GEV (block maxima) model.
##'
##' @param period A vector of periods for which the return levels will
##' be computed.
##'
##' @param credintType The type of credible interval wanted. See
##' \code{\link{credInt}}.
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
##' \itemize{
##' \item{Period}{
##' The Return Period. This is expressed in the same unit as the block
##' duration that was given at the creation of \code{object} and which
##' is stored as \code{object$blockDuration}. So if
##' \code{blockDuration} is \code{2} (years) and \code{Period} is
##' \code{100} (years), the Return Level is for \eqn{50} blocks.
##' }
##' \item{Prob}{
##' The probability of exceedance of the Return Level for the
##' considered period. This is \eqn{T / w^\star}{T / wStar} where \eqn{T}
##' is the return period and \eqn{w^\star}{wStar} is the block duration.
##' }
##' \item{Level}{The credible level in formated form, e.g. \code{"95\%"} for
##' a provided level of \code{0.95}.}
##' \item{Mode, Median, Mean}{
##' While \code{Median} and \code{Mean} are the median and mean of the
##' return levels \eqn{\rho(T;\,\boldsymbol{\theta}^{[i]})}{\rho(T,
##' \theta[i, ])} corresponding to the MCMC iterates
##' \eqn{\boldsymbol{\theta}^{[i]}}{\theta[i, ]} of the GEV parameter
##' vector \eqn{\boldsymbol{\theta}}{\theta}, the mode is obtained
##' by plugging the MAP of the GEV parameter into the return level
##' \eqn{\rho(T;\,\boldsymbol{\theta})}{\rho(T, \theta)}.  The
##' corresponding Return Level curve can be called "modal".  If the
##' MAP is not available in \code{object}, the corresponding column
##' will contain \code{NA}.
##' }
##' \item{L, U}{
##' The Lower and Upper bounds of the credible interval.
##' }
##' }
##'
##' Note that when \eqn{m} is a small integer \eqn{>1} and \eqn{T = m w^\star}{T = m *
##' wStar}, the given probability is not the
##' probability that the maximum over \eqn{m} blocks with duration
##' \eqn{w^\star}{wStar} exceeds the given level. This only holds when
##' \eqn{m} is large. 
##' 
##' @references
##' Chapter 3 of
##'
##' Coles S. (2001) \emph{An Introduction to Statistical Modeling of Extreme 
##' Values}. Springer-Verlag.
##' 
RL.GEVBayes0 <- function(object,
                         period = NULL,
                         level = 0.70,
                         credintType = c("HPD", "eqtail"),
                         smooth = missing(period),
                         ...) {

    credintType <- match.arg(credintType)
    eps <-  1e-4
    fLevel <- formatLevel(level)
    
    nm <- c("loc", "scale", "shape")
    if ((ncol(object$MCMC) != 3) || (!identical(colnames(object$MCMC), nm))) {
        stop("'MCMC' must be a matrix with 3 columns and  ",
             " colamnes c(\"loc\", \"scale\", \"rate\")")
    }
    
    if (is.null(period)) {
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
        LU <- credInt(x, level = level, type = credintType)
        res[i, c("L", "U")] <- LU
    }
    
    if (smooth) {
        for (nm in c("L", "U")) {
            ss <- smooth.spline(x = log(res[ , "Period"]), y = res[ , nm], df = 5)
            res[ , nm] <- predict(ss)$y
        }
    }
    
    res <- as.data.frame(res)
    res <- data.frame(Period = res$Period,
                      Prob = object$blockDuration / res$Period,
                      Level = fLevel,
                      Mode = res$Mode, Median = res$Median, Mean = res$Mean,
                      L = res$L, U = res$U)

    for (nm in c("potData", "estimate", "blockDuration")) {
        attr(res, nm) <- object[[nm]]
    }
    
    class(res) <- c("RL.GEVBayes", "data.frame")
       
    res
       

}
