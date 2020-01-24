coef.Bayes0 <- function(object,
                           type = c("all", "mean", "median", "mode"),
                           ...) {

    type <-  match.arg(type)
    if (type == "mean") return(object$meanPost)
    else if (type == "median") return(object$medianPost)
    else if (type %in% c("mode", "MAP")) return(object$MAP)
    else if (type == "all") {
        return(rbind("mean"= object$meanPost,
                     "median" = object$medianPost,
                     "mode" = object$MAP))
    }
    
    
}

## ****************************************************************************
##' Covariance matrix for the estimated parameters of an \code{Bayes0}
##' object.
##'
##' The covariance is usually computed by taking the covariance of the
##' MCMC iterates as sored in \code{object}.
##' 
##' @title Covariance Matrix for a Fitted EV Model
##'
##' @param object An object whith class \code{"Bayes0"}.
##'
##' @param ... Not used yet.
##' 
##' @return A covariance matrix, usually of size 3.
##' 
vcov.Bayes0 <- function(object, ...) {
    object$covPost
}

summary.Bayes0 <- function(object, ...) {
    
    ans <- object
    
    class(ans) <-  "summary.Bayes0"
    
    ans
    
}

print.summary.Bayes0 <- function(x, ...) {
    cat(sprintf("%s Model Bayesian Inference\n", x$model))
    if (!is.null(x$blockDuration)) {
        cat(sprintf("o Block duration : %s\n",
                    format(x$blockDuration, digits = 2)))
    }
    ## automatically skipped when the slots do not exist
    cat(sprintf("o Number of OT observations: %d\n", x$nOT))
    if (!is.null(x$obsDuration)) {
        cat(sprintf("o Observation duration : %s\n",
                    format(x$obsDuration, digits = 2)))
    }
    cat(sprintf("o Number of blocks used: %d\n", x$nMax))    
    cat(sprintf("o Number of MCMC iterates: %d\n", nrow(x$MCMC)))
    ## cat(sprintf("   o MCMC iterates : %4.1f\n", x$blockDuration))
    co <- x$meanPost
    coText <- paste(format(co, digits = 2),
                    sprintf("[%s]", format(x$sdPost, digits = 2)))
    names(coText) <- names(co)
    cat("o Posterior mean [sd]:\n")
    print(noquote(coText))
    cat(sprintf("%s\n", x$comment))
}

print.Bayes0 <- function(x, ...) {
    print(summary(x))
}
