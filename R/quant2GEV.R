## ****************************************************************************
##' Find the vector of the three GEV parameters corresponding to a
##' vector of three quantiles for three given distinct probabilities.
##'
##' @details Given three distinct probabilities, there is a one-to-one
##' correspondance between the vector of the three corresponding
##' quantiles and the vector of the three GEV parameters so the vector
##' of quantiles can be used to re-parameterise the GEV distribution.
##' The quantile parameterisation can be preferred to define
##' informative priors based on expert knowledge, see Coles and Tawn.
##' 
##' @title Find GEV Parameters from Given Quantiles 
##'
##' @param q A numeric vector with length \eqn{3} containing distinct
##' values for the GEV quantiles.
##'
##' @param p A numeric vector with length \eqn{3} containing distinct
##' values for the probabilities.
##'
##' @param lower.tail Logical. If \code{TRUE} the values in \code{p}
##' are for the probability of exceedance, else they are for the
##' probability of non-exceedance.
##'
##' @param plot Logical. If \code{TRUE} a simple plot will illustrate
##' the zero-finding used to find the GEV shape \eqn{\xi}.
##'
##' @param eps A small numeric number used to decide when the shape
##' \eqn{xi} is close enough to zero.
##'
##' @param trace Integer level of verbosity.
##' 
##' @return A named numeric vector containing the GEV parameters.
##'
##' @references
##'
##' Coles S. and Tawn J. (1996). A Bayesian Analysis of Extreme
##' Rainfall Data \emph{Appl. Statist.} 45 (4), pp. 463-478.
##' 
##' @examples
##'
##' co <- quant2GEV(p = c(0.1, 0.01, 0.001),
##'                 q = c(60, 80, 120), lower.tail = FALSE)
##'
##' ## check the result
##' nieve::qGEV(p = c(0.1, 0.01, 0.001),
##'             loc = co["loc"], scale = co["scale"], shape = co["shape"],
##'             lower.tail = FALSE)
##' 
quant2GEV <- function(q, p, lower.tail = TRUE,
                      plot = FALSE, eps = 1e-7,
                      trace = 0) {

    if ((length(q) != 3) || (length(p) != 3) ||
        any(duplicated(p)) || any(duplicated(q))) {
        stop("'q' and 'p' must be of length 3 with no duplicated ",
             "value")
    }
    
    if (!lower.tail) p <- 1 - p
    ind <- order(p)
    p <- p[ind]
    q <- q[ind]
    r <- -log(p)
    s <- log(r)
    Cst <- (q[1] - q[2]) / (q[1] - q[3])
    
    F <- function(xi) {
        if (abs(xi) > eps) {
            nxi <- -xi
            -Cst + (r[1]^nxi - r[2]^nxi) / (r[1]^nxi - r[3]^nxi) 
        } else {
            -Cst + (s[1] - s[2]) / (s[1] -s[3]) * (1 + (s[3] - s[2]) * xi / 2)
        }
    }

    res <- uniroot(F, c(-1, 1))
    if (trace) cat("check ", res$f.root, "\n") 
   
    if (plot) {
        xis <- seq(from = -1, to =  1, by = 0.002)
        plot(xis, sapply(xis, F), pch = 16, cex = 0.4, col = "orangered")
        abline(v = res$root, h = 0)
    }
    
    xi <- res$root
    nxi <- -xi
    scale <- xi * (q[1] - q[2]) / (r[1]^nxi - r[2]^nxi)
    
    c("loc" = q[1] - (scale / xi) * (r[1]^nxi - 1),
      "scale" = scale,
      "shape" = xi)
   
}
