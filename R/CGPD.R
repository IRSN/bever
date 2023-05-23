## ***********************************************************************
##' Density, distribution function, quantile function and random
##' generation for the five-parameter Compound Generalized Pareto
##' Distibution (CGPD).
##'
##' @name CGPD
##' @rdname CGPD
##' @aliases dCGPD pCGPD qCGPD rCGPD
##' 
##' @usage
##' dCGPD(x, loc = 0.0, scale = 1.0, shape = 0.0,
##'       scaleN, shapeN, EN, IDN, log = FALSE)
##'
##' pCGPD(q, loc = 0.0, scale = 1.0, shape = 0.0,
##'       scaleN, shapeN, EN, IDN, lower.tail = TRUE)
##' 
##' qCGPD(p, loc = 0.0, scale = 1.0, shape = 0.0,
##'       scaleN, shapeN, EN, IDN, lower.tail = TRUE)
##' 
##' rCGPD(n, loc = 0.0, scale = 1.0, shape = 0.0,
##'       scaleN, shapeN, EN, IDN)
##' 
##' 
##' @title Density, Distribution Function, Quantile Function and
##' Random Generation for the Five-Parameter Compound Generalized
##' Pareto Distribution (CGPD)
##'
##' @param loc Location parameter. Numeric vector of length one.
##' 
##' @param scale Scale parameter. Numeric vector of length one.
##'
##' @param shape Shape parameter. Numeric vector of length one.
##'
##' @param scaleN Scale of the GPD for the \eqn{N} part. Along with
##' \code{shapeN} it provides the parameterisation for the
##' Binomial-Poisson-Negative Binomial familly.
##'
##' @param shapeN Shape of the GPD for the \eqn{N} part. Along with
##' \code{scaleN} it provides the parameterisation for the
##' Binomial-Poisson-Negative Binomial familly.
##'
##' @param EN Expectation of \eqn{N}. Along with \code{IDN} it
##' provides an alternative parameterisation for the \eqn{N} part.
##'
##' @param IDN Index of Dispersion of \eqn{N}. Along with \code{EN} it
##' provides an alternative parameterisation for the \eqn{N} part.
##' 
##' @param log Logical; if \code{TRUE}, densities \code{p} are
##' returned as \code{log(p)}.
##'
##' @param x,q Vector of quantiles.
##'
##' @param p Vector of probabilities.
##'
##' @param n Sample size.
##' 
##' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
##' P[X <= x], otherwise, P[X > x].
##'
##' @return A numeric vector with length equal to the length of the
##' first argument or of the parameters. 
##'
##' @details This distribution is that of the maximum \eqn{M} of
##' \eqn{N} i.i.d. r.vs \eqn{X_i} with distribution
##' \eqn{\textrm{GPD}(\mu,\,\sigma,\,\xi)}{GPD(\mu, \sigma, \xi)}
##' where \eqn{N} is a r.v. with non-negative integer values,
##' independent of the sequence \eqn{X_i}, and having a
##' \emph{Binomial}, \emph{Poisson} or \emph{Negative Binomial}
##' distribution. The distribution of \eqn{N} can be parameterized by
##' using two parameters \eqn{\mu_N}{\mu_N} and
##' \eqn{\sigma_N}{\sigma_N} in a GPD style, or alternatively by using
##' the two parameters \eqn{\mathrm{E}(N)}{EN} and
##' \eqn{\textrm{ID}(N)}{IDN} representing the expectation and the
##' index of dispersion of \eqn{N}. The three cases \emph{Binomial},
##' \emph{Poisson} and \emph{Negative Binomial} correspond to
##' \eqn{\textrm{ID}_N < 0}{IDN < 0}, \eqn{\textrm{ID}_N = 0}{IDN = 0}
##' and \eqn{\textrm{ID}_N > 0}{IDN > 0}.
##'
##' @section Caution: This distribution is of \emph{mixed-type}. It
##' has a probability mass at \eqn{-\infty}{-Inf} corresponding to the
##' possibility that \eqn{N=0} in which case \eqn{M} is the maximum of
##' an empty set, taken as \eqn{-\infty}{-Inf} corresponding to
##' \code{max(mumeric(0))}. Consequently a sample drawn by using
##' \code{rCGPD} contains \code{-Inf} values with positive
##' probability.
##' 
##' @author Yves Deville
##'
##' @references
##'  
##' Yves Deville (2019)  "Bayesian Return Levels in Extreme-Value Analysis"
##' \emph{IRSN technical report}.
##' 
##' @examples
##' set.seed(1)
##' ExpN <- runif(1)
##' IDN <- rexp(1, rate = 1)
##' scaleN <- 1 / ExpN
##' shapeN <- (IDN - 1) / ExpN
##' loc <- rnorm(1, mean = 0, sd = 10); scale <-  rexp(1)
##' shape <- rnorm(1 , mean = 0, sd = 0.1)
##' mass <- pCGPD(-Inf, scaleN = scaleN, shapeN = shapeN,
##'            loc = loc, scale = scale, shape = shape) 
##' q <- qCGPD(p = c(mass + 0.001, 0.999), scaleN = scaleN, shapeN = shapeN,
##'            loc = loc, scale = scale, shape = shape)
##' x <- seq(from = q[1] - 1, to = q[2], length.out = 200)
##' F <- pCGPD(x, scaleN = scaleN, shapeN = shapeN,
##'            loc = loc, scale = scale, shape = shape)
##' plot(x, F, type = "l", xlab = "", ylab = "", ylim = c(0, 1),
##'      col = "orangered")
##' abline(h = mass, col = "red")
##' f <- dCGPD(x, scaleN = scaleN, shapeN = shapeN,
##'            loc = loc, scale = scale, shape = shape)
##' plot(x, f, type = "l", col = "SteelBlue3", xlab = "", ylab = "")
##' title(main = sprintf(paste("ExpN = %4.1f IDN = %4.2f,",
##'                            "loc = %4.1f, scale = %4.2f, shape = %4.2f"),
##'                      ExpN, IDN, loc, scale, shape))
##' 
dCGPD <- function(x, loc = 0.0, scale = 1.0, shape = 0.0,
                  scaleN, shapeN,
                  EN, IDN, log = FALSE) {

    if (!missing(EN)) {
        if (!missing(scaleN) || !missing(shapeN)) {
            stop("when 'EN' is given, none of the two parameters ",
                 "'scaleN' and 'shapeN' should be given")
        }
        if (missing(shapeN)) {
            stop("when 'EN' is given 'IDN' must also be given")
        }
        scaleN <- 1.0 / EN
        shapeN <- (IDN - 1.0) / EN
        
    }

    ## allow each arg to be a vector with length > 1
    n <- max(c(length(x),
               length(loc), length(scale), length(shape),
               length(scaleN), length(shapeN)))
    x <- rep(x, length.out = n)
    loc <- rep(loc, length.out = n)
    scale <- rep(scale, length.out = n)
    shape <- rep(shape, length.out = n)
    scaleN <- rep(scaleN, length.out = n)
    shapeN <- rep(shapeN, length.out = n)
    
    res <- rep(NA, length(x))
    x1 <- x - loc
    ind <- (x1 < 0.0)
    if (any(ind)) {
        if (log) res[ind] <- -Inf
        else res[ind] <- 0.0
    }
    ## this can be improved!!!
    if (any(!ind)) {
        SY <- nieve::pGPD2(x1[!ind], scale = scale[!ind], shape = shape[!ind], lower.tail = FALSE)
        fY <- nieve::dGPD2(x1[!ind], scale = scale[!ind], shape = shape[!ind])
        res[!ind] <- nieve::dGPD2(SY, scale = scaleN[!ind], shape = shapeN[!ind]) * fY
        if (log) res[!ind] <- log(res[!ind])
    }
    res
}

pCGPD <- function(q, loc = 0.0, scale = 1.0, shape = 0.0,
                  scaleN, shapeN,
                  EN, IDN, lower.tail = TRUE) {
    
    if (!missing(EN)) {
        if (!missing(scaleN) || !missing(shapeN)) {
            stop("when 'EN' is given, none of the two parameters ",
                 "'scaleN' and 'shapeN' should be given")
        }
        if (missing(shapeN)) {
            stop("when 'EN' is given 'IDN' must also be given")
        }
        scaleN <- 1.0 / EN
        shapeN <- (IDN - 1.0) / EN   
    }

    ## allow each arg to be a vector with length > 1
    n <- max(c(length(q),
               length(loc), length(scale), length(shape),
               length(scaleN), length(shapeN)))
    q <- rep(q, length.out = n)
    loc <- rep(loc, length.out = n)
    scale <- rep(scale, length.out = n)
    shape <- rep(shape, length.out = n)
    scaleN <- rep(scaleN, length.out = n)
    shapeN <- rep(shapeN, length.out = n)
    
    mass <- nieve::pGPD2(1.0, scale = scaleN, shape = shapeN, lower.tail = FALSE)
    res <- rep(NA, length(q))
    q1 <-  q - loc
    ind <- (q1 < 0.0)

    if (any(ind)) res[ind] <- mass[ind]
    
    if (any(!ind)) {
        S <- nieve::pGPD2(q1[!ind], scale = scale[!ind], shape = shape[!ind],
                         lower.tail = FALSE)
        res[!ind] <- nieve::pGPD2(S, scale = scaleN[!ind], shape = shapeN[!ind],
                                 lower.tail = FALSE)
    }
    
    if (!lower.tail) res <- 1.0 - res
    res
    
}

qCGPD <- function(p, loc = 0.0, scale = 1.0, shape = 0.0,
                  scaleN, shapeN,
                  EN, IDN, lower.tail = TRUE) {

    if (!missing(EN)) {
        if (!missing(scaleN) || !missing(shapeN)) {
            stop("when 'EN' is given, none of the two parameters ",
                 "'scaleN' and 'shapeN' should be given")
        }
        if (missing(shapeN)) {
            stop("when 'EN' is given 'IDN' must also be given")
        }
        scaleN <- 1.0 / EN
        shapeN <- (IDN - 1.0) / EN
    }

    ## allow each arg to be a vector with length > 1
    n <- max(c(length(p),
               length(loc), length(scale), length(shape),
               length(scaleN), length(shapeN)))
    p <- rep(p, length.out = n)
    loc <- rep(loc, length.out = n)
    scale <- rep(scale, length.out = n)
    shape <- rep(shape, length.out = n)
    scaleN <- rep(scaleN, length.out = n)
    shapeN <- rep(shapeN, length.out = n)
    
    mass <- nieve::pGPD2(1.0, scale = scaleN, shape = shapeN, lower.tail = FALSE)
    if (!lower.tail) p <- 1.0 - p
    res <- rep(NA, length(p))
    ind <- (p <= mass)
    if (any(ind)) res[ind] <- -Inf
    if (any(!ind)) {
        q1 <- nieve::qGPD2(p[!ind], scale = scaleN[!ind], shape = shapeN[!ind],
                    lower.tail = FALSE)
        res[!ind] <- loc[!ind] + nieve::qGPD2(q1, scale = scale[!ind],
                                              shape = shape[!ind],
                                              lower.tail = FALSE)
    }
    res

}

rCGPD <- function(n, loc = 0.0, scale = 1.0, shape = 0.0,
                  scaleN, shapeN,
                  EN, IDN) {

    if (length(n) > 1 ||
        length(loc) > 1 || length(scale) > 1 || length(shape) > 1 ||
        length(scaleN) > 1 || length(shapeN) > 1 ||
        length(EN) > 1 || length(IDN) > 1) {
        stop("no arg of 'rGCPD' can have a length > 1")
    }
    
    if (!missing(EN)) {
        if (!missing(scaleN) || !missing(shapeN)) {
            stop("when 'EN' is given, none of the two parameters ",
                 "'scaleN' and 'shapeN' should be given")
        }
        if (missing(shapeN)) {
            stop("when 'EN' is given 'IDN' must also be given")
        }
        scaleN <- 1.0 / EN
        shapeN <- (IDN - 1.0) / EN
        
    }
    p <- runif(n)
    qCGPD(p, scaleN = scaleN, shapeN = shapeN,
          loc = loc, scale = scale, shape = shape)
}
