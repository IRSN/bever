## ****************************************************************************
##' Credible interval for a unimodal density, usually a posterior
##' distribution in a Bayesian framework.
##'
##' In some applications, \code{dist} will be a specific distribution
##' within a parametric family. This can occur for instance when a
##' conjugacy property holds. 
##' 
##' @title Find a Credible Interval for a Given Unimodal Probability
##' Distribution
##' 
##' @param dist Either a character describing a distribution such as
##' \code{"norm"} or \code{"gamma"}, or an object with class
##' \code{"density"} or a numeric vector assumed to provide a sample
##' of the distribution.
##'
##' @param level Numeric value specifying the probability level.  A
##' vector of length \eqn{> 1} is not accepted yet.
##' 
##' @param type When set to \code{"HPD"}, the credible interval is
##' Highest Posterior Density interval. When set to \code{"eqtail"}
##' the bounds are found by using the cumulative distribution function
##' and by allocating the same probability to the left and right
##' tails.
##' 
##' @param plot Logical. If \code{TRUE} a plot will be shown.
##' 
##' @param ... Other arguments to be passed to the distribution
##' functions: density, quantile.
##'
##' @return A numeric vector with the two bounds or an array when
##' \code{type} is \code{"both"}.
##'
##' @note When \code{dist} is intended to describe a distribution
##' within a parametric family, it must be taken as the character
##' string which is prefixed by \code{"d"}, \code{"p"} and \code{"q"}
##' to get the density, the distribution and the quantile
##' functions. The values of the parameters can passed through the
##' dots \dots See \bold{Examples}.
##' 
##' @examples
##' credInt("gumbel", level = 0.95, plot = TRUE)
##' credInt("gamma", level = 0.95, plot = TRUE, shape = 4)
##' credInt("gumbel", level = 0.95, plot = TRUE, type = "both")
##' x <- rgamma(n = 6000, shape = 4)
##' credInt(x, level = 0.95, plot = TRUE)
##' ## nearly the same
##' d <- density(x)
##' credInt(x, level = 0.95, plot = TRUE)
credInt <- function(dist, level = 0.95,
                    type = c("HPD", "eqtail", "both"),
                    plot = FALSE, ...) {

    if (length(level) > 1) {
        stop("'level' must be of length one")
    }
    
    check <- FALSE
    alpha <- 1 - level
    type <- match.arg(type)
    
    if (class(dist) == "density") {
        
        d <- function(x) {
            approx(x = dist$x, y = dist$y, xout = x)$y
        }
        
        py <- cumsum(dist$y)* diff(dist$x[1:2])
        
        q <- function(prob) {
            approx(x = py, y = dist$x, xout = prob)$y 
        }

        Ldots <- list()
        
    } else if (class(dist) == "character") {

        q <- match.fun(paste0("q", dist))
        d <- match.fun(paste0("d", dist))
        Ldots <- list(...)
        
    } else if (class(dist) == "numeric") {
        
        q <- function(prob) {
            quantile(x = dist, prob = prob) 
        }

        dx <- density(dist)
        d <- function(x) {
            approx(x = dx$x, y = dx$y, xout = x)$y
        }
        Ldots <- list(...)
        
    } else {
        stop("'x' must be of class \"density\", \"character\"",
             " or \"numeric\".")
    }

    resHPD <- reseqtail <- numeric(0)
    
    if (type != "eqtail") {
        g <- function(alphaL) {
            q(1 - alpha + alphaL, ...) - q(alphaL, ...)
        }
        
        opt <- optimize(f = g, interval = c(0, alpha), maximum = FALSE)
        
        resHPD <- c(q(opt$minimum, ...), q(1 - alpha + opt$minimum, ...))
    }

    if (type != "HPD") {
        reseqtail <- q(c(alpha / 2, 1.0 - alpha/2), ...)
    }

    res <- rbind("HPD" = resHPD, "eqtail" = reseqtail)
    colnames(res) <- c("L", "U")
    
    if (plot) {

        cols <- c("HPD" = "orangered", "eqtail" = "SpringGreen4")
        
        if (check) {
            cat("value of the density at bounds\n") 
            print(dcheck <- do.call(d, c(list(x = res), Ldots)))
            cat("\n")
        }
        x <- seq(from = q(1e-6, ...), to = q(1 - 1e-6, ...),
                 length.out = 300)
        dx <- do.call(d, c(list(x = x), Ldots))
        mdx <-  max(dx, na.rm = TRUE)

        plot(x, dx, type = "n", ylab = "", xlab = "",
             ylim = c(0, 1.2 * mdx), yaxt = "n")
        
        if (type != "both") {
            ind <- (x >= res[1] & x <= res[2])
            polygon(x = c(res[1], x[ind], res[2]),
                    y = c(0, dx[ind], 0),
                    border = NA, col = "lavender")
            yt <- 1.1 * mdx
            names(yt) <- names(type) <- type
        } else {
            type <- c("HPD" = "HPD", "eqtail" = "eqtail")
            yt <- c("HPD" = 1.15, "eqtail" = 1.05) * mdx
        }
        for (tt in type) {
            xt <- max(res[ , 2])
            segments(x0 = res[tt, 1], x1 = res[tt, 2],
                     y0 = yt[tt], y1 = yt[tt],
                     lwd = 3, col = cols[tt])
            text(x = xt, y =  yt[tt],
                 labels = type[tt], pos = 4,  col = cols[tt])
            abline(v = res[tt, ], col = cols[tt])
        }
        
        abline(h = 0)
        lines(x, dx, col = "SteelBlue3", lwd = 3)
        if (check) abline(h = dcheck[1], lty = "dotted")
    }
    
    res
}

## ***************************************************************************
##' Format levels of confidence or probability e.g. to build variable
##' names.
##'
##' @title Format Levels of Confidence or Probabilities
##' 
##' @param level Numeric vector of levels of confidence or
##' probabilities between \code{0} and \code{1}.
##'
##' @return A character vector.
##'
##' @examples
##' formatLevel(c(0.70, 0.95))
formatLevel <- function(level) {
    paste(gsub("\\.0$", "", sprintf("%4.1f", 100 * level)), "%", sep = "")
}

.parNames <- list(GEV = list(names = c("loc", "scale", "shape"),
                      greek = c("mu", "sigma", "xi")),
                  GPD = list(names = c("loc", "scale", "shape"),
                      greek = c("mu", "sigma", "xi")),
                  GPD2 = list(names = c("scale", "shape"),
                      greek = c("sigma", "xi"),
                      greek2 = c("sigma[u]", "xi")),
                  pois = list(names = "rate",
                      greek = "lambda"),
                  poisGP = list(names = c("lambda", "scale", "shape"),
                      greek = c("lambda", "sigma", "xi"),
                      names2 = c("rate", "scale", "shape")))
                  
