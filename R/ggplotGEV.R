## ****************************************************************************
##' Autoplot an object with call \code{"GEVBayes0"}.
##'
##' @title \code{autoplot} Method for the \code{"GEVBayes0"} S3 Class
##'
##' @param object An object with class \code{"GEVBayes0"} describing
##'     inference results for a GEV model.
##'
##' @param which Not used yet. Could in the future define the type of
##'     plot to build. For now only the RL plot is available.
##'
##' @param level The level of confidence to display.
##'
##' @param points The plotting positions to use to show the points if
##'     the object has data shiipped with it. The value \code{"H"}
##'     corresponds to the \emph{H-points} or Nelson's plotting
##'     positions and the value \code{"p"} correponds to the
##'     \emph{p-points} of the \code{\link[stats]{ppoints}}
##'     function. See \code{\link[potomax]{RP.potData}} in the package
##'     \strong{potomax}.
##'
##' @param a Value of the argument \code{a} of the
##'     \code{\link[stats]{ppoints}} function if \code{points} is
##'     \code{"p"}. It will be ignored otherwise.
##' 
##' @param xVar The variable to be mapped to the \eqn{x} axis. The
##'     value \code{"T"} corresponds to a period: a multiple of the
##'     block duration attached to the object. The value \code{"p"}
##'     corresponds to an exccedance probability related to the block
##'     duration. The relation between the two is simply \eqn{p = 1 /
##'     T}.
##'
##' @param trace Integer level of verbosity.
##' 
##' @param ... Further arguments passed to \code{RL}.
##' 
##' @return An object with class \code{"gg"} inheriting from
##'     \code{"ggplot"} which is typically used through the
##'     \code{print} method.
##'
##' @importFrom dplyr filter
##'
##' @method autoplot GEVBayes0
##' @export
##' 
autoplot.GEVBayes0 <- function(object,
                               which = "RL",
                               level = 0.70,
                               points = c("p", "none"),
                               a = 0.5,
                               xVar = c("T", "p"),
                               ## allPoints = FALSE,
                               trace = 0,
                               ...) {
    which <- match.arg(which)
    points <- match.arg(points)
    xVar <- match.arg(xVar)
    
    ## avoid warning...
    Type <- Period <- Quantile <- Level <- NULL
    
    ## RL <- as.data.frame(RL(object))
    RL <- RL(object, level = level, ...)
    RL2 <- tidyr::gather(RL, key = Type, value = Quantile, -c(Period, Level))

    if (points != "none") {
        if (is.null(object$yMax)) {
            warning("No 'yMax' in 'object'")
            points <- "none"
        } else {
            nMax <- length(object$yMax)
            p <- (1 - ppoints(nMax, a = a))
            df <- data.frame(T = object$blockDuration / p,
                             p = p,
                             x = sort(object$yMax),
                             source = "data")
                ## df <- data.frame(Period = exp(Hpoints(object$nOT)) / object$estimate[1],
            ##                  z = object$threshold + sort(object$yOT))           
        }
    }
    
    RL <- within(RL, Level <- as.factor(paste("Cred.", Level)))

    gg <- autoplot(RL, xVar = xVar) 
   
        
    ## Add sample points (stored through specific SLOTS of 'object'
    if (points != "none") {
        gg <- gg +
            geom_point(data = dplyr::filter(df, T > 1),
                       mapping = aes_string(x = xVar, y = "x"),
                       alpha = 0.7)
    }

    if (FALSE) {
    
        gg <- gg + scale_colour_manual(
            name = "Quantile",
        values = c("SteelBlue3", "SpringGreen4", "purple", "orangered", 
            "black", "orangered", "SpringGreen1", "SteelBlue1", "magenta"))
        
        gg <- gg + scale_linetype_manual(
            name = "Quantile",
            values = c(c("dashed", "longdash", "twodash", "solid") ,
                   rep("blank", 5)))
        if (xVar == "T") {
            gg <- gg + scale_x_continuous(trans = .gumbel_trans_m,
                                          breaks = .gumBreaks_m,
                                          minor_breaks = .gumBreaks_m) +
                xlab("period / block duration")
        } else {
            gg <- gg + scale_x_continuous(trans = .gumbel_trans_p,
                                          breaks = .gumBreaks_p,
                                      minor_breaks = .gumBreaks_p) +
                xlab("prob of exceedance")
        }
    }
    
    gg
    
    
}

## ****************************************************************************
##' Autoplot an object with class \code{"RL.GEVBayes"} containing
##' Return Levels and related inference results.
##'
##' @method autoplot RL.GEVBayes
##' 
##' @usage
##' 
##' \method{autoplot}{RL.GEVBayes}(object, 
##'         xVar = c("T", "p"),
##'         trace = 0,... ) 
##' 
##' @title \code{autoplot} Method for the \code{"RL.GEVBayes"} S3
##' Class
##'
##' @param object An object with class \code{"RL.GEVBayes"} describing
##'     Return Levels and related Bayesian inference results for a GEV
##'     model. This object is generated by \code{\link{RL.GEVBayes0}}.
##'
##' ## @param level Credible level.
##'
##' @param xVar The variable used to match the \eqn{x}-axis. The value
##'     \code{"T"} will lead to use the column named \code{"Period"}
##'     of \code{object} while the value \code{"p"} will lead to use
##'     the column named \code{"Prob"}.
##'
##' @param trace Integer level of verbosity.
##' 
##' @param ... Further arguments passed to \code{RL}.
##' 
##' @return An object with class \code{"gg"} inheriting from
##'     \code{"ggplot"}, typically used through the \code{print}
##'     method.
##'  
##' @seealso \code{\link{autoplot.GEVBayes0}}.
##'
##' @importFrom tidyr gather
##' @method autoplot RL.GEVBayes
##' @export
##' 
autoplot.RL.GEVBayes <- function(object,
                                 ## level = 0.70,
                                 xVar = c("T", "p"),
                                 trace = 0,... ) {
    
    ## avoid NOTEs on check
    Period <- Type <- Quantile <- Level <- Prob <- NULL

    xVar <- match.arg(xVar)
    xVar <- c("T" = "Period", "p" = "Prob")[xVar]  
    
    aL <- attributes(object)
    
    RL2 <- tidyr::gather(object, key = Type, value = Quantile,
                         -c(Period, Level, Prob))

    gg <- ggplot()
    
    ## Credible intervals as a ribbon
    gg <- gg +
        geom_ribbon(
            data = object,
            mapping = aes_string(x = xVar, ymin = "L", ymax = "U"),
            fill = "SteelBlue2", alpha = 0.4)
    
    gg <- gg +
        geom_line(
            data = object,
            mapping = aes_string(x = xVar, y = "L", group = "Level",
                linetype = "Level", colour = "Level"))
    
    gg <- gg +
        geom_line(
            data = object,
            mapping = aes_string(x = xVar, y = "U", group = "Level",
                linetype = "Level", colour = "Level"))
   
    gg <- gg +
      geom_line(
            data = dplyr::filter(RL2, Type %in% c("Mean", "Median", "Mode")),
            mapping = aes_string(x = xVar, y = "Quantile",
                group = "Type", linetype = "Type", colour = "Type",
                shape = NULL),
            size = 0.8)

    if (xVar == "Period") {
        gg <- gg + scale_x_continuous(trans = .gumbel_trans_m,
                                      breaks = .gumBreaks_m,
                                      minor_breaks = .gumBreaks_m)
    } else {
        gg <- gg + scale_x_continuous(trans = .gumbel_trans_p,
                                      breaks = .gumBreaks_p,
                                      minor_breaks = .gumBreaks_p)
    }
        
    gg <- gg + scale_colour_manual(
        name = "Quantile",
        values = c("SteelBlue3", "SpringGreen4", "purple", "orangered", 
            "black", "orangered", "SpringGreen1", "SteelBlue1", "magenta"))
    
    gg <- gg + scale_linetype_manual(
        name = "Quantile",
        values = c(c("dashed", "longdash", "twodash", "solid") ,
            rep("solid", 5)))
    
    gg <- gg + xlab("Period") + ylab("Quantile")
    
    gg
    
}

