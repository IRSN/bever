## ****************************************************************************
##' Autoplot an object with call \code{"poisGPBayes"}.
##'
##' @title \code{autoplot} Method for the \code{"poisGPBayes"} S3
##' Class
##'
##' @param object An object with class \code{"poisGPBayes"} describing
##' inference results for a Poisson-GP model.
##' 
##' @param ... Further arguments passed to \code{RL}.
##' 
##' @return An object with class \code{"gg"} inheriting from
##' \code{"ggplot"} which is typically used through the \code{print}
##' method.
##' 
##' @examples
##' ## ========================================================================
##' ## Use the Garonne data from Renext
##' ## ========================================================================
##' fitBayes <- poisGPBayes(y = Garonne$OTdata$Flow,
##'                         threshold = 2500, duration = 64)
##' 
##' ## ========================================================================
##' ## S3 methods: Return Levels, ...
##' ## ========================================================================
##' RLBayes <- RL(fitBayes)
##' autoplot(fitBayes)            ## nearly the same plot
autoplot.poisGPBayes <- function(object, ... ) {
    
    ## avoid warning...
    Type <- Period <- Quantile <- Level <- NULL
    
    ## RL <- as.data.frame(RL(object))
    RL <- RL(object, ...)
    RL2 <- tidyr::gather(RL, key = Type, value = Quantile, -c(Period, Level))
    
    df <- data.frame(Period = exp(Hpoints(object$nOT)) / object$estimate[1],
                     z = object$threshold + sort(object$yOT))    
    
    gg <- ggplot()

    ## Credible intervals as a ribbon
    gg <- gg +
        geom_ribbon(
            data = RL,
            mapping = aes_string(x = "Period", ymin = "L", ymax = "U"),
            fill = "SteelBlue2", alpha = 0.4)
    
    ## Quantiles (Mean, Median, Mode) as lines with a legend
    gg <- gg +
        geom_line(
            data = dplyr::filter(RL2, Type %in% c("Mean", "Median", "Mode")),
            mapping = aes_string(x = "Period", y = "Quantile", group = "Type",
                colour = "Type"), size = 0.8)
    
    ## Add sample points (stored through specific SLOTS of 'object'
    gg <- gg +
        geom_point(data = dplyr::filter(df, Period > 1),
                   mapping = aes_string(x = "Period", y = "z"),
                   size = 1, alpha = 0.7)
    
    gg <- gg + scale_x_log10() 
  
    gg
    
    
}

## ****************************************************************************
##' Autoplot an object with class \code{"RL.poisGPBayes"} containing
##' Return Levels and related inference results.
##'
##' @title \code{autoplot} Method for the \code{"RL.poisGPBayes"} S3
##' Class
##'
##' @param object An object with class \code{"RL.poisGPBayes"}
##' describing Return Levels and related Bayesian inference results
##' for a Poisson-GP model.
##' 
##' @param ... Further arguments passed to \code{RL}.
##' 
##' @return An object with class \code{"gg"} inheriting from
##' \code{"ggplot"}, typically used through the \code{print} method.
##'  
##' @seealso \code{\link{autoplot.poisGPBayes}}.
##' 
autoplot.RL.poisGPBayes <- function(object, ... ) {

    ## avoid NOTEs on check
    Period <- Type <- Quantile <- Level <- NULL

    aL <- attributes(object)
    
    RL2 <- tidyr::gather(object, key = Type, value = Quantile, -c(Period, Level))
 
    df <- data.frame(Period = exp(Hpoints(aL$nOT)) / aL$estimate[1],
                     z = aL$threshold + sort(aL$yOT))    
    
    gg <- ggplot()
    
    ## Credible intervals as a ribbon
    gg <- gg +
        geom_ribbon(
            data = object,
            mapping = aes_string(x = "Period", ymin = "L", ymax = "U"),
            fill = "SteelBlue2", alpha = 0.4)
    
    ## Quantiles (Mean, Median, Mode) as lines with a legend
    gg <- gg +
        geom_line(
            data = dplyr::filter(RL2, Type %in% c("Mean", "Median", "Mode")),
            mapping = aes_string(x = "Period", y = "Quantile", group = "Type",
                colour = "Type"), size = 0.8)
    
    ## Add sample points (stored through specific ATTRIBUTES of 'object'
    gg <- gg +
        geom_point(data = dplyr::filter(df, Period > 1),
                   mapping = aes_string(x = "Period", y = "z"),
                   size = 1, alpha = 0.7)
    
    gg <- gg + scale_x_log10() 

    gg
    
}
