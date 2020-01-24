## ****************************************************************************
##' Autoplot an object with call \code{"GEVBayes"}.
##'
##' @title \code{autoplot} Method for the \code{"GEVBayes"} S3
##' Class
##'
##' @param object An object with class \code{"GEVBayes"} describing
##' inference results for a GEV model.
##' 
##' @param ... Further arguments passed to \code{RL}.
##' 
##' @return An object with class \code{"gg"} inheriting from
##' \code{"ggplot"} which is typically used through the \code{print}
##' method.
##' 
autoplot.GEVBayes <- function(object, ... ) {
    
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
##' Autoplot an object with class \code{"RL.GEVBayes"} containing
##' Return Levels and related inference results.
##'
##' @title \code{autoplot} Method for the \code{"RL.GEVBayes"} S3
##' Class
##'
##' @param object An object with class \code{"RL.GEVBayes"}
##' describing Return Levels and related Bayesian inference results
##' for a GEV model.
##' 
##' @param ... Further arguments passed to \code{RL}.
##' 
##' @return An object with class \code{"gg"} inheriting from
##' \code{"ggplot"}, typically used through the \code{print} method.
##'  
##' @seealso \code{\link{autoplot.GEVBayes}}.
##' 
autoplot.RL.GEVBayes <- function(object, ... ) {

    ## avoid NOTEs on check
    Period <- Type <- Quantile <- Level <- NULL

    aL <- attributes(object)
    
    RL2 <- tidyr::gather(object, key = Type, value = Quantile, -c(Period, Level))

    if (!is.na(aL$nMax)) {    
        df <- data.frame(Period = exp(Hpoints(aL$nMax)),
                         z = sort(aL$yMax))    
    }
    
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
    if (!is.na(aL$nMax)) {    
        gg <- gg +
            geom_point(data = dplyr::filter(df, Period > 1),
                       mapping = aes_string(x = "Period", y = "z"),
                       size = 1, alpha = 0.7)
    }
    
    ## gg <- gg + scale_x_log10() 
    gg <- gg + scale_x_continuous(trans = .gumbel_trans_m,
                                  breaks = .gumBreaks_m,
                                  minor_breaks = .gumBreaks_m) 
    gg
    
}
