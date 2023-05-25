## ****************************************************************************
##' Autoplot an object with call \code{"poisGPBayes"}.
##'
##' @title \code{autoplot} Method for the \code{"poisGPBayes"} S3
##'     Class
##'
##' @method autoplot poisGPBayes0
##' 
##' @usage 
##' \method{autoplot}{poisGPBayes0}(object,
##'          which = c("RL", "pp"),
##'          level = 0.70,
##'          points = c("H", "p", "none"), a = 0.5,
##'          allPoints = FALSE,
##'          trace = 0, ...)
##' 
##' @param object An object with class \code{"poisGPBayes"} describing
##'     inference results for a Poisson-GP model.
##' 
##' @param which The type of plot wanted.
##' 
##' @param level Credible level(s).
##'
##' @param points,a,allPoints \code{\link[potomax]{autoplot.poisGP}}.
##'
##' @param trace Integer level of verbosity.
##' 
##' @param ... Further arguments passed to \code{RL}.
##' 
##' @return An object with class \code{"gg"} inheriting from
##'     \code{"ggplot"} which is typically used through the
##'     \code{print} method.
##'
##' @method autoplot poisGPBayes0
##' @export
##' 
##' @examples
##' ## ========================================================================
##' ## Use the Garonne data from Renext
##' ## ========================================================================
##' fitBayes <- poisGPBayes(data = Garonne$OTdata$Flow,
##'                         threshold = 2500,
##'                         effDuration = Garonne$OTinfo$effDuration)
##' 
##' ## ========================================================================
##' ## S3 methods: Return Levels, ...
##' ## ========================================================================
##' RLBayes <- RL(fitBayes)
##' autoplot(fitBayes)            ## nearly the same plot
autoplot.poisGPBayes0 <- function(object,
                                  which = c("RL", "pp"),
                                  level = 0.70,
                                  points = c("H", "p", "none"), a = 0.5,
                                  allPoints = FALSE,
                                  trace = 0,
                                  ...) {
    
    ## avoid warning...
    Type <- Period <- Quantile <- Level <- NULL

    which <- match.arg(which)
    points <- match.arg(points)
    
    if ((points != "none") && is.null(object$data)) {
        points <- "none"
    }   

    if (which == "RL") {
        
        ## RL <- as.data.frame(RL(object))
        RLs <- RL(object, level = level, trace = trace)

        gg <- autoplot(RLs, ...)
        if (!allPoints)  gg <- gg + ylim(object$threshold, NA)
        
        if (points != "none") {

            df <- RP(object$data, points = points, a = a)$data
            
            if (nlevels(df$block) < 6 ) {
                group <- "block"
            } else {
                group <- "source"
            }
            
            if (group == "block") {
                gg <- gg + geom_point(data = df,
                                      mapping = aes_string(x = "T", y = "x",
                                                           group = "block",
                                          ## linetype = NULL, ## fill = "block"
                                          ## colour = "block",
                                                           shape = "block"),
                                      ...) 
            } else {
                gg <- gg + geom_point(data = df,
                                      mapping = aes_string(x = "T", y = "x",
                                                           group = "source",
                                          ## linetype = NULL, ## fill = "source"
                                          ## colour = "source",
                                                           shape = "source"),
                                      ...)
            }
        }

        ## gg <- gg + guides(shape = FALSE)
        
        if (FALSE) {
            gg <- gg + scale_colour_manual(
                name = "Quantile",
                values = c("SteelBlue3", "SpringGreen4", "purple", "orangered", 
                    "black", "orangered", "SpringGreen1", "SteelBlue1", "magenta"))
        
            gg <- gg + scale_linetype_manual(
                name = "Quantile",
                values = c(c("dashed", "longdash", "twodash", "solid") ,
                    rep("blank", 5)))
        }
        
        gg <- gg + scale_x_log10() 
        
        gg
        
    } else {
        cat("Not implemented yet!")
        
    }

    gg
    
}

## ****************************************************************************
##' Autoplot an object with class \code{"RL.poisGPBayes"} containing
##' Return Levels and related inference results.
##'
##' @title \code{autoplot} Method for the \code{"RL.poisGPBayes"} S3
##'     Class
##'
##' @param object An object with class \code{"RL.poisGPBayes"}
##'     describing Return Levels and related Bayesian inference
##'     results for a Poisson-GP model.
##' 
##' @param ... Further arguments passed to \code{RL}.
##' 
##' @return An object with class \code{"gg"} inheriting from
##'     \code{"ggplot"}, typically used through the \code{print}
##'     method.
##'  
##' @seealso \code{\link{autoplot.poisGPBayes0}}.
##'
##' @method autoplot RL.poisGPBayes
##' @export
autoplot.RL.poisGPBayes <- function(object, ... ) {

    ## avoid NOTEs on check
    Period <- Type <- Quantile <- Level <- NULL

    aL <- attributes(object)
    
    RL2 <- tidyr::gather(object, key = Type, value = Quantile, -c(Period, Level))
    RL2 <- within(RL2, Level <- as.factor(paste("Cred.", Level)))

    leg <- FALSE
    
    if (leg) {
        gg <- ggplot(mapping = aes_string(x = "T", y = "x",
                         group = "Type",
                         shape = "Type", colour = "Type", fill = "Type"))
    } else {
        gg <- ggplot() 
    }
        
    ## ========================================================================
    ## Credible intervals as a ribbon and two lines for the limits
    ## ========================================================================
    
    if (!leg) {
        gg <- gg +
            geom_ribbon(
                data = object,
                mapping = aes_string(x = "Period", ymin = "L", ymax = "U"),
                fill = "SteelBlue2", alpha = 0.4)
        gg <- gg +
            geom_line(
                data = object,
                mapping = aes_string(x = "Period", y = "L",
                    group = "Level", linetype = "Level", colour = "Level",
                shape = NULL))
        gg <- gg +
            geom_line(
                data = object,
                mapping = aes_string(x = "Period", y = "U",
                    group = "Level", linetype = "Level", colour = "Level",
                    shape = NULL))   
    }
    
    ## ========================================================================
    ## Quantiles (Mean, Median, Mode) as lines with a legend
    ## ========================================================================
    gg <- gg +
        geom_line(
            data = dplyr::filter(RL2, Type %in% c("Mean", "Median", "Mode")),
            mapping = aes_string(x = "Period", y = "Quantile",
                group = "Type", linetype = "Type", colour = "Type",
                shape = NULL),
            size = 0.8)
    
    ## gg <- gg + scale_shape_manual(
    ##     name = "Quantile",
    ##     values = c(rep(NA, 4), c(21, 22, 24, 23, 25)))
    
    gg <- gg + scale_colour_manual(
        name = "Quantile",
        values = c("SteelBlue3", "SpringGreen4", "purple", "orangered", 
            "black", "orangered", "SpringGreen1", "SteelBlue1", "magenta"))
    
    gg <- gg + scale_linetype_manual(
        name = "Quantile",
        values = c(c("dashed", "longdash", "twodash", "solid") ,
            rep("blank", 5)))
    
    gg <- gg + scale_x_log10() 

    gg
    
}




