## ****************************************************************************
##' Predictive Return Level plot using an object containing the Return
##' Levels.
##'
##' The argument \code{object} contains a table of computed predictive
##' Return Levels, in correspondence with an exceedance
##' probability. The plot shows the tail-quantile function for the
##' predictive distribution.  
##'
##' @title \code{autoplot} Method for Predictive Return Levels
##' 
##' @param object An object with class \code{"predRL"} corresponding
##' to predictive Return Levels.
##' 
##' @param ... Not used yet.
##'
##' @return An object with class \code{"gg"} inheriting from
##' \code{"ggplot"}. 
##'
##' @note In practice the predictive distribution for an Extreme-Value
##' model has usually an heavy tail. Therefore at least for small
##' probabilities of exceedance \eqn{p} the curve is convex (upward
##' concave).
##'
##' @seealso \code{\link{autoplot.predRLList}} to compare predictive
##' plots: for several predicted durations, several priors, several
##' datasets, ...
##' 
autoplot.predRL <- function(object, ...) {

    g <- ggplot(data = object)
    g <-  g + geom_line(mapping = aes_string(x = "Prob", y = "Quant"))
    g <- g + scale_x_continuous(trans = .gumbel_trans_p,
                                breaks = .gumBreaks_p,
                                minor_breaks = .gumBreaks_p) +
        xlab("prob of exceedance")
    g <- g + theme_gray()
    g
    
}


## ****************************************************************************
##' Autoplot a list of prediction results, usually to compare them.
##'
##' The list must be given the S3 class \code{"predRLList"} in order
##' to call the \code{autoplot} method. The method can be used to draw
##' several curves on the same plot or to use facets thanks to the
##' functions \code{facet_grid} or \code{facte_wrap}.
##' 
##' @title \code{autoplot} method for a list of prediction results.
##'
##' @param object A named list with all its elements of class
##' \code{"predRL"}.
##'
##' @param groupName A string used to label the legend.
##' 
##' @param ... Not used yet.
##'
##' @return An object of class \code{"gg"} inheriting from
##' \code{"ggplot"}.
##'
##' @examples
##' prior <- set_prior(prior = "flatflat", model = "gev")
##' post <- rpost_rcpp(n = 10000, model = "gev", prior = prior, data = portpirie)
##' pfitGEV0 <- GEVBayes0(MCMC = post$sim_vals, yMax = portpirie)
##' pL <- list("1 year" = predict(pfitGEV0),
##'            "10 year"= predict(pfitGEV0, newDuration = 10))
##' class(pL) <- "predRLList"
##' g <- autoplot(pL, groupName = "period") +
##'     ggtitle("Predictive RL plot : Bayesian GEV for portpirie") + theme_gray()
##' g
##' g + facet_grid(NewDuration ~ ., scales = "free_y")
##' 
autoplot.predRLList <- function(object, groupName = "", ...) {
    
    oNames <- names(object)
    for (i in seq_along(object)) {
        object[[i]] <- within(object[[i]], NewDuration <- as.numeric(NewDuration))
        if (i == 1) {
            df <- data.frame(Name = oNames[[i]], object[[i]], stringsAsFactors = FALSE)
        } else {
            df <-
                dplyr::bind_rows(df,
                                 data.frame(Name = oNames[[i]],
                                            object[[i]], stringsAsFactors = FALSE))
            
        }
    }
    df <- within(df, Name <- as.factor(Name))
    
    g <- ggplot(data = df)
    g <-  g + geom_line(mapping = aes_string(x = "Prob", y = "Quant",
                            group = "Name", colour = "Name"),
                        size = 0.8, alpha = 0.8)
    
    ## g <- g + facet_grid(NewDuration ~ ., scales = "free_y")
    g <- g + labs(colour = groupName)
    
    g <- g + scale_x_continuous(trans = .gumbel_trans_p,
                                breaks = .gumBreaks_p,
                                minor_breaks = .gumBreaks_p) +
        xlab("prob of exceedance")
    
    g <- g + theme_gray()
    g
    
   
}

