## ****************************************************************************
##' Compute Return Levels and confidence or credible intervals on these.
##'
##' 
##' @title Compute Return Levels using an Extreme-Value Model
##'
##' @param object An object representing an extreme-value model, i.e. 
##' the estimation and inference results.
##' 
##' @param ... Further arguments to be passed to method such as a
##' vector of period, the type of inferenc.
##' 
##' @return An object containing the return levels such as an object
##' inheriting from \code{data.frame}.
##' 
## RL <- function(object, ...) {
##     UseMethod("RL", object)
## }

## *****************************************************************************
##' This method is typically used with a \code{poisGPBayes} object to
##' access the method written for Poor Man's objects.
##'
##' @title Coerce an Object into a \code{poisGPBayes0} Object.
##'
##' @param object The object to be coerced. 
##' 
##' @param ... Not used yet.
##'
##' @return An object with class \code{"poisGPBayes0"} inheriting from
##' \code{"Bayes0"}.
##' 
as.poisGPBayes0 <- function(object, ...) {
    UseMethod("as.poisGPBayes0", object)
}
