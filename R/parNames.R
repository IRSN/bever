
.parNames <- list(GEV = list(names = c("loc", "scale", "shape"),
                      greek = c("mu", "sigma", "xi")),
                  GPD = list(names = c("loc", "scale", "shape"),
                      greek = c("mu", "sigma", "xi")),
                  GP = list(names = c("loc", "scale", "shape"),
                      greek = c("mu", "sigma", "xi")),
                  PP = list(names = c("loc", "scale", "shape"),
                      greek = c("mu", "sigma", "xi")),
                  GPD2 = list(names = c("scale", "shape"),
                      greek = c("sigma", "xi"),
                      revdbayes = c("sigma[u]", "xi")),
                  pois = list(names = "rate",
                      greek = "lambda"),
                  poisGP = list(names = c("lambda", "scale", "shape"),
                      greek = c("lambda", "sigma", "xi"),
                      revdbayes = c("lambda", "sigma[u]", "xi"),
                      names2 = c("rate", "scale", "shape")))

## ****************************************************************************
##' Given a candidate vector of parameter names, check that these
##' match a known naming rule suitable for the given model.
##' 
##' @title Check a Candidate Vector or Parameter Names for an
##' Extreme-Value Model
##' 
##' @param parNames Candidate vector of parameter names.
##' 
##' @param model A possible Extreme-Value model. This is
##' case-insensitive.
##'
##' @param all Logical. If \code{TRUE}, the vector is assumed to
##' contain all the parameter names, possibly in a different order.
##'
##' @param trace Integer level of verbosity.
##'
##' @return A list with the following items
##' 
##' \item{rule}{
##'
##' The name of the identified rule for naming the parameters. For
##' instance, \code{"names"} corresponds to the rule used in this
##' package, \code{"greek"} corresponds using the names of the greek
##' letters associated to the parameters.
##'
##' }
##' \item{inParNames}{
##'
##' A copy of the vector \code{parNames} that was given on input.
##'
##' }
##' \item{parNames}{
##'
##' The vector of names corresponding to the rule \code{"names"} hence
##' the convention used in this package.
##' 
##' }
##' \item{indIn, indOut}{
##' 
##' Vectors of integer indices relating to the conventional vector of
##' names. \code{indIn} corresponds to the given names, possibly
##' translated into the conventional rule. \code{indOut} corresponds
##' to the omitted names, and can only have a positive length when
##' \code{all} is \code{FALSE} on input.
##' 
##' }
##'
##' @importFrom Matrix invPerm
##' @export
##' 
##' @examples
##' checkParNames(parNames = c("sigma", "mu", "xi"), model = "GEV")
##' checkParNames(parNames = c("sigma", "xi"), model = "GEV", all = FALSE)
##' checkParNames(parNames = c("lambda", "sigma", "xi"), model = "poisgp")
checkParNames <- function(parNames, model,
                          all = TRUE,
                          trace = 0) {

    model <- toupper(model)
    alloModels <- toupper(names(.parNames))
    i <-  match(model, alloModels)
    
    if ((length(model) != 1) || is.na(i)) {
        stop("'model' does not match any allowed models")
    }

    alloParNames <- .parNames[[i]]
    
    if (all) {

        test <- sapply(alloParNames,
                       function(item) setequal(x = item, y = parNames))
        
        if (!any(test)) {
            stop("'parNames' does not match a known naming rule ",
                 "for the parameters")
        }
        
        rule <- names(alloParNames)[test]
        indIn <- match(parNames, alloParNames[[rule]])
        indOut <- setdiff(1:length(alloParNames), indIn)
        ip <- Matrix::invPerm(indIn)
        names(ip) <- alloParNames[[rule]]
        
        list(rule = rule,
             inParNames = alloParNames[[rule]],
             parNames = alloParNames[["names"]],
             indIn = indIn,
             indOut = integer(0),
             pos = ip)
    } else {
        
        test <- sapply(alloParNames, function(item) all(parNames %in% item))
        if (!any(test)) {
            stop("'parNames' does not match a known naming rule ",
                 "for the parameters")
        }

        rule <- names(alloParNames)[test]
        indIn <- match(parNames, alloParNames[[rule]])
        indOut <- setdiff(1:length(alloParNames[[rule]]), indIn)
        ip <- Matrix::invPerm(c(indIn, indOut))
        names(ip) <- alloParNames[[rule]]
        
        list(rule = rule,
             inParNames = alloParNames[[rule]],
             parNames = alloParNames[["names"]],
             indIn = indIn,
             indOut = indOut,
             pos = ip)
    }
        
}
