# ****************************************************************************
##' Create an object representing a "Poor Man's" Posterior for a
##' Poisson-GP model, mainly using MCMC iterates.
##'
##' @details The Poisson-GP model involves three parameters: one is
##' the rate \code{lambda} of the Poisson process describing the
##' exceedances over the threshold. The two other parameters are the
##' scale and the shape of the GP Distribution used for the excesses
##' over the threshold. The object MCMC can have either three columns
##' or two columns only. In the second case, the column corresponding
##' to the rate \code{lambda} is omitted; it is then assumed that
##' \code{lambda} is a posteriori independent of the GP parameters and
##' that it has a gamma posterior margin. The shape of the posterior
##' distribution is \eqn{a_n := a_0 + n_{\textrm{OT}}}{an = a0 + nOT}
##' and the rate is \eqn{b_n := b_0 + w_{\textrm{OT}}}{an = a0 + wOT},
##' where \eqn{w_{\textrm{OT}}}{wOT} and \eqn{w_{\textrm{OT}}}{wOT}
##' are the number of exceedances and the duration of the observation
##' period.
##'
##' @usage
##' poisGPBayes0(MCMC, threshold,
##'              data = NULL, effDuration,
##'              MAX.data = NULL, MAX.effDuration = NULL,
##'              OTS.data = NULL, OTS.threshold = NULL, OTS.effDuration = NULL,
##'              MAP = NULL,
##'              nOT = NA,
##'              a0 = 1.0, b0 = 0.0)
##' 
##' @title Create a Posterior for a Poisson-GP Model
##'
##' @param MCMC An object that can be coerced into a matrix containing
##' the MCMC iterates. It should have the burnin period removed and be
##' thinned if necessary.
##' 
##' @param data An optional structure containing the observations Over
##' the Threshold. The user should make sure that these data are
##' identical to those used to obtain \code{MCMC}. It can be a numeric
##' vector or an object inheriting from one of the two classes
##' \code{"Rendata"} (from \strong{Renext}) and \code{"potData"} (from
##' \strong{potomax}) in which case it will be coerced into the class
##' \code{"potData"}.
##'
##' @param threshold the POT threshold.
##' 
##' @param effDuration The observation effective duration. When
##' \code{data} can be coerced into \code{potData}, the provided
##' value for \code{effDuration} will be ignored and the value
##' of \code{object$OT$effDuration} will be used instead.
##'
##' @param nOT An optional integer giving the number of exceedances
##' during the observation period. This number is needed when
##' \code{MCMC} does not have a column for the rate of the Poisson
##' process. It should be consistent with \code{MCMC}.  When
##' \code{data} can be coerced into \code{potData}, the provided value
##' for \code{effDuration} will be ignored and the value of
##' \code{length(object$OT$data)} will be used instead. However in
##' this case no \code{MAX} or \code{OTS} data should exist in
##' \code{data}.
##'
##' @param MAX.data,MAX.effDuration Optional censored data (type block
##' maxima or \eqn{r}-largest) as used in
##' \code{\link[potomax]{poisGP}}. Ignored if \code{data} inherits
##' from one of the two classes \code{"Rendata"} and \code{"potData"}
##'
##' @param OTS.data,OTS.threshold,OTS.effDuration Optional censored
##' data as used in \code{\link[potomax]{poisGP}}. Ignored if
##' \code{data} inherits from one of the two classes \code{"Rendata"}
##' and \code{"potData"}
##'
##' 
##' @param MAP An optional vector of Maximum A Posteriori for the
##' parameter vector. Should be named with names matching the
##' colnames of \code{MCMC}.
##'
##' @param nOT An optional integer giving the number of exceedances
##' during the observation period. This number is needed when
##' \code{MCMC} does not have a column for the rate of the Poisson
##' process. It should be consistent with \code{MCMC}.
##'
##' @param a0,b0 The shape and the rate of the prior for \code{lambda}
##' (used only when \code{MCMC} has two columns).
##' 
##' @return An object with class \code{"GEVBayes0"} inheriting from
##' \code{"Bayes0"}. This object can be used to produce RL plots.
##'
##' @seealso \code{\link{RL}} method to generate a data frame of
##' "classical" return levels (as shown on a classical RL plot),
##' \code{\link{predict.poisGPBayes0}} to generate a data frame of
##' predictive return levels (as shown on a predictive RL plot).
##'
##' @section Caution: The user must be sure that the data provided by
##' using the dedicated arguments are consistent with the MCMC
##' iterates provided in \code{MCMC}. The 'data' dedicated arguments
##' are: \code{data}, \code{effDuration} or/and \code{MAX.data},
##' \code{MAX.effDuration} \code{OTS.data}, \code{OTS.threshold},
##' \code{OTS.effDuration}. Since attaching data to the object is
##' essentially for graphical purpose, it can be simpler not to attach
##' the data and to create a separated object with class
##' \code{"potData"}. This object can be used with its
##' \code{autolayer} method.
##' 
##' @examples
##' require(revdbayes)
##' ## ========================================================================
##' ## use the 'rainfall data': daily rainfalls over 57 years
##' ## ========================================================================
##' data(rainfall)
##' rainfall2 <- rainfall[!is.na(rainfall)]
##' u <- 40  ## threshold (mm)
##' w <- 57  ## obs duration (year)
##' nOT <- sum(rainfall2 > u)
##' prior <- set_prior(prior = "flatflat", model = "gp")
##' nSim <- 10000
##'
##' ## ========================================================================
##' ## get a posterior for the "excess part" of the model
##' ## ========================================================================
##' postGP <- rpost(n = 10000, model = "gp", prior = prior,
##'                 data = rainfall2, thresh = u)
##' MCMCGP <- postGP$sim_vals
##' MAP <- postGP$f_mode
##' names(MAP) <- c("scale", "shape")
##' 
##' ## ========================================================================
##' ## assuming a "flat gamma prior" for 'lambda', add MCMC iterates
##' ## for 'lambda' to those existing for GP (new matrix column)
##' ## ========================================================================
##' MCMCpoisGP <- cbind(lambda = rgamma(nSim, shape = 1 + nOT, rate = 0 + w),
##'                     MCMCGP)
##'
##' ## ========================================================================
##' ## create the object
##' ## ========================================================================
##' post0 <- poisGPBayes0(MCMC = MCMCpoisGP, threshold = u, effDuration = w) 
##' summary(post0)
##' coef(post0)
##'
##' ## ========================================================================
##' ## 'lambda' and the GP parameters are a posteriori independent
##' ## ========================================================================
##' cov2cor(vcov(post0))
##'
##' ## ========================================================================
##' ## create an object assuming 'lambda' independent of the GP param
##' ## ========================================================================
##' post1 <- poisGPBayes0(MCMC = MCMCGP, threshold = u,
##'                       effDuration = w, nOT = nOT, MAP = MAP)
##' 
poisGPBayes0 <- function(MCMC, threshold,
                         data = NULL, effDuration,
                         MAX.data = NULL, MAX.effDuration = NULL,
                         OTS.data = NULL, OTS.threshold = NULL,
                         OTS.effDuration = NULL, 
                         MAP = NULL,
                         nOT = NA,
                         a0 = 1.0, b0 = 0.0) {

    cd <- class(data)
    if (inherits(data, "Rendata") || inherits(data, "potData")) {

        pd <- as.potData(data)
        
        ## if (!missing(effDuration)) {
        ##     if (effDuration != pd$OT$effDuration) {
        ##         stop("Two different effective durations are provided in ",
        ##              "'effDuration' and in 'data'")  
        ##     }
        ## }
        if (!missing(effDuration) || !missing(nOT) || !missing(MAX.data) ||
            !missing(MAX.effDuration) || !missing(OTS.data) || 
            !missing(OTS.threshold) || !missing(OTS.effDuration)) {
            warning("Since 'data' has class \"", cd, "\" the formal arguments ",
                    "'effDuration', 'nOT', 'MAX.*' and 'OTS.*' are ignored")
        }
        
        effDuration <- pd$OT$effDuration
        nOT <- length(pd$OT$data)
        
    } else {

        if (!is.finite(effDuration) || effDuration < 0.0) {
            stop("'effDuration' must be a finite postive number")
        }
        
        if (!is.null(data)) {

            ## 'data' should be numeric here
            
            pd <- try(potomax::potData(data = data,
                                       effDuration = effDuration,
                                       MAX.data = MAX.data,
                                       MAX.effDuration = MAX.effDuration,
                                       OTS.data = OTS.data,
                                       OTS.threshold = OTS.threshold,
                                       OTS.effDuration = OTS.effDuration),
                      silent = TRUE)
            
            if (inherits(pd, "try-error")) {
                warning("No valid data to attach to the 'poisGPBayes0' object")
                pd <- NULL
            }
            
        } else {
            pd <- NULL
        }
        
    }
    
    if (!is.null(pd)) {
        pd <- threshData(threshold = threshold, data = pd, exceed = FALSE)   
    }
    
    cpn <-  checkParNames(parNames = colnames(MCMC), model = "poisGP",
                          all = FALSE)
    colnames(MCMC) <- parNames <- cpn$parNames[cpn$indIn]
    nSim <- nrow(MCMC)

    if (missing(MAP)) {
        MAP <- rep(NA, length(parNames))
        names(MAP) <- parNames           
    } else {
        
        if (length(MAP) != ncol(MCMC)) {
            stop("'MAP is expected to have length ncol(MCMC)")
        }
    }
    
    if (length(cpn$indOut) != 0) {
        
        if ((length(cpn$indOut) != 1) ||
            (cpn$parNames[cpn$indOut] != "lambda")) {
            stop("'mcmc' must have one column for each of the three 'poisgp'",
                 " parameters, possibly excepting the rate of the Poisson",
                 " process")
        }

        if (is.na(nOT)) {
            stop("When 'MCMC' does not embed the rate of the Possion Process ",
                 "the number 'nOT' of threshold exceedances must be given. ",
                 "'lambda' is assumed to be a posteriori independent of the ",
                 "GP parameters with gamma marginal posterior.")
        }


        if (!is.null(pd) && (pd$MAX$flag || pd$OTSflag)) {
            stop("Since the provided data embed MAX or OTS blocks ",
                 "the rate 'lambda' can not be a posteriori independent ",
                 "of the GP parameters so it must have a column in 'MCMC'.")
        }
        
        lambdaOut <- TRUE

        an <- a0 + nOT
        bn <- b0 + effDuration
        
        meanPost <- c("lambda" = an / bn, apply(MCMC, 2, mean))
        sdPost <- c("lambda" = sqrt(an) / bn, apply(MCMC, 2, sd))
        medianPost <- c("lambda" = qgamma(0.5, shape = an, rate = bn),
                        apply(MCMC, 2, median))
        
        if (length(MAP) == 2) {
            MAP <- c("lambda" = (an - 1) / bn, MAP)
        }
        
        covPost <- cov(MCMC)
        covPost <- cbind("lambda" = c(an / bn / bn, 0, 0),
                         rbind("lambda" = c(0, 0), covPost))
        comment <-  paste0("'lambda' and the GP parameter are a posteriori",
                           " independent") 
    } else {

        an <- bn <- NA
        lambdaOut <- FALSE
        meanPost <- apply(MCMC, 2, mean)
        sdPost <- apply(MCMC, 2, sd)
        medianPost <- apply(MCMC, 2, median)
        covPost <- cov(MCMC)
        comment <- NULL
        
    }
    
    res <- list(MCMC = MCMC,
                threshold = threshold,
                effDuration = effDuration,
                model = "poisGP",
                nOT = nOT,
                data = pd,
                an = an,
                bn = bn,
                meanPost = meanPost,
                sdPost = sdPost,
                medianPost = medianPost,
                MAP = MAP,
                covPost = covPost,
                lambdaOut = lambdaOut,
                comment = comment)
    
    class(res) <- c("poisGPBayes0", "Bayes0")
    res
    
}

# ****************************************************************************
##' Prediction of an Extreme-Value model of type Poisson-GP using MCMC
##' iterates. 
##'
##' The rate \code{lambda} can be omitted in the MCMC iterates of
##' \code{object}, which will be reported via the logical flag
##' \code{lambdaOut} set to \code{TRUE} among the elements of the
##' \code{object}. Then, the parameter \code{lambda} is assumed to be
##' a posteriori independent of the GP parameters. In this case, the
##' predictive distribution can be computed by using the predictive
##' distribution of the number of exceedances on the new period. This
##' specific treatement can be by-passed by using \code{addLamba =
##' TRUE}; in this case a column of MCMC iterates for \code{lambda} is
##' added to the matrix of MCMC iterates for the GP parameters, see
##' the \strong{Examples} section of the \code{\link{poisGPBayes0}}
##' page.
##'
##' @usage 
##' \method{predict}{poisGPBayes0}(object, newDuration = 1.0, prob,
##'         type = "RL",
##'         approx = FALSE,
##'         addLambda = !object$lambdaOut,
##'         trace = 0, ...)
##' 
##' @title Predictive Quantiles or Return Levels for a EV Model
##' of type Poisson-GP
##'
##' @param object An \code{poisGPBayes0} object, usually created by
##' using the eponymous creator function \code{\link{poisGPBayes0}}.
##'
##' @param newDuration The duration of the 'new' period for which the
##' maximum is to be predicted.
##' 
##' @param prob A vector of exceedance probabilities. The default
##' value contains probabilities such as \code{0.01} and \code{0.001}.
##' 
##' @param type The type of prediction wanted. Remind that the
##' \code{predict} method of the \strong{revdbayes} package uses this
##' argument to allow several types of predictions: density
##' \code{"d"}, quantile \code{"q"}, ...
##'
##' @param approx Logical. For the default \code{FALSE}, each value of
##' the tail quantile function is computed by zero-finding. For
##' \code{TRUE}, the quantiles are computed by using a fine of pairs
##' (argument, value) of the distribution function. This is likely to
##' be faster than \code{approx = FALSE} when \code{length(prob)} is
##' large.
##'
##' @param addLambda Logical. If the description of \code{object} that
##' \code{'lambda'} is idependent of the GP parameters, the exact
##' posterior of \code{'lambda'} will by default be used to compute
##' the predictive distribution. But if \code{addLambda} is then
##' passed with ist value set to \code{TRUE}, a new colum
##' corresponding to \code{'lambda'} is simply added to the MCMC
##' iterates and the computation is carried over by ignoring the
##' specific independence property. In practice this seems to make
##' little difference.
##'
##' @param trace Integer level of verbosity.
##' 
##' @param ... Not used yet.
##'
##' @seealso \code{\link[revdbayes]{predict.evpost}} in the
##' \strong{revdbayes} package , and the documentation of the creator
##' \code{\link{poisGPBayes0}}.
##' 
##' @return A data frame with the following columns
##' \item{NewDuration}{
##'
##' The duration of the "new" period on which the maximum is
##' predicted.
##'
##' }
##' \item{Prob}{
##'
##' A probability of exceedance.
##'
##' }
##' \item{Quantile}{
##'
##'  The return level corresponding to the probability.
##' 
##' }
##' 
##' The dataframe is given the S3 class \code{"predRL"} and it
##' receives several attributes such as the names of the factor
##' columns.
predict.poisGPBayes0 <- function(object,
                                 newDuration = 1.0,
                                 prob,
                                 type = "RL",
                                 approx = FALSE,
                                 addLambda = !object$lambdaOut,
                                 trace = 0,
                                 ...) {
    
    type <-  match.arg(type)
    nSim <- nrow(mcmc)
    
    if (missing(prob)) {
        prob <- c(0.800, 0.600, 0.500, 0.250, 0.200, 0.150, 0.100,
                  0.050, 0.030, 0.025, 0.020, 0.015,
                  0.010, 0.008, 0.005, 0.004, 0.003, 0.002, 0.001)
    } else {
        if (any(prob < 0.0) || any(prob > 1.0)) {
            stop("'prob' must be a numeric vector with values between ",
                 "0.0 and 1.0")
        }
        prob <- sort(prob, decreasing = TRUE)
    }
    nProb <- length(prob)
    T <- 1 / prob
    
    newDuration <- sort(newDuration)
    threshold <- object$threshold
    
    if (object$lambdaOut) {
        
        if (!addLambda) {
            
            if (trace) {
                cat("Using the representation as a mixture of CGP",
                    " Distributions\n")
            }
            
            ## =================================================================
            ## Find a interval to localise the largest and the smallest
            ## quantiles.
            ##
            ## Note that the posterior probability mass located at
            ## -Inf depends on the value of 'newDuration' but not of the
            ## GPD parameters. With the probability 'mass', there will
            ## be no event in the future period, so the exceedance
            ## probability corresponding to a finite Return Level must
            ## be < 1 - mass. For instance if the mass is 0.30 then
            ## the return levels corresponding to exceedances
            ## probabilities > 0.7 are infinite or NA.
            ## =================================================================

            for (iY in length(newDuration):1) {

                if (trace) {
                    cat(sprintf("newDuration = %5.0f\n", newDuration[iY]))
                }
                
                quant <- rep(NA, nProb)
                nu <- object$an
                p <- newDuration[iY] / (object$bn +  newDuration[iY])
                sigmaN <- (1 - p) / nu / p 
                xiN <- 1 / nu
                EN <- 1 / sigmaN
                IDN <- 1 + xiN / sigmaN
                mass <- nieve::pGPD2(1.0, scale = sigmaN, shape = xiN,
                                     lower.tail = FALSE)
                
                iProbMax <- min(seq_along(prob)[prob < 1 - mass])
                probMax <- prob[iProbMax]
                probTestU <- c(0.900, 0.950, 0.975, 0.980, 0.990, 0.999, 0.9999)

                if (trace) {
                    cat(sprintf("The probability mass is %5.2f iProbMax =%d\n",
                                mass, iProbMax))
                }
                
                ## =============================================================
                ## This is the posterior cumulative distribution function
                ## =============================================================
                
                FTilde <- function(q) {
                    res <- rep(NA, length(q))
                    for (i in seq_along(res)) {
                        res[i] <- mean(pCGPD(q = q[i], loc = threshold,
                                             scale = object$MCMC[ , 1],
                                             shape = object$MCMC[ , 2],
                                             scaleN = sigmaN, shapeN = xiN))
                    }
                    res
                }
                if (trace) {
                    cat("Finding quantile 'qU' corresponding to the smallest",
                        " exceedance probability, typically 0.001.\n")
                }
                OKU <-  FALSE
                yL <-  threshold
                for (i in seq_along(probTestU)) {
                    
                    if (!OKU) {
                        yU <- quantile(qCGPD(p = min(prob),
                                             loc = threshold,
                                             scale = object$MCMC[ , 1],
                                             shape = object$MCMC[ , 2],
                                             scaleN = sigmaN,
                                             shapeN = xiN,
                                             lower.tail = FALSE),
                                       prob = probTestU[i])
                        
                        qU <- try(uniroot(f = function(q){FTilde(q) - 1 + min(prob)},
                                          interval = c(yL, yU)),
                                  silent = TRUE)
                        if (!inherits(qU, "try-error") && qU$estim.prec < 1e-4) {
                            OKU <-  TRUE
                            qU <- qU$root
                            if (trace) {
                                cat("At trial number ", i, "'qU' is successfully",
                                    "localised.\n") 
                            }
                        } else {
                            if (trace) {
                                cat("The interval (", yL, ", ", yU, ") does not",
                                    "contain 'qU'\n") 
                            }
                        }
                    }
                    
                    if (OKU) break
                }
                
                qL <- try(uniroot(f = function(q) { FTilde(q) - 1 + probMax },
                                  interval = c(yL, yU)), silent = TRUE)
                
                if (!inherits(qL, "try-error") && qL$estim.prec < 1e-4) {
                    qL <- qL$root
                    if (trace) {
                        cat("'qL' is successfully localised.\n") 
                    }
                } else {
                    if (trace) {
                        cat("The interval (", yL, ", ", yU, ") does not contain",
                            "qL\n")
                    }
                }
                
                quant[iProbMax] <-  qL
                quant[nProb] <- qU
                
                if (nProb - iProbMax > 2) {
                    ind <- (nProb - 1):(iProbMax + 1)
                    if (approx) {
                        
                        yGrid <- seq(from = qL, to = qU, length.out = 200)
                        pGrid <- FTilde(yGrid)
                        quant[ind] <- approx(x = pGrid, y = yGrid, xout = 1 - prob[ind])$y
                    } else {
                        
                        for (i in ind) {
                            q <- try(uniroot(f = function(q) { FTilde(q) - 1 + prob[i] },
                                             interval = c(qL, qU)),
                                     silent = TRUE)
                            if (!inherits(q, "try-error") &&
                                q$estim.prec < 1e-2) {
                                OKL <-  TRUE
                                quant[i] <- q$root
                            } else {
                                stop("problem in quantile determination. Try",
                                     " 'approx = TRUE'")
                            }
                        }
                    }
                }
                    
                if (iY == length(newDuration)) {
                    res <- data.frame(NewDuration = newDuration[iY],
                                      Prob = prob,
                                      Quant = quant)
                } else {
                    res <-
                        dplyr::bind_rows(res,
                                         data.frame(NewDuration= newDuration[iY],
                                                    Prob = prob,
                                                    Quant = quant))
                }
            }
            
            res <- within(res, NewDuration <- factor(NewDuration)) 
            attr(res, "model") <- "poisGP"
            attr(res, "threshold") <- threshold
            attr(res, "obsDuration") <-  object$obsDuration
            attr(res, "blockDuration") <- NULL
            class(res) <-  c("predRL", "data.frame")
            
            return(res)
            
        }  else { 

            object$MCMC <- cbind(lambda = rgamma(n = nrow(object$MCMC),
                                     shape = object$an,
                                     rate = object$bn),
                                 object$MCMC)
            
        }
        
    } 
       
    ## ========================================================================
    ## Now, the standard case: there are three columns in mcmc,
    ## corresponing to the standard names c("lambda", "scale",
    ## "shape")
    ## ========================================================================

    thetaStar <- t(apply(object$MCMC,
                         MARGIN = 1, Renext::Ren2gev,
                         threshold = threshold,
                         w = 1.0,
                         distname.y = "GPD",
                         jacobian = 0))

    GEVpost <-  GEVBayes0(MCMC = thetaStar,
                          blockDuration = 1.0)
    
    res <- predict.GEVBayes0(GEVpost,
                             newDuration = newDuration,
                             type = "RL",
                             prob = prob,
                             approx = approx)
    
    ## ========================================================================
    ## change the attributes to avoid misunderstanding
    ## ========================================================================
    
    attr(res, "model") <- "poisgGP"
    attr(res, "threshold") <- threshold
    attr(res, "obsDuration") <- object$obsDuration
    attr(res, "block") <- NULL
   
    return(res)    
   
}

as.poisGPBayes0.poisGPBayes0 <- function(object, ...) {
    object
}
