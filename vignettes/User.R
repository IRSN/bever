## ----global_options, include=FALSE---------------------------------------
library(knitr)
knitr::opts_chunk$set(
    cache = FALSE,
    fig.width = 6, fig.height = 4, fig.path = "Rgraphics/")
KABLE <- TRUE
library(bever)
library(ismev)


## ----PoorMan_GEV1, message=FALSE-----------------------------------------
priorGEV1 <- set_prior(prior = "flatflat", model = "gev")
postGEV1 <- rpost_rcpp(n = 10000, model = "gev", prior = priorGEV1, data = portpirie)
class(postGEV1)
head(postGEV1$sim_vals, n = 3)


## ----PoorMan_GEV2, message=FALSE-----------------------------------------
postGEV0 <- GEVBayes0(MCMC = postGEV1$sim_vals)
class(postGEV0)


## ----PoorMan_GEV3, message=FALSE-----------------------------------------
summary(postGEV0)
coef(postGEV0)
vcov(postGEV0)


## ----PoorMan_poisGP_1, message=FALSE-------------------------------------
data(rainfall)
rainfall2 <- rainfall[!is.na(rainfall)]
u <- 40  ## threshold (mm)
w <- 57  ## obs duration (year)
nOT <- sum(rainfall2 > u)
nSim <- 10000
priorGP1 <- set_prior(prior = "flatflat", model = "gp")
postGP1 <- rpost_rcpp(n = nSim, model = "gp", prior = priorGP1,
                      data = rainfall2, thresh = u)
class(postGP1)
head(postGP1$sim_vals, n = 3)


## ----PoorMan_poisGP_2, message=FALSE-------------------------------------
MCMCpoisGP <- cbind(lambda = rgamma(nSim, shape = 1 + nOT, rate = 0 + w),
                    postGP1$sim_vals)
postGP0 <- poisGPBayes0(MCMC = MCMCpoisGP, threshold = u, effDuration = w) 
summary(postGP0)
coef(postGP0)


## ----PoorMan_poisGP_3, message=FALSE-------------------------------------
postGP0b <- poisGPBayes0(MCMC = postGP1$sim_vals, threshold = u,
                         effDuration = w, nOT = nOT) 
summary(postGP0b)
coef(postGP0b)


## ----GaronneML, message=FALSE--------------------------------------------
library(bever)
y <- Garonne$OTdata$Flow
u <- Garonne$OTinfo$threshold
w <- Garonne$OTinfo$effDuration
gfitML <- poisGP(data = y, threshold = u, effDuration = w)
class(gfitML)
coef(gfitML)
confint(gfitML)
autoplot(gfitML) +
    ggtitle("Fitted Poisson-GP model for \"Garonne\"") + theme_gray()



## ----GaronneML_RL1, message=FALSE----------------------------------------
gRLML <- RL(gfitML, confintMethod = "proflik", trace = 0)
class(gRLML)
tail(gRLML, n = 6)
autoplot(gRLML) +
    ggtitle("Classical RL plot : frequentist POT for Garonne") + theme_gray()



## ----GaronneBayes_RL2, message=FALSE-------------------------------------
gfitBayes <- poisGPBayes(data = y, threshold = u, effDuration = w)
class(gfitBayes)
coef(gfitBayes)


## ----GaronneBayes_RL3, message=FALSE-------------------------------------
gRLBayes <- RL(gfitBayes, trace = 0)
class(gRLBayes)
tail(gRLBayes, n = 6)
autoplot(gRLBayes) +
    ggtitle("Classical RL plot : Bayesian POT for Garonne") + theme_gray()



## ----portpirieML_RL1, message=FALSE, warning=FALSE-----------------------
fitp <- poisGP(MAX.data = as.list(portpirie),
               MAX.effDuration=rep(1, length(portpirie)))
autoplot(fitp) + theme_gray()


## ----portpirieBayes_RL1, message=FALSE-----------------------------------
prior <- set_prior(prior = "flatflat", model = "gev")
post <- rpost_rcpp(n = 10000, model = "gev", prior = prior, data = portpirie)
MAPp <- post$f_mode
names(MAPp) <- c("loc", "scale", "shape")
pfitGEV0 <- GEVBayes0(MCMC = post$sim_vals, yMax = portpirie, MAP = MAPp)


## ----portpirieBayes_RL2, message=FALSE-----------------------------------
myRL <- RL(pfitGEV0) 
autoplot(myRL) + ggtitle("Classical RL plot: Bayesian GEV for portpirie") +
    theme_gray()


## ----gBayes_predict, message=FALSE---------------------------------------
gpredBayes <- predict(gfitBayes)
autoplot(gpredBayes) +
    ggtitle("Predictive RL plot : Bayesian POT for Garonne") + theme_gray()



## ----pBayes_predict, message=FALSE---------------------------------------
ppredGEV0 <- predict(pfitGEV0)
autoplot(ppredGEV0) +
    ggtitle("Predictive RL plot : Bayesian GEV for portpirie") + theme_gray()


