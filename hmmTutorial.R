#modelling using HMM

require(depmixS4)

#ok, spróbujemy stworzyć prosty model dwustanowy - rzuty moneta
# musi być w formie response tzn musi implementowac metody dens, predict, (?fit)

setClass("coinToss",contains="response")
setGeneric("coinToss", function(y, pstart = NULL, fixed = NULL, ...) standardGeneric("coinToss"))

setMethod("coinToss",
          signature(y="ANY"),
          function(y,pstart=NULL,fixed=NULL, ...) {
            y <- matrix(y,length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 1
            if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
            if(!is.null(pstart)) {
              if(length(pstart)!=npar) stop("length of ’pstart’ must be ",npar)
              parameters$heads <- pstart[1]
              #parameters$sigma <- log(pstart[2])
              #parameters$nu <- log(pstart[3])
            }
            mod <- new("coinToss",parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
            mod
          }
    )

setMethod("show","coinToss",
          function(object) {
            cat("Model of coin tossing \n")
            cat("Parameters: \n")
            cat("Prob of heads: ", object@parameters$heads, "\n")
          }
)
setMethod("dens","coinToss",
          function(object,log=FALSE) {
            object@parameters$heads^object@y*(1-object@parameters$heads)^(1-object@y)
          }
        )
setMethod("getpars","response",
          function(object,which="pars",...) {
            switch(which,
                   "pars" = {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            return(pars)
          }
        )

setMethod("setpars","coinToss",
          function(object, values, which="pars", ...) {
            npar <- npar(object)
            if(length(values)!=npar) stop("length of ’values’ must be",npar)
            # determine whether parameters or fixed constraints are being set
            nms <- names(object@parameters)
            switch(which,
                   "pars"= {
                     object@parameters$heads <- values[1]
                     #object@parameters$sigma <- values[2]
                     #object@parameters$nu <- values[3]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          }
)
setMethod("fit","coinToss",
          function(object,w) {
            if(missing(w)) w <- NULL
            y <- object@y
            object
          }
)
setMethod("predict","coinToss",
          function(object) {
            if(object@parameters$heads>0.5){
              ret <- 1
            }
            else ret <- 0
            return(ret)
          }
)


observations <- c(0,0,1,1,1,0,0,0,1)

coinModel <- coinToss(observations, c(0.3))
dens(coinModel)
summary(coinModel)
str(coinModel)
coinModel@parameters

rModels <- list(
  list(
    coinToss(observations,pstart=c(0.9),fixed=c(TRUE))
  ),
  list(
    coinToss(observations,pstart=c(0.1), fixed=c(TRUE)))
)
rModels[[1]][[1]]@fixed
#transition probs
transition <- list()
transition[[1]] <- transInit(~1, nstates=2, family=multinomial("identity"), pstart=c(0.2, 0.8), fixed=c(TRUE,TRUE))
transition[[2]] <- transInit(~1,nstates=2, family=multinomial("identity"), pstart=c(0.7,0.3),fixed=c(TRUE,TRUE))
transition[[1]]@parameters
transition[[2]]@parameters

instart=c(0.5,0.5)
inMod <- transInit(~1,ns=2,ps=instart,family=multinomial("identity"), data=data.frame(1), fixed=c(TRUE, TRUE))
mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=c(length(observations)))
mod@transition
fitted.mod <- fit(mod,emc=em.control(rand=FALSE))
fitted.mod@transition
fitted.mod@homogeneous
posterior(fitted.mod)
