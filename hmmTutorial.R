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
setMethod("fit","exgaus",
          function(object,w) {
            if(missing(w)) w <- NULL
            y <- object@y
            fit <- gamlss(y~1,weights=w,family=exGAUS(),
                          control=gamlss.control(n.cyc=100,trace=FALSE),
                          mu.start=object@parameters$mu,
                          sigma.start=exp(object@parameters$sigma),
                          nu.start=exp(object@parameters$nu))
            pars <- c(fit$mu.coefficients,fit$sigma.coefficients,fit$nu.coefficients)
            object <- setpars(object,pars)
            object
          }
)
setMethod("predict","exgaus",
          function(object) {
            ret <- object@parameters$mu
            return(ret)
          }
)


coinModel <- coinToss(c(1,0,1,1), c(0.3))
dens(coinModel)
summary(coinModel)
str(coinModel)
coinModel@parameters

fm <- fit(coinModel, emc=em.control(rand=FALSE))
summary(fm)

data(speed)
summary(speed)
speed$corr
glm(speed$rt~1)