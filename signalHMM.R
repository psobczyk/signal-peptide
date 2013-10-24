#modelling signal peptides using HMM
setwd("~/Dropbox//Doktorat/sekwencje_sygnalowe/")
require(depmixS4)
source("signalClass.R")

#modeling part

#tu muszą zostać wczytane dane
numStates <- 4

state1 <- c(0.1, 0.4, 0.2, 0.3)
state2 <- c(0.1, 0.1, 0.2, 0.6)
state3 <- c(0.1, 0.4, 0.2, 0.3)
state4 <- c(0.1, 0.4, 0.2, 0.3)

rModels <- list(
  list(
    signal(observations,pstart=state1,fixed=c(TRUE))
  ),
  list(
    signal(observations,pstart=state2, fixed=c(TRUE))
  ),
  list(
    signal(observations,pstart=state3,fixed=c(TRUE))
  ),
  list(
    signal(observations,pstart=state4, fixed=c(TRUE)))
)

#transition probs
transition <- list()
transition[[1]] <- transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0.7, 0.3, 0, 0))
transition[[2]] <- transInit(~1,nstates=numStates, family=multinomial("identity"), pstart=c(0, 0.8, 0.2,0))
transition[[3]] <- transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0.7, 0.3, 0, 0))
transition[[4]] <- transInit(~1,nstates=numStates, family=multinomial("identity"), pstart=c(0, 0.8, 0.2,0))
transition

#zaczynamy od konkretnego stanu
instart=c(1,0,0,0)

inMod <- transInit(~1,ns=2,ps=instart,family=multinomial("identity"), data=data.frame(1))
mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=c(length(observations)))

fitted.mod <- fit(mod,  emc=em.control(rand=FALSE))
class(fitted.mod)

posterior(fitted.mod)
