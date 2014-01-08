#hmm model - trzy stany zwiazane z regionami n,h i c

uruchom_model2 <- function(probka){
  #tu muszą zostać wczytane dane
  numStates <- 3
  liczba_odpowiedzi <- length(unique(probka))
  responses <- 1:16 %in% probka
  #szanse wyprodukowania sygnalow
  state1 <- transition1[responses] #n-region
  state2 <- transition2[responses]  #h-region
  state3 <- transition3[responses]  #c-region
  #state1 <- transition1[1:liczba_odpowiedzi] #n-region
  #state2 <- transition2[1:liczba_odpowiedzi]  #h-region
  #state3 <- transition3[1:liczba_odpowiedzi]  #c-region
  
  rModels <- list(
    list(
      GLMresponse(formula=probka~1, family=multinomial(), pstart=state1, fixed=rep(TRUE,100)[1:liczba_odpowiedzi])
    ),
    list(
      GLMresponse(formula=probka~1, family=multinomial(), pstart=state2, fixed=rep(TRUE,100)[1:liczba_odpowiedzi])
    ),
    list(
      GLMresponse(formula=probka~1, family=multinomial(), pstart=state3, fixed=rep(TRUE,100)[1:liczba_odpowiedzi])
    )
  )
  
  #transition probs, wyliczone z dlugosci regionów
  transition <- list()
  transition[[1]] <- transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0.82, 0.18, 0), fixed=c(TRUE,TRUE,TRUE))
  transition[[2]] <- transInit(~1,nstates=numStates, family=multinomial("identity"), pstart=c(0, 0.9, 0.1), fixed=c(TRUE,TRUE,TRUE))
  transition[[3]] <- transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0, 0, 1), fixed=c(TRUE,TRUE,TRUE))
  
  #zaczynamy od konkretnego stanu
  instart=c(1,0,0)
  inMod <- transInit(~1,ns=3,ps=instart,family=multinomial("identity"), data=data.frame(1), fixed=c(TRUE,TRUE,TRUE, TRUE))
  
  mod <- makeDepmix(response=rModels,transition=transition,prior=inMod)
  #summary(mod)
  #pars <- c(unlist(getpars(mod)))
  #length(pars)
  free <- c(0,0,0,rep(1,9), rep(c(1,1,1), liczba_odpowiedzi))
  fitted <- fit(mod, fixed=!free)
  #fitted@response
  #cbind(probka,posterior(fitted))
  #all_nhc[numer_probki,]
  #summary(fitted)
  return(fitted)
}
