#model hmm do sekwencji sygnalowej
#mamy cztery typy sygnałów
#1) aminokwasy z dodatnim ładunkiem K,R,H
#2) aminokwasy hydrofobowe V,I,L,M,F,W,C
#3) aminokwasy polarne bez ładunku S,T,N,Q
#4) reszta D,E,Y,G

#mamy cztery stany ukryte
#1) n-region
#2) h-region
#3) c-region
#4) region odcięcia

#stany ukryte produkują sygnały zgodnie z rozkładem wielomianowym
#zmiana stanów zgodnie z zadaną macierzą przejść (dofitowane parametry rozkładu wykładniczego)

#w takim modelu pakiet depmixS4 powinien sobie sam poradzić
require(depmixS4)

numer_probki <- 3
probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)[1:all_nhc[numer_probki,4]])
liczba_odpowiedzi <- length(unique(probka))

#tu muszą zostać wczytane dane
numStates <- 3
#szanse wyprodukowania sygnalow
state1 <- (t1/sum(t1))[1:liczba_odpowiedzi] #n-region
state2 <- (t2/sum(t2))[1:liczba_odpowiedzi]  #h-region
state3 <- (t3/sum(t3))[1:liczba_odpowiedzi]  #c-region

rModels <- list(
  list(
      GLMresponse(formula=probka~1, family=multinomial(), pstart=state1, fixed=c(TRUE, TRUE, TRUE, TRUE, TRUE))
  ),
  list(
    GLMresponse(formula=probka~1, family=multinomial(), pstart=state2, fixed=c(TRUE, TRUE, TRUE, TRUE, TRUE))
  ),
  list(
    GLMresponse(formula=probka~1, family=multinomial(), pstart=state3, fixed=c(TRUE, TRUE, TRUE, TRUE, TRUE))
  )
)

#transition probs, wyliczone z dlugosci regionów
transition <- list()
transition[[1]] <- transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0.814, 0.186, 0), fixed=c(TRUE,TRUE,TRUE))
transition[[2]] <- transInit(~1,nstates=numStates, family=multinomial("identity"), pstart=c(0, 0.91, 0.09), fixed=c(TRUE,TRUE,TRUE))
transition[[3]] <- transInit(~1, nstates=numStates, family=multinomial("identity"), pstart=c(0, 0, 1), fixed=c(TRUE,TRUE,TRUE))
transition

#zaczynamy od konkretnego stanu
instart=c(1,0,0)
inMod <- transInit(~1,ns=3,ps=instart,family=multinomial("identity"), data=data.frame(1), fixed=c(TRUE,TRUE,TRUE))

mod <- makeDepmix(response=rModels,transition=transition,prior=inMod)
summary(mod)
pars <- c(unlist(getpars(mod)))
pars
free <- c(0,0,0,rep(1,21))
fitted <- fit(mod, fixed=!free)
fitted@response
posterior(fitted)
all_nhc[numer_probki,]
summary(fitted)



