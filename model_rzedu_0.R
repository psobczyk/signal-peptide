#model markowa rzędu 0 (czestości wystepowania aminikwasów)

source("main.R")

#uruchamiamy ukryte łańcuchy
require(depmixS4)
source("run_model.R")

procent_rozpoznania <- NULL
testowane_bialka <- sample(1:length(analized_sequences),100, replace=FALSE)
for(numer_probki in testowane_bialka){
  probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)[1:(all_nhc[numer_probki,4]+30)])
  fitted.model <- uruchom_model(probka)
  viterbi_path <- fitted.model@posterior[1:all_nhc[numer_probki,4],1]
  expected <- c(rep(1,all_nhc[numer_probki,2]-1),rep(2,all_nhc[numer_probki,3]-all_nhc[numer_probki,2]),rep(3,all_nhc[numer_probki,4]-all_nhc[numer_probki,3]+1))
  procent_rozpoznania <- c(procent_rozpoznania, sum(viterbi_path==expected)/length(viterbi_path))
}

#bardzo prosta statystyka opisowa ;)
mean(procent_rozpoznania)

wynik <- cbind(probka, fitted.model@posterior)
