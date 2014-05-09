#skrypt
#libraries -------
library(hsmm)

#skrypt
pcname <- Sys.info()['nodename'] 
if(pcname=="piotr-tobit")
  setwd("~/Dropbox/doktorat/sekwencje_sygnalowe")

#wczytanie danych
source("wczytywanie_danych.R")
#wczytanie funkcji przygotowujacych dane treningowe
source("get_sig.R")
#algorytm viterbiego
source("myViterbi.R")

#building training set ----
analized_sequences <- speuk[with_sig]
euk_not <- read.fasta("euk_not.fasta", seqtype = "AA")

n_region <- NULL
h_region <- NULL
c_region <- NULL
reszta <- NULL
for(i in 1:length(analized_sequences)){
  n_region <- c(n_region, analized_sequences[[i]][1:lengths[i,1]]) 
  h_region <- c(h_region, analized_sequences[[i]][all_nhc[i,2]:(all_nhc[i,3]-1)])
  c_region <- c(c_region, analized_sequences[[i]][all_nhc[i,3]:(all_nhc[i,4]-1)])
  reszta <- c(reszta, analized_sequences[[i]][all_nhc[i,4]:(length(analized_sequences[[i]]))])
}

aa5 = list(positively.charged=c("K", "R", "H"),
           hydrofobic=c("V","I","L","M","F","W","C"),
           polar.uncharged=c("S", "T", "N", "Q"),
           rest=c("D","E","A","P","Y","G"))

t1 <- table(degenerate(n_region, aa5))
t2 <- table(degenerate(h_region, aa5))
t3 <- table(degenerate(c_region, aa5))
t4 <- table(degenerate(reszta, aa5))

overall <- table(degenerate(unlist(analized_sequences), aa5))
overall.probs <- overall/sum(overall)          
overall.probs.log = log(overall.probs) #for viterbi

#setting params for hmm -------
additional.aminoacids = 10 #aminoacids choosen after cleavage site
pipar <- c(1,0,0,0)
tp <- matrix(c(0.814, 0.186, 0, 0,
               0, 0.91, 0.09, 0,
               0, 0, 0.78, 0.22,
               0, 0, 0, 1), 4, byrow = TRUE)

od <- matrix(c((t1/sum(t1))[1:4],
               (t2/sum(t2))[1:4],
               (t3/sum(t3))[1:4],
               (t4/sum(t4))[1:4]), 4, byrow = TRUE)

# comparison of two models ------
numb.trials <- 400
wyniki <- NULL
testowane_bialka <- sample(1:length(analized_sequences), numb.trials, replace=FALSE)
for(numer_probki in testowane_bialka){
  probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)
                       [1:(all_nhc[numer_probki,4]+additional.aminoacids)])
  probka <- na.omit(probka)
  viterbi.res <- viterbi(probka, pipar, tp, od)
  viterbi_path <- viterbi.res$path
  prob.signal <- viterbi.res$viterbi[length(probka), viterbi_path[length(probka)]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka, 0)
  wyniki <- rbind(wyniki, c(prob.signal, prob.non))
}

#negative ----
max.length = 50
wyniki.not <- NULL
testowane_bialka <- sample(1:length(euk_not),numb.trials, replace=FALSE)
for(numer_probki in testowane_bialka){
  probka <- as.numeric(degenerate(euk_not[[numer_probki]], aa5)[1:max.length])
  viterbi.res <- viterbi(probka, pipar, tp, od)
  viterbi_path <- viterbi.res$path
  prob.signal <- viterbi.res$viterbi[length(probka), viterbi_path[length(probka)]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka, 0)
  wyniki.not <- rbind(wyniki.not, c(prob.signal, prob.non))
}

# results ----------
cheat = 3.5 #we give handicap to signal peptides
sum(wyniki[,1]+cheat>wyniki[,2])/numb.trials
sum(wyniki.not[,1]+cheat<wyniki.not[,2])/numb.trials

#' TO DO
#' 0. Przeanalizowanie czy wszystko jest poprawnie przeskalowane
#' 1. Model nieparametryczny z hsmm -> wyjąć prawdopodobieństwa
#' 2. Więcej modeli negatywnych
#' 
