#skrypt
#libraries -------
library(hsmm)
library(pROC)
#skrypt
pcname <- Sys.info()['nodename'] 
if(pcname=="piotr-tobit")
  setwd("~/Dropbox/doktorat/sekwencje_sygnalowe")
if(pcname=="MICHALKOMP")
  setwd("C:/Users/Michal/Dropbox/sekwencje_sygnalowe")
#wczytanie danych
source("wczytywanie_danych.R")
#wczytanie funkcji przygotowujacych dane treningowe
source("get_sig.R")
#algorytm viterbiego
source("myViterbi.R")

#test function

bihmm <- function(list_prot, aa_group) {
  t(vapply(list_prot, function(prot) {
    probka <- as.numeric(degenerate(prot, aa_group)[1:max.length])
    viterbi.res <- viterbi(probka, pipar, tp, od)
    viterbi_path <- viterbi.res[["path"]]
    c.site <- max(which(viterbi_path==3))
    if(c.site==-Inf) 
      c.site=length(probka)
    prob.signal <- viterbi.res$viterbi[length(probka), viterbi_path[length(probka)]]
    prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka, 0)
    prob.signal <- viterbi.res$viterbi[c.site, viterbi_path[c.site]]
    prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
    c(prob.signal = prob.signal, prob.non = prob.non, c.site = c.site)
  }, c(0, 0, 0)))
}


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

overall <- t4 #table(degenerate(unlist(analized_sequences), aa5))
overall.probs <- overall/sum(overall)          
overall.probs.log = log(overall.probs) #for viterbi

#setting params for hmm -------
additional.aminoacids = 10 #aminoacids choosen after cleavage site
pipar <- c(1,0,0,0)
tp <- matrix(c(0.814, 0.186, 0, 0,
               0, 0.91, 0.09, 0,
               0, 0, 0.78, 0.22,
               0, 0, 0, 1-1/additional.aminoacids), 4, byrow = TRUE)

od <- matrix(c((t1/sum(t1))[1:4],
               (t2/sum(t2))[1:4],
               (t3/sum(t3))[1:4],
               (t4/sum(t4))[1:4]), 4, byrow = TRUE)

# comparison of two models ------
numb.trials <- 400
max.length = 50
wyniki <- NULL
cuts <- NULL
testowane_bialka <- sample(1:length(analized_sequences), numb.trials, replace=FALSE)
for(numer_probki in 1:length(analized_sequences)){
  #probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)
  #                     [1:(all_nhc[numer_probki,4]+additional.aminoacids)])
  probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)[1:max.length])
  probka <- na.omit(probka)
  viterbi.res <- viterbi(probka, pipar, tp, od)
  viterbi_path <- viterbi.res$path
  c.site <- max(which(viterbi_path==3))
  if(c.site==-Inf) c.site=length(probka)
  cuts <- c(cuts, c.site)
  prob.signal <- viterbi.res[["viterbi"]][length(probka), viterbi_path[length(probka)]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka, 0)
  prob.signal <- viterbi.res$viterbi[c.site, viterbi_path[c.site]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
  wyniki <- rbind(wyniki, c(prob.signal, prob.non))
}

#negative ----
wyniki.not <- NULL
cuts.non <- NULL
testowane_bialka <- sample(1:length(euk_not), numb.trials, replace=FALSE)
tmp <- bihmm(euk_not[testowane_bialka], aa5)

# results ----------
cheat = 6 #we give handicap to signal peptides
sum(wyniki[,1]+cheat*cuts/max.length>wyniki[,2])/length(analized_sequences) #numb.trials
sum(wyniki.not[,1]+cheat*cuts.non/max.length<wyniki.not[,2])/length(euk_not) #numb.trials

wyniki <- bihmm(euk_not[testowane_bialka], aa5)

#naive way to compute AUC - does it make any sense?
a <- wyniki[,1]-wyniki[,2]
b <- wyniki.not[,1]-wyniki.not[,2]
standardized.probability <- exp(c(a,b)-max(b)) #possible bad idea but nothing better yet

auc(response=c(rep(0,length(a)), rep(1, length(b))), predictor=exp(c(a,b)-max(b)))
roc(response=c(rep(0,length(a)), rep(1, length(b))), predictor=exp(c(a,b)-max(b)), plot=T)
#' TO DO
#' 0. Przeanalizowanie czy wszystko jest poprawnie przeskalowane (czyli wyrugorwanie części handicapu)
#' 1. Model nieparametryczny z hsmm -> wyjąć prawdopodobieństwa
#' 2. Więcej modeli negatywnych
#' 


