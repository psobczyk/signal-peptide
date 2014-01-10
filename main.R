#skrypt
pcname <- Sys.info()['nodename'] 
if(pcname=="piotr-tobit")
  setwd("~/Dropbox/doktorat/sekwencje_sygnalowe")

#wczytanie danych
source("wczytywanie_danych.R")
#wczytanie funkcji przygotowujacych dane treningowe
source("get_sig.R")
#połączone regiony
n_region <- NULL
h_region <- NULL
c_region <- NULL
reszta <- NULL
analized_sequences <- speuk[with_sig]

for(i in 1:3644){
  n_region <- c(n_region, analized_sequences[[i]][1:lengths[i,1]]) 
  h_region <- c(h_region, analized_sequences[[i]][all_nhc[i,2]:(all_nhc[i,3]-1)])
  c_region <- c(c_region, analized_sequences[[i]][all_nhc[i,3]:(all_nhc[i,4]-1)])
  reszta <- c(reszta, analized_sequences[[i]][all_nhc[i,4]:(length(analized_sequences[[i]]))])
}

#oszacowanie parametrow rozkladow wykladniczych
#library("fitdistrplus") 
#fitdist(lengths[,1], "exp", method ="mle")
#fitdist(lengths[,2], "exp", method ="mle")
#fitdist(lengths[,3], "exp", method ="mle")

aa5 = list(positively.charged=c("K", "R", "H"),
           hydrofobic=c("V","I","L","M","F","W","C"),
           polar.uncharged=c("S", "T", "N", "Q"),
           rest=c("D","E","A","P","Y","G"))

t1 <- table(degenerate(n_region, aa5))
t1/sum(t1)
t2 <- table(degenerate(h_region, aa5))
t2/sum(t2)
t3 <- table(degenerate(c_region, aa5))
t3/sum(t3)
t4 <- table(degenerate(reszta, aa5))
t4/sum(t4)

