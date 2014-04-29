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

params <- NULL
max.length = 32

n_region_lengths <- table(lengths[,1])
max.nRegion.length = which.min(1:max(as.numeric(names(n_region_lengths))) %in% names(n_region_lengths))-1
nRegion.lengths = n_region_lengths[1:max.nRegion.length]
params <- nRegion.lengths/sum(nRegion.lengths)

h_region_lengths <- table(lengths[,2])[-1]
max.hRegion.length = which.min(1:max(as.numeric(names(h_region_lengths))) %in% names(h_region_lengths))-1
max.hRegion.length = 25
hRegion.lengths = h_region_lengths[1:max.hRegion.length]
params <- cbind(params, c(hRegion.lengths/sum(hRegion.lengths), rep(0, max.length-max.hRegion.length)))

c_region_lengths <- table(lengths[,3])
max.cRegion.length = which.min(2:max(as.numeric(names(c_region_lengths))) %in% names(c_region_lengths))-1
cRegion.lengths = c_region_lengths[1:max.cRegion.length]
params <- cbind(params, c(0,cRegion.lengths/sum(cRegion.lengths), rep(0, max.length-max.cRegion.length-1)))
params <- cbind(params, rep(1/max.length, max.length))

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
t2 <- table(degenerate(h_region, aa5))
t3 <- table(degenerate(c_region, aa5))
t4 <- table(degenerate(reszta, aa5))

