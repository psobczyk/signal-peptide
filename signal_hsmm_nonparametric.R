#skrypt
#libraries -------
library(hsmm)
library(pROC)
#skrypt
pcname <- Sys.info()['nodename'] 
if(pcname=="piotr-tobit")
  setwd("~/Dropbox/doktorat/sekwencje_sygnalowe")

#wczytanie danych
source("wczytywanie_danych.R")
#wczytanie funkcji przygotowujacych dane treningowe
source("get_sig.R")
source("mydurationViterbi.R")

#building training set ----
analized_sequences <- speuk[with_sig]
euk_not <- read.fasta("euk_not.fasta", seqtype = "AA")

#połączone regiony
n_region <- NULL
h_region <- NULL
c_region <- NULL
reszta <- NULL

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

#setting params -------
additional_margin = 10
pipar <- c(1,0,0,0)
tpmpar <- matrix(c(0, 1, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 0, 1,
                   0, 0, 0, 0), 4, byrow = TRUE)
od <- matrix(c((t1/sum(t1))[1:4],
               (t2/sum(t2))[1:4],
               (t3/sum(t3))[1:4],
               (t4/sum(t4))[1:4]), 4, byrow = TRUE)


# comparison of two models ------
numb.trials <- 100
max.length = 50
wyniki <- NULL
cuts <- NULL
testowane_bialka <- sample(1:length(analized_sequences), numb.trials, replace=FALSE)
for(numer_probki in testowane_bialka){ #1:length(analized_sequences)){
  #probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)
  #                     [1:(all_nhc[numer_probki,4]+additional.aminoacids)])
  probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)[1:max.length])
  probka <- na.omit(probka)
  viterbi.res <- duration.viterbi(probka, pipar, tpmpar, od, params=params)
  viterbi_path <- viterbi.res$path
  c.site <- max(which(viterbi_path==3))
  if(c.site==-Inf) c.site=length(probka)
  cuts <- c(cuts, c.site)
  prob.signal <- viterbi.res$viterbi[length(probka), viterbi_path[length(probka)]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka, 0)
  prob.signal <- viterbi.res$viterbi[c.site, viterbi_path[c.site]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
  wyniki <- rbind(wyniki, c(prob.signal, prob.non))
}

plot.Viterbi.results(viterbi_path, analized_sequences[[numer_probki]])

#negative ----
wyniki.not <- NULL
cuts.non <- NULL
testowane_bialka <- sample(1:length(euk_not),numb.trials, replace=FALSE)
for(numer_probki in testowane_bialka){ #1:length(euk_not)){
  probka <- as.numeric(degenerate(euk_not[[numer_probki]], aa5)[1:max.length])
  viterbi.res <- duration.viterbi(probka, pipar, tpmpar, od, params=params)
  viterbi_path <- viterbi.res$path
  c.site <- max(which(viterbi_path==3))
  if(c.site==-Inf) c.site=length(probka)
  cuts.non <- c(cuts.non, c.site)
  prob.signal <- viterbi.res$viterbi[length(probka), viterbi_path[length(probka)]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka, 0)
  prob.signal <- viterbi.res$viterbi[c.site, viterbi_path[c.site]]
  prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
  wyniki.not <- rbind(wyniki.not, c(prob.signal, prob.non))
}

# results ----------
cheat = 5 #we give handicap to signal peptides
sum(wyniki[,1]+cheat*cuts/max.length>wyniki[,2])/numb.trials #length(analized_sequences) 
sum(wyniki.not[,1]+cheat*cuts.non/max.length<wyniki.not[,2])/numb.trials #length(euk_not) 

#naive way to compute AUC - does it make any sense?
a <- wyniki[,1]-wyniki[,2]
b <- wyniki.not[,1]-wyniki.not[,2]
standardized.probability <- exp(c(a,b)-max(b)) #possible bad idea but nothing better yet

auc(response=c(rep(0,length(a)), rep(1, length(b))), predictor=exp(c(a,b)-max(b)))
roc(response=c(rep(0,length(a)), rep(1, length(b))), predictor=exp(c(a,b)-max(b)), plot=T)

