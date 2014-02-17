#skrypt
#libraries -------
library(hsmm)

#loading data -------
pcname <- Sys.info()['nodename'] 
if(pcname=="piotr-tobit")
  setwd("~/Dropbox/doktorat/sekwencje_sygnalowe")

source("wczytywanie_danych.R")
source("get_sig.R") #wczytanie funkcji przygotowujacych dane treningowe

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

aa5 = list(positively.charged=c("K", "R", "H"),
           hydrofobic=c("V","I","L","M","F","W","C"),
           polar.uncharged=c("S", "T", "N", "Q"),
           rest=c("D","E","A","P","Y","G"))
t1 <- table(degenerate(n_region, aa5))
t2 <- table(degenerate(h_region, aa5))
t3 <- table(degenerate(c_region, aa5))
t4 <- table(degenerate(reszta, aa5))

#distribution of regions lengths
fit_poisson <- function(x){
  xfit<-floor(seq(min(x),max(x),length=60))
  yfit<-dpois(xfit,lambda=mean(x)) 
  h <- hist(x, probability=T)
  par(mfrow=c(2,1))
  plot(h$mids, h$density, type="l", xlab="długość regionu", ylab = "gęstość", col="red", 
       ylim=c(0, max(yfit, h$density)), main="Zgodność z rozkładem Poissona")
  lines(xfit, yfit, col="blue", lwd=2, xaxt="n",yaxt="n",xlab="",ylab="") 
  legend("topright",col=c("red","blue"),lty=1,legend=c("data","fitted"),cex=0.6)
  x = table(x)
  k=as.numeric(names(x))
  plot(k,log(x)+lfactorial(k), col="red",pch="o", ylab="częstość", main="Poissoness plot")
  a <-lm(log(x)+lfactorial(k)~k)
  abline(coef=a$coefficients, col="blue")
  legend("topleft",lty=1,pch=c("o", ""), col=c("red","blue"),legend=c("data","fitted"),cex=0.6)
}
fit_poisson(lengths[,1])
fit_poisson(lengths[,2])
fit_poisson(lengths[,3])

#setting params -------
additional_margin = 10
pipar <- c(1,0,0,0)
tpmpar <- matrix(c(0, 1, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 0, 1,
                   0.01, 0, 0, 0), 4, byrow = TRUE)
rdpar <- list(lambda = c(5.360593, 11.56915, 4.554336,additional_margin))
odpar <- list(m = matrix(c((t1/sum(t1))[1:4],
                           (t2/sum(t2))[1:4],
                           (t3/sum(t3))[1:4],
                           (t4/sum(t4))[1:4]), 4, byrow = TRUE))

#computing viterbi path ------
procent_rozpoznania <- NULL
testowane_bialka <- sample(1:length(analized_sequences),100, replace=FALSE)
for(numer_probki in testowane_bialka){
  probka <- as.numeric(degenerate(analized_sequences[[numer_probki]], aa5)[1:(all_nhc[numer_probki,4]+additional_margin)])
  viterbi_path <- hsmm.viterbi(probka,od = "mult", rd = "pois",
                                 pi.par = pipar, tpm.par = tpmpar,
                                 od.par = odpar, rd.par = rdpar)
  viterbi_path <- viterbi_path$path[1:all_nhc[numer_probki,4]]
  expected <- c(rep(1,all_nhc[numer_probki,2]-1),rep(2,all_nhc[numer_probki,3]-all_nhc[numer_probki,2]),rep(3,all_nhc[numer_probki,4]-all_nhc[numer_probki,3]+1))
  procent_rozpoznania <- c(procent_rozpoznania, sum(viterbi_path==expected)/length(viterbi_path))
}

#results ----------
mean(procent_rozpoznania)
hist(procent_rozpoznania,probability=T,breaks=20)
