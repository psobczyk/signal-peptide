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
if(pcname=="MICHALKOMP")
  setwd("D:/michal/doktorat/signal-peptide")

#wczytanie funkcji przygotowujacych dane treningowe
source("get_sig.R")
#algorytm viterbiego
source("myViterbi.R")

bihmm <- function(list_prot, aa_group) {
  t(vapply(list_prot, function(prot) {
    probka <- as.numeric(degenerate(prot, aa_group)[1:max.length])
    probka <- na.omit(probka)
    viterbi.res <- viterbi(probka, pipar, tp, od)
    viterbi_path <- viterbi.res[["path"]]
    c.site <- ifelse(any(viterbi_path==3), max(which(viterbi_path==3)), length(probka))
    prob.signal <- viterbi.res[["viterbi"]][c.site, viterbi_path[c.site]]
    prob.non <- Reduce(function(x, y) x + overall.probs.log[y], probka[1:c.site], 0)
    c(prob.signal = prob.signal, prob.non = prob.non, c.site = c.site)
  }, c(0, 0, 0)))
}

aa1 = list('1' = c("G", "A", "P", "V", "L", "I", "M"), 
           '2' = c("F", "W", "Y"), 
           '3' = c("K", "R", "H"), 
           '4' = c("D", "E"), 
           '5' = c("S", "T", "C", "N", "Q"))

aa2 = list('1' = c("G", "A", "P", "V", "L", "I", "M", "F"), 
           '2' =c("K", "R", "H"), 
           '3' =c("D", "E"), 
           '4' =c("S", "T", "C", "N", "Q", "Y", "W"))

aa3 = list('1' = c("G", "A", "P", "V", "L", "I", "M", "F", "W", "S", "T", "C", "N", "Q", "Y"), 
           '2' = c("D", "E"), 
           '3' = c("K", "R", "H"))

aa4 = list('1' = c("G", "A", "V", "L", "I"), 
           '2' = c("M", "C"), 
           '3' = c("F", "W", "Y"), 
           '4' = c("S", "T", "N", "Q"), 
           '5' = c("K", "R", "H"), 
           '6' = c("D", "E"), 
           '7' ="P")

aa5 = list('1' = c("K", "R", "H"),
           '2' = c("V","I","L","M","F","W","C"),
           '3' = c("S", "T", "N", "Q"),
           '4' = c("D","E","A","P","Y","G"))

aa6 <- list(`1` = c("A", "C", "Q", "E"), 
            `2` = c("R", "H",  "K"), 
            `3` = c("N", "D", "G", "P", "S", "T"), 
            `4` = c("I", "L",  "M", "F", "W", "Y", "V"))

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

numb.trials <- 400
testowane_bialka <- sample(1:length(analized_sequences), numb.trials, replace=FALSE)
testowane_bialka <- sample(1:length(euk_not), numb.trials, replace=FALSE)
prot_pos = analized_sequences
prot_neg = euk_not[testowane_bialka]

cheat <- 6 #we give handicap to signal peptides
test_pos <- lapply(read.fasta("test_pos.fasta", seqtype = "AA"), toupper)
test_neg <- lapply(read.fasta("test_neg.fasta", seqtype = "AA"), toupper)
test_neg <- test_neg[sapply(test_neg[1:100], length) > 100][1:100]
real_cl <- sapply(read_uniprot("test_pos.txt", euk = TRUE), function(i) 
  attr(i, "sig")[2])


############tu start funkcji
bihmm_test <- function(prot_pos, prot_neg, test_pos, test_neg, real_cl, cheat) {

t1 <- table(degenerate(n_region, aa_group))
t2 <- table(degenerate(h_region, aa_group))
t3 <- table(degenerate(c_region, aa_group))
t4 <- table(degenerate(reszta, aa_group))

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
wyniki <- bihmm(prot_pos, aa5)

#negative ----

wyniki.not <- bihmm(prot_neg, aa5)

# results ----------

#sum(wyniki[,1] + cheat * cuts/max.length > wyniki[,2])/length(analized_sequences) #numb.trials
#sum(wyniki.not[,1] + cheat * cuts.non/max.length <wyniki.not[,2])/length(euk_not) #numb.trials

#naive way to compute AUC - does it make any sense?
a <- wyniki[,1] - wyniki[,2]
b <- wyniki.not[,1] - wyniki.not[,2]
#standardized.probability <- exp(c(a,b) - max(b)) #possible bad idea but nothing better yet

#auc(response=c(rep(0,length(a)), rep(1, length(b))), predictor=exp(c(a,b)-max(b)))
#roc(response=c(rep(0,length(a)), rep(1, length(b))), predictor=exp(c(a,b)-max(b)), plot=T)
#' TO DO
#' 0. Przeanalizowanie czy wszystko jest poprawnie przeskalowane (czyli wyrugorwanie części handicapu)
#' 1. Model nieparametryczny z hsmm -> wyjąć prawdopodobieństwa
#' 2. Więcej modeli negatywnych
#' 
#

# quick test


wyniki <- bihmm(test_pos, aa5)
wyniki.not <- bihmm(test_neg, aa5)
a <- wyniki[,1]-wyniki[,2]
b <- wyniki.not[,1]-wyniki.not[,2]
standardized.probability <- exp(c(a,b)-max(b)) #possible bad idea but nothing better yet
list(auc = auc(response=c(rep(0,length(a)), rep(1, length(b))), predictor=exp(c(a,b)-max(b))),
     roc = roc(response=c(rep(0,length(a)), rep(1, length(b))), predictor=exp(c(a,b)-max(b)), plot=FALSE),
     mse_cl = sqrt(mean((real_cl[rownames(wyniki)] - wyniki[, 3])^2, na.rm = TRUE)))
}

