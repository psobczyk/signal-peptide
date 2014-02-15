#setwd("~/Dropbox/Doktorat/sekwencje_sygnalowe/")
#source("wczytywanie_danych.R")

con <- function(x, n){
  result <- NULL
  for (i in 1:n){
    result[i] <- sum(x==i)/length(x) 
  }
  result
} #computing probs for n groups (it must be exclusive)

distance <- function(con1, con2){
  sum((con1-con2)^2/con1)
} #distance we want to maximise while choping amino acid sequence into two

distribution <- function(rawSeqs, aagroups){
  N <- length(rawSeqs)
  seqs <- vector(mode="list", length=N)
  for(j in 1:N){
    a = rawSeqs[j][[1]]
    aas = NULL
    for (i in 1:length(a)){
      aas = c(aas, a[i])
    }
    seqs[[j]] <- aas
  }

  #now we have list of vectors representing sequences of aminoacids
  seqs <- lapply(seqs, FUN=function(x) degenerate(toupper(x), aagroups))
  #now our sequences are degenarated ti certain aminoacid groups
  #print(seqs[50:100])
  cutoff <- NULL
  for(seq in seqs){
    n1 <- 15
    n2 <- 100
    probs <- NULL
    for (i in n1:(min(n2, length(seq)-2))){
      con1 = con(seq[1:i], length(aagroups)) + 0.1
      con2 = con(seq[(i+1):length(seq)], length(aagroups)) + 0.1
      probs <- c(probs, distance(con1, con2))
    } #we compute the most likely place for end of signal sequence
    #print(n1+which.max(unlist(probs)))
    #print(probs)
    cutoff <- c(cutoff, n1-1+which.max(probs))
    #print(paste((n1+which.max(na.omit(probs))), " ", length(seq)))
  }
  #print(cutoff)
  good <- (cutoff>31 & cutoff<150)
  #for(i in 1:length(good)) print(paste(i, " ", length(seqs[[i]]), " ", cutoff[i]))
  goodProt <- seqs[good]
  goodCut <- cutoff[good]
  learningPositive <- vector(mode="list", length=length(good))
  learningNegative <- vector(mode="list", length=length(good))
  for (i in 1:sum(good)){
    learningPositive[i] = list(goodProt[[i]][1:goodCut[i]])
    learningNegative[i] = list(goodProt[[i]][(1+goodCut[i]):length(goodProt[[i]])] )
    #if(sum(is.na(learningPositive[[i]]))>0){
    #  print(goodProt[[i]])
    #  print(goodCut[i])
    #}
  }
  #print(unlist(learningPositive))
  list(con(unlist(learningPositive),length(aagroups)), con(unlist(learningNegative),length(aagroups)))
}

degenerate <- function(seq, aa_group) {
  for (i in 1L:length(aa_group)) {
    seq[seq %in% aa_group[[i]]] <- i
  }
  seq 
}

